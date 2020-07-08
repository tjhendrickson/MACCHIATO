#!/opt/Miniconda3/bin/python

from __future__ import print_function
import argparse
import os
from subprocess import Popen, PIPE
import subprocess
from multiprocessing import Pool, Lock
from utils.find_bolds import create_bold_lists 
from utils.workflow_logger import workflow_log
from workflow.core import execute_MACCHIATO_instances

# function which actually launches processes to underlying system
def run(command, env={}, cwd=None):
    merged_env = os.environ
    merged_env.update(env)
    merged_env.pop("DEBUG", None)
    print(command)
    process = Popen(command, stdout=PIPE, stderr=subprocess.STDOUT,
                    shell=True, env=merged_env, cwd=cwd)
    while True:
        line = process.stdout.readline()
        print(line)
        line = str(line)[:-1]
        if line == '' and process.poll() != None:
            break
    if process.returncode != 0:
        raise Exception("Non zero return code: %d"%process.returncode)

            

# specify arguments that MACCHIATO accepts
parser = argparse.ArgumentParser(description='')
parser.add_argument('input_dir', help='The directory where the preprocessed derivative needed live')
parser.add_argument('output_dir', help='The directory where the output files should be stored.')
parser.add_argument('analysis_level', help='Level of the analysis that will be performed. ' 
                    'If "group" is specified the tool will run through all available participants independently and output to same output directory. '
                    'If "participant" is specified one must include arguments to "--subject_label" and --session_label".',
                    choices=['group','participant'])
parser.add_argument('--subject_label', help='The label of the participant that should be analyzed. The label '
                    'corresponds to sub-<subject_label> from the BIDS spec '
                    '(so it does not include "sub-"). If this parameter is not '
                    'provided all subjects should be analyzed. Multiple '
                    'participants can be specified with a space separated list.',
                    nargs="+")
parser.add_argument('--session_label', help='The label of the session that should be analyzed. The label '
                    'corresponds to ses-<session_label> from the BIDS spec '
                    '(so it does not include "ses-"). If this parameter is not '
                    'provided all sessions should be analyzed. Multiple '
                    'sessions can be specified with a space separated list.',
                    nargs="+")

parser.add_argument('--preprocessing_type', help='BIDS-apps preprocessing pipeline run on data. Choices include "HCP" and "fmriprep". ',choices=['HCP','fmriprep'],default='HCP')
parser.add_argument('--use_ICA_outputs',help='Use ICA (whether FIX or AROMA) outputs for network matrix estimation. Choices include "Y/yes" or "N/no".',choices=['Yes','yes','No','no'],default='Yes')
parser.add_argument('--combine_resting_scans',help='If multiple of the same resting state BIDS file type exist should they be combined? Choices include "Y/yes" or "N/no".',choices=['Yes','yes','No','no'],default='No')
parser.add_argument('--parcellation_file', help='The CIFTI label file to use or used to parcellate the brain. ', default='NONE')
parser.add_argument('--parcellation_name', help='Shorthand name of the CIFTI label file. ', default='NONE')
# parser.add_argument('--motion_confounds',help='What type of motion confounds to use, if any. Choices are "Movement_Regressors" (motion rotation angles and translations in mm), '
#                                         ' "Movement_Regressors_dt" (detrended motion rotation angles and translations in mm), "Movement_Regressors_demean" (demeaned motion rotation angles and translations in mm) "Movement_RelativeRMS" (RMS intensity difference of volume N to the reference volume), '
#                                         ' "Movement_RelativeRMS_mean" (square of RMS intensity difference of volume N to the reference volume), "Movement_AbsoluteRMS" (absolute RMS intensity difference of volume N to the reference volume, '
#                                         ' "Movement_AbsoluteRMS_mean" (square of absolute RMS intensity difference of volume N to the reference volume), "dvars" ( RMS intensity difference of volume N to volume N+1 (see Power et al, NeuroImage, 59(3), 2012)), '
#                                         ' "fd" ( frame displacement (average of rotation and translation parameter differences - using weighted scaling, as in Power et al.))',
#                                         choices = ['NONE','Movement_Regressors','Movement_Regressors_dt','Movement_RelativeRMS','Movement_RelativeRMS_mean','Movement_AbsoluteRMS','Movement_AbsoluteRMS_mean','dvars','fd'],default='NONE')
parser.add_argument('--reg_name',help='What type of registration do you want to use? Choices are "MSMAll_2_d40_WRN" and "NONE"',choices = ['NONE','MSMAll_2_d40_WRN'],default='MSMAll_2_d40_WRN')
parser.add_argument('--apply_Fishers_r_to_z_transform', help="For correlation outputs, should Fisher's r-to-z transformation be applied? Choises are 'Yes' or 'No'.", choices = ['Yes','YES','yes','No','NO','no'],default='Yes')
parser.add_argument('--network_matrix_calculation', help="What method to employ for network matrix estimation. "
                                                    " Choices are 'All', 'correlation','partial_correlation', "
                                                    " 'dynamic_time_warping', 'tangent', 'covariance', 'sparse_inverse_covariance', "
                                                    "  'precision', 'sparse_inverse_precision'. NOTE: Specifying sparse matrices or dynamic time warping will result in increased computation time. For more information on available methods https://nilearn.github.io/connectivity/index.html#functional-connectivity-and-resting-state", 
                                                    choices =['All','all','correlation','partial_correlation',
                                                              'dynamic_time_warping','covariance',
                                                              'precision','sparse_inverse_precision',
                                                              'sparse_inverse_covariance'], default='correlation',nargs='+')
parser.add_argument('--graph_theory',help="Whether or not to output graph theoretical measures. Choices are 'All', 'clustering_coefficient','local_efficiency','strength','node_betweenness_centrality', 'edge_betweenness_centrality', 'eigenvector_centrality', and 'None'. If there are multiple measures that are of interest but not all are, separate as many choices as interested with a space.", choices = ['All','all','clustering_coefficient','local_efficiency','strength',
'node_betweenness_centrality', 'edge_betweenness_centrality', 'eigenvector_centrality','None','none','NONE'], default = 'NONE',nargs='+')
parser.add_argument('--wavelet',help="Whether or not to output wavelet matrices which includes measures of entropy, strength, and diversity. Which includes functional connectivity range scales: Scale 1 (0.23-0.45 Hz), Scale 2 (0.11-0.23 Hz), Scale 3 (0.06-0.11 Hz), Scale 4 (0.03-0.06 Hz), Scale 5 (0.01-0.03 Hz), Scale 6 (0.007-0.01 Hz). Choices are 'Yes' or 'No'.", 
                    choices = ['Yes','YES','yes','No','NO','no'],default='No')
parser.add_argument('--num_cpus', help='How many concurrent CPUs to use',default=1)
args = parser.parse_args()

# now parse arguments and print to standard output (STDOUT)
layout,ICA_outputs,graph_theory,network_matrix_calculation,combine_resting_scans,wavelet = workflow_log(preprocessing_type=args.preprocessing_type,
                                                                                                    #motion_confounds=args.motion_confounds,
                                                                                                    ICA_outputs=args.use_ICA_outputs,
                                                                                                    graph_theory=args.graph_theory,
                                                                                                    network_matrix_calculation=args.network_matrix_calculation,
                                                                                                    input_dir=args.input_dir,
                                                                                                    parcel_name=args.parcel_name,
                                                                                                    parce_file=args.parcel_file,
                                                                                                    fishers_r_to_z_transform=args.apply_Fishers_r_to_z_transform,
                                                                                                    selected_reg_name=args.reg_name,
                                                                                                    wavelet=args.wavelet,
                                                                                                    combine_resting_scans=args.combine_resting_scans)

# set up multiprocessing/parallelization allocation
l = Lock()
multiproc_pool = Pool(int(args.num_cpus))

# create list of bold scans to be run
bolds = create_bold_lists(layout=layout,
                          combine_resting_scans=args.combine_resting_scans,
                          preprocessing_type=args.preprocessing_type,
                          ICA_outputs=ICA_outputs,
                          msm_all_reg_name="MSMAll_2_d40_WRN",
                          selected_reg_name=args.reg_name)

execute_MACCHIATO_instances(preprocessing_type=args.preprocessing_type,
            parcel_file=args.parcel_file,
            parcel_name=args.parcel_name,
            selected_reg_name=args.reg_name,
            #motion_confounds=args.motion_confounds,
            ICA_outputs=ICA_outputs,
            combine_resting_scans=combine_resting_scans,
            wavelet=wavelet,
            network_matrix_calculation=args.network_matrix_calculation,
            output_dir=args.output_dir, 
            bold=bolds[0])

# multiproc_pool.map(partial(execute_MACCHIATO_instances,preprocessing_type=args.preprocessing_type,
#             parcel_file=args.parcel_file,
#             parcel_name=args.parcel_name,
#             selected_reg_name=args.reg_name,
#             motion_confounds_filename=motion_confounds_filename,
#             ICA_outputs=ICA_outputs,
#             combine_resting_scans=combine_resting_scans,
#             wavelet=wavelet,
#             network_matrix_calculation=network_matrix_calculation,
#             output_dir=args.output_dir,
#             msm_all_reg_name = "MSMAll_2_d40_WRN"),
#             sorted(bolds))
