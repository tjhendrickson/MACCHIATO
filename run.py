#!/opt/Miniconda3/bin/python

from __future__ import print_function
import argparse
import os
from functools import partial
from subprocess import Popen, PIPE
import subprocess
from multiprocessing import Pool, Lock
from workflow.utils import MACCHIATO_setup
from tools.core import execute_MACCHIATO_instances
import pdb

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
parser.add_argument('--input_dir', help='The directory where the preprocessed derivative needed live',required=True)
parser.add_argument('--output_dir', help='The directory where the output files should be stored.',required=True)
parser.add_argument('--group', help='Whether to run this participant by participant or the entire group. Choices are "participant" or "batch". If participant by participant "--participant_label" and "--session_label" must be specified',choices = ['participant','batch'],required=True)
parser.add_argument('--participant_label', help='The label of the participant that should be analyzed. The label '
                   'corresponds to sub-<participant_label> from the BIDS spec '
                   '(so it does not include "sub-"). If this parameter is not '
                   'provided all subjects should be analyzed. Multiple '
                   'participants can be specified with a space separated list.',
                   nargs="+")
parser.add_argument('--session_label', help='The label of the session that should be analyzed. The label '
                   'corresponds to ses-<session_label> from the BIDS spec '
                   '(so it does not include "ses-"). If this parameter is not '
                   'provided all sessions within a subject should be analyzed.',
                   nargs="+")
parser.add_argument('--preprocessing_type', help='BIDS-apps preprocessing pipeline run on data. Choices include "HCP" and "fmriprep". ',choices=['HCP','fmriprep'],default='HCP')
parser.add_argument('--use_ICA_outputs',help='Use ICA (whether FIX or AROMA) outputs for network matrix estimation. Choices include "Y/yes" or "N/no".',choices=['Yes','yes','No','no'],default='Yes')
parser.add_argument('--combine_resting_scans',help='If multiple of the same resting state BIDS file type exist should they be combined? Choices include "Y/yes" or "N/no".',choices=['Yes','yes','No','no'],default='No')
parser.add_argument('--parcellation_file', help='The CIFTI label file to use or used to parcellate the brain. ', required=True)
parser.add_argument('--parcellation_name', help='Shorthand name of the CIFTI label file. ', required=True)
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


# now parse arguments and print to standard output (STDOUT), slight difference in call if done as batch or participant

if args.group == 'batch':
    setupParams = MACCHIATO_setup(group=args.group,
                    preprocessing_type=args.preprocessing_type,
                    ICA_outputs=args.use_ICA_outputs,
                    graph_theory=args.graph_theory,
                    network_matrix_calculation=args.network_matrix_calculation,
                    input_dir=args.input_dir,
                    parcel_name=args.parcellation_name,
                    parcel_file=args.parcellation_file,
                    fishers_r_to_z_transform=args.apply_Fishers_r_to_z_transform,
                    selected_reg_name=args.reg_name,
                    wavelet=args.wavelet,
                    combine_resting_scans=args.combine_resting_scans)
elif args.group == 'participant':
    setupParams = MACCHIATO_setup(group=args.group,
                    preprocessing_type=args.preprocessing_type,
                    ICA_outputs=args.use_ICA_outputs,
                    graph_theory=args.graph_theory,
                    network_matrix_calculation=args.network_matrix_calculation,
                    input_dir=args.input_dir,
                    parcel_name=args.parcellation_name,
                    fishers_r_to_z_transform=args.apply_Fishers_r_to_z_transform,
                    parcel_file=args.parcellation_file,
                    selected_reg_name=args.reg_name,
                    wavelet=args.wavelet,
                    combine_resting_scans=args.combine_resting_scans,
                    participant_label=args.participant_label,
                    session_label=args.session_label)

#transform "self" variables from MACCHIATO_setup to dictionary
setupParams = setupParams.__dict__
# set up multiprocessing/parallelization allocation
l = Lock()
multiproc_pool = Pool(int(args.num_cpus))

execute_MACCHIATO_instances(preprocessing_type=args.preprocessing_type,
            parcel_file=args.parcellation_file,
            parcel_name=args.parcellation_name,
            selected_reg_name=args.reg_name,
            ICA_outputs=setupParams['ICA_outputs'],
            combine_resting_scans=setupParams['combine_resting_scans'],
            wavelet=setupParams['wavelet'],
            network_metric=setupParams['network_matrix_calculation'],
            output_dir=args.output_dir, 
            graph_theory_metric=setupParams['graph_theory'],
            bold=setupParams['bolds'][0])

# multiproc_pool.map(partial(execute_MACCHIATO_instances,
#                             preprocessing_type=args.preprocessing_type,
#                             parcel_file=args.parcellation_file,
#                             parcel_name=args.parcellation_name,
#                             selected_reg_name=args.reg_name,
#                             ICA_outputs=setupParams['ICA_outputs'],
#                             combine_resting_scans=setupParams['combine_resting_scans'],
#                             wavelet=setupParams['wavelet'],
#                             network_matrix_calculation=setupParams['network_matrix_calculation'],
#                             output_dir=args.output_dir, 
#                             graph_theory=setupParams['graph_theory']),
#                             sorted(setupParams['bolds']))
