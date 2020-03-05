#!/opt/Miniconda3/bin/python

from __future__ import print_function
import argparse
import os
import nibabel
import nibabel.cifti2 as ci
from glob import glob
from subprocess import Popen, PIPE
import subprocess
from bids.grabbids import BIDSLayout
from functools import partial
from collections import OrderedDict
import pdb
from rsfMRI_network_metrics import NetworkIO
from multiprocessing import Pool, Lock
import time

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

parser = argparse.ArgumentParser(description='')
parser.add_argument('input_dir', help='The directory where the preprocessed derivative needed live')
parser.add_argument('output_dir', help='The directory where the output files should be stored.')
parser.add_argument('--preprocessing_type', help='BIDS-apps preprocessing pipeline run on data. Choices include "HCP" and "fmriprep". ',choices=['HCP','fmriprep'],default='HCP')
parser.add_argument('--use_ICA_outputs',help='Use ICA (whether FIX or AROMA) outputs for network matrix estimation. Choices include "Y/yes" or "N/no".',choices=['Yes','yes','No','no'],default='Yes')
parser.add_argument('--combine_resting_scans',help='If multiple of the same resting state BIDS file type exist should they be combined? Choices include "Y/yes" or "N/no".',choices=['Yes','yes','No','no'],default='No')
parser.add_argument('--parcellation_file', help='The CIFTI label file to use or used to parcellate the brain. ', default='NONE')
parser.add_argument('--parcellation_name', help='Shorthand name of the CIFTI label file. ', default='NONE')
parser.add_argument('--motion_confounds',help='What type of motion confounds to use, if any. Choices are "Movement_Regressors" (motion rotation angles and translations in mm), '
                                        ' "Movement_Regressors_dt" (detrended motion rotation angles and translations in mm), "Movement_Regressors_demean" (demeaned motion rotation angles and translations in mm) "Movement_RelativeRMS" (RMS intensity difference of volume N to the reference volume), '
                                        ' "Movement_RelativeRMS_mean" (square of RMS intensity difference of volume N to the reference volume), "Movement_AbsoluteRMS" (absolute RMS intensity difference of volume N to the reference volume, '
                                        ' "Movement_AbsoluteRMS_mean" (square of absolute RMS intensity difference of volume N to the reference volume), "dvars" ( RMS intensity difference of volume N to volume N+1 (see Power et al, NeuroImage, 59(3), 2012)), '
                                        ' "fd" ( frame displacement (average of rotation and translation parameter differences - using weighted scaling, as in Power et al.))',
                                        choices = ['NONE','Movement_Regressors','Movement_Regressors_dt','Movement_RelativeRMS','Movement_RelativeRMS_mean','Movement_AbsoluteRMS','Movement_AbsoluteRMS_mean','dvars','fd'],default='NONE')
parser.add_argument('--reg_name',help='What type of registration do you want to use? Choices are "MSMAll_2_d40_WRN" and "NONE"',choices = ['NONE','MSMAll_2_d40_WRN'],default='MSMAll_2_d40_WRN')
parser.add_argument('--apply_Fishers_r_to_z_transform', help="For correlation outputs, should Fisher's r-to-z transformation be applied? Choises are 'Yes' or 'No'.", choices = ['Yes','YES','yes','No','NO','no'],default='Yes')
parser.add_argument('--network_matrix_calculation', help="What method to employ for network matrix estimation. "
                                                    " Choices are 'All', 'correlation','partial_correlation', "
                                                    " 'dynamic_time_warping', 'tangent', 'covariance', 'sparse_inverse_covariance', "
                                                    "  'precision', 'sparse_inverse_precision'. ", 
                                                    choices =['All','all','correlation','partial_correlation',
                                                              'dynamic_time_warping','tangent','covariance',
                                                              'precision','sparse_inverse_precision',
                                                              'sparse_inverse_covariance'], default='correlation',nargs='+')
parser.add_argument('--graph_theory',help="Whether or not to output graph theoretical measures. Choices are 'All', 'clustering_coefficient','local_efficiency','strength','node_betweenness_centrality', 'edge_betweenness_centrality', 'eigenvector_centrality', and 'None'. If there are multiple measures that are of interest but not all are, separate as many choices as interested with a space.", choices = ['All','all','clustering_coefficient','local_efficiency','strength',
'node_betweenness_centrality', 'edge_betweenness_centrality', 'eigenvector_centrality','None','none','NONE'], default = 'NONE',nargs='+')
parser.add_argument('--wavelet',help="Whether or not to output wavelet matrices which includes measures of entropy, strength, and diversity. Which includes functional connectivity range scales: Scale 1 (0.23-0.45 Hz), Scale 2 (0.11-0.23 Hz), Scale 3 (0.06-0.11 Hz), Scale 4 (0.03-0.06 Hz), Scale 5 (0.01-0.03 Hz), Scale 6 (0.007-0.01 Hz). Choices are 'Yes' or 'No'.", 
                    choices = ['Yes','YES','yes','No','NO','no'],default='No')
parser.add_argument('--num_cpus', help='How many concurrent CPUs to use',default=1)
args = parser.parse_args()

# global variables
output_dir = args.output_dir
parcel_file = args.parcellation_file
parcel_name = args.parcellation_name
selected_reg_name = args.reg_name
msm_all_reg_name = "MSMAll_2_d40_WRN"
preprocessing_type = args.preprocessing_type
motion_confounds = args.motion_confounds
combine_resting_scans = args.combine_resting_scans

if preprocessing_type == 'HCP':
    if not motion_confounds == 'NONE':
        motion_confounds_dict = {'Movement_Regressors': 'Movement_Regressors.txt',
        'Movement_Regressors_dt': 'Movement_Regressors_dt.txt',
        'Movement_Regressors_demean': 'Movement_Regressors_demean.txt',
        'Movement_RelativeRMS': 'Movement_RelativeRMS.txt',
        'Movement_RelativeRMS_mean': 'Movement_RelativeRMS_mean.txt',
        'Movement_AbsoluteRMS': 'Movement_AbsoluteRMS.txt',
        'Movement_AbsoluteRMS_mean': 'Movement_AbsoluteRMS_mean.txt',
        'dvars': 'Movement_dvars.txt',
        'fd': 'Movement_fd.txt'}
        motion_confounds_filename = motion_confounds_dict[motion_confounds]
    else:
        motion_confounds_filename = 'NONE'
        motion_confounds_filepath = 'NONE'
elif preprocessing_type == 'fmriprep' and motion_confounds != 'NONE':
    pass

 # use ICA outputs
if args.use_ICA_outputs == 'yes' or args.use_ICA_outputs == 'Yes':
    ICAoutputs = 'YES'
else:
    ICAoutputs = 'NO'
# if all graph theory metrics are requested, transform arg from string to list
if args.graph_theory == 'all' or args.graph_theory == 'All':
    args.graph_theory = ['clustering_coefficient','local_efficiency','strength',
'node_betweenness_centrality', 'edge_betweenness_centrality', 'eigenvector_centrality']
if args.network_matrix_calculation == 'all' or args.network_matrix_calculation == 'All':
    args.network_matrix_calculation = ['correlation','partial_correlation','dynamic_time_warping',
                                       'tangent','covariance', 'precision',
                                       'sparse_inverse_precision','sparse_inverse_covariance']
# check on arguments
layout = BIDSLayout(os.path.join(args.input_dir))
print("Running MACCHIATO ")
if parcel_file == 'NONE':
    print('\n')
    raise ValueError('Parcellating output selected but no parcel file specified after argument "--parcellation_file". Exiting.')
else:
    print('\t-Parcellation file to be used to parcellate outputs: %s' %str(parcel_file))
if parcel_name == 'NONE':
    print('\n')
    raise ValueError('Parcellating output selected but no parcel name specified after argument "--parcellation_name". Exiting.')
else:
    print('\t-Short hand parcellation name to be used: %s' %str(parcel_name))
print('\t-Network matrix metric/s to compute: %s' %(args.network_matrix_calculation))
print("\t-Whether or not to compute Fisher's r-to-z transform to network matrices: %s" %(args.apply_Fishers_r_to_z_transform))
print('\t-Whether or not to perform wavelet entropy on matrices: %s' %(args.wavelet))
print('\t-Graph theory metric/s to compute: %s' %(args.graph_theory))
print('\t-Input registration file to be used: %s' %str(selected_reg_name))
print('\t-Whether motion confounds will be used for output: %s' %str(motion_confounds))
print('\t-The preprocessing pipeline that the input comes from: %s' %str(preprocessing_type))
print('\t-Use ICA outputs: %s' %str(ICAoutputs))
print('\t-Use mixed effects if multiple of same acquisition: %s' %str(args.combine_resting_scans))
print('\n')

l = Lock()
multiproc_pool = Pool(int(args.num_cpus))

def run_MACCHIATO(bold,ICAstring,preprocessing_type,parcel_file,parcel_name,
                  selected_reg_name,motion_confounds,ICAoutputs,network_matrix_calculation,
                  combine_resting_scans,output_dir):
    if combine_resting_scans == 'No' or combine_resting_scans == 'no':
        fmritcs = bold
        level_2_foldername = 'NONE'
        if 'ses' in fmritcs:
            subject_label = fmritcs.split('sub-')[1].split('/')[0]
            ses_label = fmritcs.split('ses-')[1].split('/')[0]
            # set output folder path
            outdir=output_dir + "/sub-%s/ses-%s" % (subject_label, ses_label)
        else:
            subject_label = fmritcs.split('sub-')[1].split('/')[0]
            outdir=output_dir + "/sub-%s" % (subject_label)
        if not os.path.isfile(os.path.join(outdir,
                                           os.path.basename(fmritcs).split('.')[0] +
                                           '_'+parcel_name+'.ptseries.nii')):
            if preprocessing_type == 'HCP':
                if ICAoutputs == 'YES':
                    if selected_reg_name == msm_all_reg_name:
                        vol_fmritcs=fmritcs.replace('_Atlas_MSMAll_2_d40_WRN_hp2000_clean.dtseries.nii','_hp2000_clean.nii.gz')
                    else:
                        vol_fmritcs=fmritcs.replace('_Atlas_hp2000_clean.dtseries.nii','_hp2000_clean.nii.gz')
                else:
                    if selected_reg_name == msm_all_reg_name:
                        vol_fmritcs = fmritcs.replace('_Atlas_MSMAll_2_d40_WRN_hp2000.dtseries.nii','_hp2000.nii.gz')
                    else:
                        vol_fmritcs = fmritcs.replace('_Atlas_hp2000.dtseries.nii','_hp2000.nii.gz')
                shortfmriname=fmritcs.split("/")[-2]
                # create confounds if dvars or fd selected
                if motion_confounds_filename == 'Movement_dvars.txt':
                    os.system("${FSL_DIR}/bin/fsl_motion_outliers -i " + vol_fmritcs + \
                                    " -o " + outdir + "/sub-" + subject_label + "/ses-" + \
                                        ses_label + "/MNINonLinear/" + "Results/" + shortfmriname + "/" + motion_confounds_filename + " --dvars")
                elif motion_confounds_filename == 'Movement_fd.txt':
                    os.system("${FSL_DIR}/bin/fsl_motion_outliers -i " + vol_fmritcs + \
                                    " -o " + outdir + "/sub-" + subject_label + "/ses-" + \
                                        ses_label + "/MNINonLinear/" + "Results/" + shortfmriname + "/" + motion_confounds_filename + " --fd")
                # create full path to confounds file if not 'NONE'
                if motion_confounds_filename != 'NONE' and ICAoutputs == 'YES':
                    if selected_reg_name == msm_all_reg_name:
                        motion_confounds_filepath = fmritcs.replace(shortfmriname+'_Atlas_MSMAll_2_d40_WRN_hp2000_clean.dtseries.nii',motion_confounds_filename)
                    else:
                        motion_confounds_filepath = fmritcs.replace(shortfmriname+'_Atlas_hp2000_clean.dtseries.nii',motion_confounds_filename)
                elif motion_confounds_filename != 'NONE' and ICAoutputs == 'NO':
                    if selected_reg_name == msm_all_reg_name:
                        motion_confounds_filepath = fmritcs.replace(shortfmriname+'_Atlas_MSMAll_2_d40_WRN_hp2000.dtseries.nii',motion_confounds_filename)
                    else:
                        motion_confounds_filepath = fmritcs.replace(shortfmriname+'_Atlas_hp2000.dtseries.nii',motion_confounds_filename)
                fmriname = os.path.basename(fmritcs).split(".")[0]
                assert fmriname
            elif preprocessing_type == 'fmriprep':
                pass
            pdb.set_trace()
            Network_init = NetworkIO(output_dir,fmritcs,parcel_file,parcel_name)
            Network_init.create_network_matrix(method=network_matrix_calculation)
        else:
            print('data already exists within: ' +outdir)
    else:
        level_2_foldername = 'rsfMRI_combined'
        fmrinames = []
        fmritcs = bold[0]
        if 'ses' in fmritcs:
            subject_label = fmritcs.split('sub-')[1].split('/')[0]
            ses_label = fmritcs.split('ses-')[1].split('/')[0]
            # set output folder path
            outdir=output_dir + "/sub-%s/ses-%s" % (subject_label, ses_label)
        else:
            subject_label = fmritcs.split('sub-')[1].split('/')[0]
            outdir=output_dir + "/sub-%s" % (subject_label)
        if not os.path.isfile(os.path.join(outdir,'rsfMRI_combined_hp200_s4_level2.fsf')):
            for fmritcs in bold:
                pass # call to rsfMRI_network_metrics or to internal single run function
                fmriname = os.path.basename(fmritcs).split(".")[0] 
                fmrinames.append(fmriname)
            #retrieve subject and session numbers
            if 'ses' in bold[0]:
                subject_label = fmritcs.split('sub-')[1].split('/')[0]
                ses_label = fmritcs.split('ses-')[1].split('/')[0]
                # set output folder path
                outdir = output_dir + "/sub-%s/ses-%s" % (subject_label, ses_label)
            else:
                subject_label = fmritcs.split('sub-')[1].split('/')[0]
                outdir = output_dir + "/sub-%s" % (subject_label)
            if preprocessing_type == 'HCP':
                if ICAoutputs == 'YES':
                    if selected_reg_name == msm_all_reg_name:
                        vol_fmritcs = bold[0].replace('_Atlas_MSMAll_2_d40_WRN_hp2000_clean.dtseries.nii','_hp2000_clean.nii.gz')
                    else:
                        vol_fmritcs = bold[0].replace('_Atlas_hp2000_clean.dtseries.nii','_hp2000_clean.nii.gz')
                else:
                    if selected_reg_name == msm_all_reg_name:
                        vol_fmritcs = bold[0].replace('_Atlas_MSMAll_2_d40_WRN_hp2000.dtseries.nii','_hp2000.nii.gz')
                    else:
                        vol_fmritcs = bold[0].replace('_Atlas_hp2000.dtseries.nii','_hp2000.nii.gz')
                zooms = nibabel.load(vol_fmritcs).get_header().get_zooms()
                fmrires = str(int(min(zooms[:3])))
                pass # call to rsfMRI_network_metrics    
if args.combine_resting_scans == 'No' or args.combine_resting_scans == 'no':
    if preprocessing_type == 'HCP':
        # use ICA outputs
        if ICAoutputs == 'YES':
            ICAstring="_FIXclean"
            if selected_reg_name == msm_all_reg_name:
                bolds = [f.filename for f in layout.get(type='clean',extensions="dtseries.nii",task='rest') if msm_all_reg_name+'_hp2000_clean' in f.filename]
            else:
                bolds = [f.filename for f in layout.get(type='clean',extensions="dtseries.nii", task='rest') if '_hp2000_clean' and not msm_all_reg_name in f.filename]
        # do not use ICA outputs
        else:
            ICAstring=""
            if selected_reg_name == msm_all_reg_name:
                bolds = [f.filename for f in layout.get(extensions="dtseries.nii", task='rest') if msm_all_reg_name + '_hp2000' in f.filename and not 'clean' in f.filename]
            else:
                bolds = [f.filename for f in layout.get(extensions="dtseries.nii", task='rest') if '_hp2000' in f.filename and not 'clean' and not msm_all_reg_name in f.filename]
    elif preprocessing_type == 'fmriprep':
        #use ICA outputs
        if ICAoutputs == 'YES':
            ICAstring="_AROMAclean"
            bolds = [f.filename for f in layout.get(type='bold',task='rest') if 'smoothAROMAnonaggr' in f.filename]
        # do not use ICA outputs
        else:
            ICAstring=""
            bolds = [f.filename for f in layout.get(type='bold',task='rest') if 'preproc' in f.filename]
        bolds_ref = [f.filename for f in layout.get(type='boldref',task='rest')]
    pdb.set_trace()
    run_MACCHIATO(ICAstring=ICAstring, 
                               preprocessing_type=preprocessing_type,
                               parcel_file=parcel_file,
                               parcel_name=parcel_name,
                               selected_reg_name=selected_reg_name,
                               motion_confounds=motion_confounds,
                               ICAoutputs=ICAoutputs,
                               combine_resting_scans=combine_resting_scans,
                               network_matrix_calculation=args.network_matrix_calculation,
                               output_dir=output_dir, bold=bolds[0])
    """
    multiproc_pool.map(partial(run_MACCHIATO,ICAstring=ICAstring, 
                               preprocessing_type=preprocessing_type,
                               parcel_file=parcel_file,
                               parcel_name=parcel_name,
                               selected_reg_name=selected_reg_name,
                               motion_confounds=motion_confounds,
                               ICAoutputs=ICAoutputs,
                               combine_resting_scans=combine_resting_scans,
                               output_dir=output_dir),
                sorted(bolds))
    """
else:
    combined_bolds_list = []
    if layout.get_sessions() > 0:
        for scanning_session in layout.get_sessions():
            if preprocessing_type == 'HCP':
                # use ICA outputs
                if ICAoutputs == 'YES':
                    ICAstring="_FIXclean"
                    if selected_reg_name == msm_all_reg_name:
                        bolds = [f.filename for f in layout.get(type='clean',extensions="dtseries.nii",task='rest',session=scanning_session) if msm_all_reg_name+'_hp2000_clean' in f.filename]
                    else:
                        bolds = [f.filename for f in layout.get(type='clean',extensions="dtseries.nii", task='rest',session=scanning_session) if '_hp2000_clean' and not msm_all_reg_name in f.filename]
                # do not use ICA outputs
                else:
                    ICAstring=""
                    if selected_reg_name == msm_all_reg_name:
                        bolds = [f.filename for f in layout.get(extensions="dtseries.nii", task='rest',session=scanning_session) if msm_all_reg_name + '_hp2000' in f.filename and not 'clean' in f.filename]
                    else:
                        bolds = [f.filename for f in layout.get(extensions="dtseries.nii", task='rest',session=scanning_session) if '_hp2000' in f.filename and not 'clean' and not msm_all_reg_name in f.filename]
            elif preprocessing_type == 'fmriprep':
                #use ICA outputs
                if ICAoutputs == 'YES':
                    ICAstring="_AROMAclean"
                    bolds = [f.filename for f in layout.get(type='bold',task='rest',session=scanning_session) if 'smoothAROMAnonaggr' in f.filename]
                # do not use ICA outputs
                else:
                    ICAstring=""
                    bolds = [f.filename for f in layout.get(type='bold',task='rest') if 'preproc' in f.filename]
                bolds_ref = [f.filename for f in layout.get(type='boldref',task='rest')]
            if len(bolds) >= 2:
                combined_bolds_list.append(bolds)
    else:
        for scanning_session in layout.get_subjects():
            if preprocessing_type == 'HCP':
                # use ICA outputs
                if ICAoutputs == 'YES':
                    ICAstring="_FIXclean"
                    if selected_reg_name == msm_all_reg_name:
                        bolds = [f.filename for f in layout.get(type='clean',extensions="dtseries.nii",task='rest',subject=scanning_session) if msm_all_reg_name+'_hp2000_clean' in f.filename]
                    else:
                        bolds = [f.filename for f in layout.get(type='clean',extensions="dtseries.nii", task='rest',subject=scanning_session) if '_hp2000_clean' and not msm_all_reg_name in f.filename]
                # do not use ICA outputs
                else:
                    ICAstring=""
                    if selected_reg_name == msm_all_reg_name:
                        bolds = [f.filename for f in layout.get(extensions="dtseries.nii", task='rest',subject=scanning_session) if msm_all_reg_name + '_hp2000' in f.filename and not 'clean' in f.filename]
                    else:
                        bolds = [f.filename for f in layout.get(extensions="dtseries.nii", task='rest',subject=scanning_session) if '_hp2000' in f.filename and not 'clean' and not msm_all_reg_name in f.filename]
            elif preprocessing_type == 'fmriprep':
                #use ICA outputs
                if ICAoutputs == 'YES':
                    ICAstring="_AROMAclean"
                    bolds = [f.filename for f in layout.get(type='bold',task='rest',subject=scanning_session) if 'smoothAROMAnonaggr' in f.filename]
                # do not use ICA outputs
                else:
                    ICAstring=""
                    bolds = [f.filename for f in layout.get(type='bold',task='rest') if 'preproc' in f.filename]
                bolds_ref = [f.filename for f in layout.get(type='boldref',task='rest')]
            if len(bolds) >= 2:
                combined_bolds_list.append(bolds)
    multiproc_pool.map(partial(run_MACCHIATO,ICAstring=ICAstring, 
                               preprocessing_type=preprocessing_type,
                               parcel_file=parcel_file,
                               parcel_name=parcel_name,
                               selected_reg_name=selected_reg_name,
                               motion_confounds=motion_confounds,
                               ICAoutputs=ICAoutputs,
                               combine_resting_scans=combine_resting_scans,
                               network_matrix_calculation=args.network_matrix_calculation,
                               output_dir=output_dir),
                sorted(combined_bolds_list))
    