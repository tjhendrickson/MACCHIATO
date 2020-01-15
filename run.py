#!/opt/Anaconda2/bin/python

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
from sklearn import cross_decomposition, covariance
import pandas as pd

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
parser.add_argument('--use_ICA_outputs',help='Use ICA (whether FIX or AROMA) outputs in seed analysis. Choices include "Y/yes" or "N/no".',choices=['Yes','yes','No','no'],default='Yes')
parser.add_argument('--stages',help='Which stages to run. Space separated list. ',nargs="+", choices=['rsfMRISeedAnalysis', 'Generatefsf'],default=['Generatefsf', 'rsfMRISeedAnalysis'])
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
parser.add_argument('--network_matrix_estimation', help="What method to employ for network matrix estimation. "
                                                    " Choices are 'correlation' = correlation,'partial-correlation' = partial correlation,'ridge-partial-correlation' = regularized partial correlation, 
                                                    " 'dynamic-time-warping' = dynamic time warping, 'tanget' = tanget, 'covariance' = covariance', 'sparse-inverse-covariance' = sparse inverse covariance, 
                                                    " 'group-sparse-covariance' = group sparse covariance, 'precision' = precision, 'coherence' = coherence, 'mutual-information' = mutual information, "
                                                    " 'canonical-correlation' = canonical correlation. ")
parser.add_argument('--graph_theory',help="Where or not to output graph theoretical measures which include: efficiency, centrality, cost/degree, strength, density, ")
args = parser.parse_args()

# global variables
highpass = "2000"
lowresmesh = 32
highresmesh = 164
smoothing = args.smoothing
parcel_file = args.parcellation_file
parcel_name = args.parcellation_name
selected_reg_name = args.reg_name
msm_all_reg_name = "MSMAll_2_d40_WRN"
preprocessing_type = args.preprocessing_type
motion_confounds = args.motion_confounds

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

# need a subject label in order to start
if args.participant_label:
    subject_label=args.participant_label[0]
    layout = BIDSLayout(os.path.join(args.input_dir,'sub-'+subject_label))
else:
    raise ValueError('An argument must be specified for participant label. Quitting.')
# if subject label has sessions underneath those need to be outputted into different directories
if args.session_label:
    ses_to_analyze = args.session_label
else:
    ses_to_analyze = layout.get_sessions(subject=subject_label)

if ses_to_analyze:
    for ses_label in ses_to_analyze:
        # set output folder path
        outdir=args.output_dir + "/sub-%s/ses-%s" % (subject_label, ses_label)
        if preprocessing_type == 'HCP':
            # use ICA outputs
            if ICAoutputs == 'YES':
                if selected_reg_name == msm_all_reg_name:
                    bolds = [f.filename for f in layout.get(subject=subject_label, session=ses_label,type='clean',
                                                            extensions="dtseries.nii", task='rest',) if msm_all_reg_name+'_hp2000_clean' in f.filename]
                else:
                    bolds = [f.filename for f in layout.get(subject=subject_label, session=ses_label,type='clean',
                                                            extensions="dtseries.nii", task='rest',) if '_hp2000_clean' and not msm_all_reg_name in f.filename]

            # do not use ICA outputs
            else:
                if selected_reg_name == msm_all_reg_name:
                    bolds = [f.filename for f in layout.get(subject=subject_label, session=ses_label,
                                                            extensions="dtseries.nii", task='rest') if msm_all_reg_name + '_hp2000' in f.filename and not 'clean' in f.filename]
                else:
                    bolds = [f.filename for f in layout.get(subject=subject_label, session=ses_label,
                                                            extensions="dtseries.nii", task='rest') if '_hp2000' in f.filename and not 'clean' and not msm_all_reg_name in f.filename]

        elif preprocessing_type == 'fmriprep':
            #use ICA outputs
            if ICAoutputs == 'YES':
                bolds = [f.filename for f in layout.get(subject=subject_label,session=ses_label,type='bold',task='rest') if 'smoothAROMAnonaggr' in f.filename]
            # do not use ICA outputs
            else:
                bolds = [f.filename for f in layout.get(subject=subject_label,session=ses_label,type='bold',task='rest') if 'preproc' in f.filename]
            # will need bold reference images
            bolds_ref = [f.filename for f in layout.get(subject=subject_label,session=ses_label,type='boldref',task='rest')]
        for idx,fmritcs in enumerate(bolds):
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
                zooms = nibabel.load(vol_fmritcs).get_header().get_zooms()
                fmrires = str(int(min(zooms[:3])))
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
                bold_ref = "NONE"
            elif preprocessing_type == 'fmriprep':
                #reference image
                bold_ref = bolds_ref[idx]
                vol_fmritcs='NONE'
else:
    outdir=args.output_dir + "/sub-%s" % (subject_label)
    # retrieve preprocessing BIDS layout for participant specified
    if preprocessing_type == 'HCP':
        # use ICA outputs
        if ICAoutputs == 'YES':
            if selected_reg_name == msm_all_reg_name:
                bolds = [f.filename for f in layout.get(subject=subject_label, type='clean',
                                                        extensions="dtseries.nii", task='rest',) if msm_all_reg_name+'_hp2000_clean' in f.filename]
            else:
                bolds = [f.filename for f in layout.get(subject=subject_label, type='clean',
                                                        extensions="dtseries.nii", task='rest',) if '_hp2000_clean' and not msm_all_reg_name in f.filename]
        # do not use ICA outputs
        else:
            if selected_reg_name == msm_all_reg_name:
                bolds = [f.filename for f in layout.get(subject=subject_label,
                                                        extensions="dtseries.nii", task='rest') if msm_all_reg_name + '_hp2000' in f.filename and not 'clean' in f.filename]
            else:
                bolds = [f.filename for f in layout.get(subject=subject_label,
                                                        extensions="dtseries.nii", task='rest') if '_hp2000' in f.filename and not 'clean' and not msm_all_reg_name in f.filename]

    elif preprocessing_type == 'fmriprep':
        #use ICA outputs
        if ICAoutputs == 'YES':
            bolds = [f.filename for f in layout.get(subject=subject_label,type='bold',task='rest') if 'smoothAROMAnonaggr' in f.filename]
        # do not use ICA outputs
        else:
            bolds = [f.filename for f in layout.get(subject=subject_label,type='bold',task='rest') if 'preproc' in f.filename]
        bolds_ref = [f.filename for f in layout.get(subject=subject_label,session=ses_label,type='boldref',task='rest')]
    for idx,fmritcs in enumerate(bolds):
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

            zooms = nibabel.load(vol_fmritcs).get_header().get_zooms()
            fmrires = str(int(min(zooms[:3])))
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

            fmriname = fmritcs.path.basename.split(".")[0]
            assert fmriname
            bold_ref = "NONE"
        elif preprocessing_type == 'fmriprep':
            bold_ref = bolds_ref[idx]
            vol_fmritcs="NONE"
    