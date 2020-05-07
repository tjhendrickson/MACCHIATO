#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  7 15:57:19 2020

@author: timothy
"""

import os

def fmri_motion_confounds(fmritcs,vol_fmritcs,motion_confounds_filename,
                          ICA_outputs,selected_reg_name,msm_all_reg_name,
                          output_dir,subject_label,ses_label):
    shortfmriname=fmritcs.split("/")[-2]
    # create confounds if dvars or fd selected
    if motion_confounds_filename == 'Movement_dvars.txt':
        os.system("${FSL_DIR}/bin/fsl_motion_outliers -i " + vol_fmritcs + \
                        " -o " + output_dir + "/sub-" + subject_label + "/ses-" + \
                            ses_label + "/MNINonLinear/" + "Results/" + shortfmriname + "/" + motion_confounds_filename + " --dvars")
    elif motion_confounds_filename == 'Movement_fd.txt':
        os.system("${FSL_DIR}/bin/fsl_motion_outliers -i " + vol_fmritcs + \
                        " -o " + output_dir + "/sub-" + subject_label + "/ses-" + \
                            ses_label + "/MNINonLinear/" + "Results/" + shortfmriname + "/" + motion_confounds_filename + " --fd")
    # create full path to confounds file if not 'NONE'
    if motion_confounds_filename != 'NONE' and ICA_outputs == 'YES':
        if selected_reg_name == msm_all_reg_name:
            motion_confounds_filepath = fmritcs.replace(shortfmriname+'_Atlas_MSMAll_2_d40_WRN_hp2000_clean.dtseries.nii',motion_confounds_filename)
        else:
            motion_confounds_filepath = fmritcs.replace(shortfmriname+'_Atlas_hp2000_clean.dtseries.nii',motion_confounds_filename)
    elif motion_confounds_filename != 'NONE' and ICA_outputs == 'NO':
        if selected_reg_name == msm_all_reg_name:
            motion_confounds_filepath = fmritcs.replace(shortfmriname+'_Atlas_MSMAll_2_d40_WRN_hp2000.dtseries.nii',motion_confounds_filename)
        else:
            motion_confounds_filepath = fmritcs.replace(shortfmriname+'_Atlas_hp2000.dtseries.nii',motion_confounds_filename)
    else:
        motion_confounds_filepath == 'NONE'
    return motion_confounds_filepath
