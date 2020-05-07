#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  7 17:22:55 2020

@author: timothy
"""

import os
from ..utils.motion_confounds import fmri_motion_confounds
import nibabel
from .connectivity_metrics import Network
import pdb


#handles all logic of inputted data to run MACCHIATO        
def execute_MACCHIATO_instances(bold,preprocessing_type,parcel_file,parcel_name,
                  selected_reg_name,motion_confounds_filename,ICA_outputs,
                  combine_resting_scans,wavelet,network_matrix_calculation,
                  output_dir,msm_all_reg_name):
    if combine_resting_scans == 'NO':
        fmritcs = bold
        level_2_foldername = 'NONE'
        if 'ses' in fmritcs:
            subject_label = fmritcs.split('sub-')[1].split('/')[0]
            ses_label = fmritcs.split('ses-')[1].split('/')[0]
            # set output folder path
            output_dir=output_dir + "/sub-%s/ses-%s" % (subject_label, ses_label)
        else:
            subject_label = fmritcs.split('sub-')[1].split('/')[0]
            output_dir=output_dir + "/sub-%s" % (subject_label)
        if not os.path.isfile(os.path.join(output_dir,
                                           os.path.basename(fmritcs).split('.')[0] +
                                           '_'+parcel_name+'.ptseries.nii')):
            if preprocessing_type == 'HCP':
                if ICA_outputs == 'YES':
                    if selected_reg_name == msm_all_reg_name:
                        vol_fmritcs=fmritcs.replace('_Atlas_MSMAll_2_d40_WRN_hp2000_clean.dtseries.nii','_hp2000_clean.nii.gz')
                    else:
                        vol_fmritcs=fmritcs.replace('_Atlas_hp2000_clean.dtseries.nii','_hp2000_clean.nii.gz')
                else:
                    if selected_reg_name == msm_all_reg_name:
                        vol_fmritcs = fmritcs.replace('_Atlas_MSMAll_2_d40_WRN_hp2000.dtseries.nii','_hp2000.nii.gz')
                    else:
                        vol_fmritcs = fmritcs.replace('_Atlas_hp2000.dtseries.nii','_hp2000.nii.gz')
                fmri_motion_confounds(fmritcs=fmritcs,vol_fmritcs=vol_fmritcs,
                                      motion_confounds_filename=motion_confounds_filename,
                                      ICA_outputs=ICA_outputs,selected_reg_name=selected_reg_name,
                                      msm_all_reg_name=msm_all_reg_name,
                                      output_dir=output_dir,
                                      subject_label=subject_label,
                                      ses_label=ses_label)
                fmriname = os.path.basename(fmritcs).split(".")[0]
                assert fmriname
            elif preprocessing_type == 'fmriprep':
                pass
            Network_init = Network(output_dir,fmritcs,parcel_file,parcel_name)
            pdb.set_trace()
            if len(network_matrix_calculation) == 1:
                Network_init.create_network_matrix(method=network_matrix_calculation[0])
            else:
                for network_matrix in network_matrix_calculation:
                    Network_init.create_network_matrix(method=network_matrix)
                    
        else:
            print('data already exists within: ' +output_dir)
    else:
        level_2_foldername = 'rsfMRI_combined'
        fmrinames = []
        fmritcs = bold[0]
        if 'ses' in fmritcs:
            subject_label = fmritcs.split('sub-')[1].split('/')[0]
            ses_label = fmritcs.split('ses-')[1].split('/')[0]
            # set output folder path
            output_dir=output_dir + "/sub-%s/ses-%s" % (subject_label, ses_label)
        else:
            subject_label = fmritcs.split('sub-')[1].split('/')[0]
            output_dir=output_dir + "/sub-%s" % (subject_label)
        if not os.path.isfile(os.path.join(output_dir,'rsfMRI_combined_hp200_s4_level2.fsf')):
            for fmritcs in bold:
                pass # call to rsfMRI_network_metrics or to internal single run function
                fmriname = os.path.basename(fmritcs).split(".")[0] 
                fmrinames.append(fmriname)
            #retrieve subject and session numbers
            if 'ses' in bold[0]:
                subject_label = fmritcs.split('sub-')[1].split('/')[0]
                ses_label = fmritcs.split('ses-')[1].split('/')[0]
                # set output folder path
                output_dir = output_dir + "/sub-%s/ses-%s" % (subject_label, ses_label)
            else:
                subject_label = fmritcs.split('sub-')[1].split('/')[0]
                output_dir = output_dir + "/sub-%s" % (subject_label)
            if preprocessing_type == 'HCP':
                if ICA_outputs == 'YES':
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
