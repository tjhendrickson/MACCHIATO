#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  7 16:16:34 2020

@author: timothy
"""

from bids.grabbids import BIDSLayout
import os

# ensures arguments specified to MACCHIATO are parsed correctly
# then prints specified arguments to standard output (STDOUT)
def workflow_log(preprocessing_type,motion_confounds,ICA_outputs,
                 graph_theory,network_matrix_calculation,input_dir,parcel_name,
                 parcel_file,fishers_r_to_z_transform,selected_reg_name,
                 wavelet,combine_resting_scans):
    try:
        layout = BIDSLayout(os.path.join(input_dir))
    except: 
        raise TypeError('Input folder: {input_dir} does not look like a BIDS folder or {input_dir} is not a path to a folder. '.format(input_dir))
    # check on arguments
    if combine_resting_scans == 'Yes' or combine_resting_scans == 'yes':
        combine_resting_scans = 'YES'
    elif combine_resting_scans == 'No' or combine_resting_scans == 'no':
        combine_resting_scans = 'NO'
    if wavelet == 'Yes' or wavelet == 'yes':
        wavelet = 'YES'
    elif wavelet == 'No' or wavelet == 'no':
        wavelet = 'NO'
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
    elif preprocessing_type == 'fmriprep' and motion_confounds != 'NONE':
        pass
    
     # use ICA outputs
    if ICA_outputs == 'yes' or ICA_outputs == 'Yes':
        ICA_outputs = 'YES'
    else:
        ICA_outputs = 'NO'
    # if all graph theory metrics are requested, transform arg from string to list
    if graph_theory == 'all' or graph_theory == 'All':
        graph_theory = ['clustering_coefficient','local_efficiency','strength',
    'node_betweenness_centrality', 'edge_betweenness_centrality', 'eigenvector_centrality']
    if network_matrix_calculation == 'all' or network_matrix_calculation == 'All':
        network_matrix_calculation = ['correlation','partial_correlation','dynamic_time_warping',
                                           'tangent','covariance', 'precision',
                                           'sparse_inverse_precision','sparse_inverse_covariance']
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
    print('\t-Network matrix metric/s to compute: %s' %(network_matrix_calculation))
    print("\t-Whether or not to compute Fisher's r-to-z transform to network matrices: %s" %(fishers_r_to_z_transform))
    print('\t-Whether or not to perform wavelet entropy on matrices: %s' %(wavelet))
    print('\t-Graph theory metric/s to compute: %s' %(graph_theory))
    print('\t-Input registration file to be used: %s' %str(selected_reg_name))
    print('\t-Whether motion confounds will be used for output: %s' %str(motion_confounds))
    print('\t-The preprocessing pipeline that the input comes from: %s' %str(preprocessing_type))
    print('\t-Use ICA outputs: %s' %str(ICA_outputs))
    print('\t-Use mixed effects if multiple of same acquisition: %s' %str(combine_resting_scans))
    print('\n')
    
    return layout,motion_confounds_filename,ICA_outputs,graph_theory,network_matrix_calculation,combine_resting_scans,wavelet

