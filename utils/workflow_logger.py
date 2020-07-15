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
class workflow_log:
    def __init__(self,**kwargs):
        '''
        
        Parameters
        ----------
        group : [participant, batch]
            Whether to run this participant by participant or the entire group. 
        preprocessing_type : [HCP, fmriprep]
            BIDS-apps preprocessing pipeline run on data. Choices include "HCP" and "fmriprep". 
        ICA_outputs : [Yes, yes, No, no]
            Use ICA (whether FIX or AROMA) outputs for network matrix estimation.
        graph_theory : 
            Whether or not to output graph theoretical measures
        network_matrix_calculation : [All, all, correlation, partial_correlation, dynamic_time_warping, covariance, precision, sparse_inverse_precision, sparse_inverse_covariance]
            What method to employ for network matrix estimation.
        input_dir : string
            The directory where the preprocessed derivative needed live
        parcel_name : string
            Shorthand name of the CIFTI label file. 
        parcel_file : string
            The CIFTI label file to use or used to parcellate the brain.
        fishers_r_to_z_transform : [Yes, YES, yes, No, NO, no]
            Should Fisher's r-to-z transformation be applied?
        selected_reg_name : TYPE
            What type of registration do you want to use? Choices are "MSMAll_2_d40_WRN" and "NONE"'
        wavelet : TYPE
            Whether or not to output wavelet matrices 
        combine_resting_scans : TYPE
            If multiple of the same resting state BIDS file type exist should they be combined?

        Raises
        ------
        TypeError
            DESCRIPTION.
        ValueError
            DESCRIPTION.

        Returns
        -------
        layout : TYPE
            DESCRIPTION.
        motion_confounds_filename : TYPE
            DESCRIPTION.
        ICA_outputs : TYPE
            DESCRIPTION.
        graph_theory : TYPE
            DESCRIPTION.
        network_matrix_calculation : TYPE
            DESCRIPTION.
        combine_resting_scans : TYPE
            DESCRIPTION.
        wavelet : TYPE
            DESCRIPTION.

        '''
        # TODO: Need to figure out how to parse kwargs elegantly
        try:
            layout = BIDSLayout(os.path.join(input_dir))
        except: 
            raise TypeError('Input folder: {input_dir} does not look like a BIDS folder or {input_dir} is not a path to a folder. '.format(input_dir=input_dir))
        # check on arguments
        if group == 'participant':
            if not participant_label or not session_label:
                raise ValueError('If "--group participant" is specified, either --participant_label and/or --session_label must also be specified. Exiting.')

        if combine_resting_scans == 'Yes' or combine_resting_scans == 'yes':
            combine_resting_scans = 'YES'
        elif combine_resting_scans == 'No' or combine_resting_scans == 'no':
            combine_resting_scans = 'NO'
        if wavelet == 'Yes' or wavelet == 'yes':
            wavelet = 'YES'
        elif wavelet == 'No' or wavelet == 'no':
            wavelet = 'NO'
        
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
        print('\t-Parcellation file to be used to parcellate outputs: %s' %str(parcel_file))
        print('\t-Short hand parcellation name to be used: %s' %str(parcel_name))
        print('\t-Network matrix metric/s to compute: %s' %(network_matrix_calculation))
        print("\t-Whether or not to compute Fisher's r-to-z transform to network matrices: %s" %(fishers_r_to_z_transform))
        print('\t-Whether or not to perform wavelet entropy on matrices: %s' %(wavelet))
        print('\t-Graph theory metric/s to compute: %s' %(graph_theory))
        print('\t-Input registration file to be used: %s' %str(selected_reg_name))
        print('\t-The preprocessing pipeline that the input comes from: %s' %str(preprocessing_type))
        print('\t-Use ICA outputs: %s' %str(ICA_outputs))
        print('\t-Use mixed effects if multiple of same acquisition: %s' %str(combine_resting_scans))
        print('\n')
        
        return layout,ICA_outputs,graph_theory,network_matrix_calculation,combine_resting_scans,wavelet

