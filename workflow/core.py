#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  7 17:22:55 2020

@author: timothy
"""
from __future__ import print_function
import argparse
from bids.grabbids import BIDSLayout
import cifti
from datetime import date
from functools import partial
import h5py
from mpi4py import MPI
import nibabel
import numpy as np
import os
from pprint import pprint
import pandas as pd
from subprocess import Popen, PIPE
import sys

sys.path.append('../')
from tools.connectivity_metrics import NetworkIO, GraphTheoryIO

# specify arguments that MACCHIATO accepts

class MACCHIATO_setup:
    def __init__(self,args_dict):

        '''
        Parameters
        ----------
        input_dir : string
            The directory where the preprocessed derivative needed live
        output_dir : string
            The directory where the outputs from MACCHIATO will be placed
        analysis_level : [participant]
            Run participant by participant (part of BIDS specification) 
        preprocessing_type : [HCP, fmriprep]
            BIDS-apps preprocessing pipeline run on data. Choices include "HCP" and "fmriprep". 
        denoised_outputs : [Yes, yes, No, no]
            Use denoised (whether FIX, AROMA, or scrubbed depending on pipeline) outputs for network matrix estimation.
        combine_resting_scans : string
            If multiple of the same resting state BIDS file type exist should they be combined?
        parcellation_file : string
            The CIFTI label file to use or used to parcellate the brain.
        parcellation_name : string
            Shorthand name of the CIFTI label file. 
        reg_name : string
            What type of registration do you want to use? Choices are "MSMAll_2_d40_WRN" and "NONE"'
        fishers_r_to_z_transform : [Yes, YES, yes, No, NO, no]
            Should Fisher's r-to-z transformation be applied?
        network_matrix_calculation : [All, all, correlation, partial_correlation, dynamic_time_warping, covariance, precision, sparse_inverse_precision, sparse_inverse_covariance]
            What method to employ for network matrix estimation.


        Raises
        ------
        TypeError
            DESCRIPTION.
        ValueError
            DESCRIPTION.

        Returns a dictionary with keys and values labeled below
        -------
        layout : bids.grabbids.BIDSLayout object describing BIDS dataset
        denoised_outputs : TYPE
            Parsed denoised_outputs string
        network_matrix_calculation : string or list
            Parsed network metrics to calculate downstream
        combine_resting_scans : string
            Parsed combine_resting_scan string
        bolds: list or list of lists
            Path/s to bold timeseries data which will be processed downstream
        
        Example usage
        -------
        >>> parameters = MACCHIATO_setup(group = group,
            preprocessing_type = preprocessing_type,
            denoised_outputs = use_denoised_outputs,
            network_matrix_calculation = network_matrix_calculation,
            input_dir = input_dir,
            parcel_name = parcel_name,
            parcel_file = parcel_file,
            fishers_r_to_z_transform = apply_Fishers_r_to_z_transform,
            selected_reg_name = reg_name,
            combine_resting_scans = combine_resting_scans)
        >>> parameters['denoised_outputs']
        >>> 'YES'

        '''
        # place arguments internal to class
        self.input_dir = args_dict.get('--input_dir')
        self.output_dir = args_dict.get('--output_dir')
        self.analysis_level = args_dict.get('--analysis_level')
        self.preprocessing_type = args_dict.get('--preprocessing_type')
        self.denoised_outputs = args_dict.get('--use_denoised_outputs')
        self.combine_resting_scans = args_dict.get('--combine_resting_scan')
        self.parcellation_file = args_dict.get('--parcellation_file')
        self.parcellation_name = args_dict.get('--parcellation_name')
        self.selected_reg_name = args_dict.get('--reg_name')
        self.apply_Fishers_r_to_z_transform = args_dict.get('--apply_Fishers_r_to_z_transform')
        self.network_matrix_calculation = args_dict.get('--network_matrix_calculation')
        self.msm_all_reg_name="MSMAll_2_d40_WRN"
        if '--participant_label' in args_dict:
            self.participant_label = args_dict.get('--participant_label')
        else:
            self.participant_label = None    
        if '--session_label' in args_dict:
            self.session_label = args_dict.get('--session_label')
        else:
            self.session_label = None

        # run workflow logger
        self.workflow_logger()
        
        # create bold lists
        self.create_bold_lists()

        # generate HDF5 outputs
        self.execute_MACCHIATO_instances()
        
    def workflow_logger(self):
        '''
        Ensures arguments specified to MACCHIATO are parsed correctly
        Prints specified arguments to standard output (STDOUT)
        '''
        try:
            self.layout = BIDSLayout(os.path.join(self.input_dir))
        except: 
            raise TypeError('Input folder: {input_dir} does not look like a BIDS folder or {input_dir} is not a path to a folder. '.format(input_dir=self.input_dir))
        if self.combine_resting_scans == 'Yes' or self.combine_resting_scans == 'yes':
            self.combine_resting_scans = 'YES'
        elif self.combine_resting_scans == 'No' or self.combine_resting_scans == 'no':
            self.combine_resting_scans = 'NO'
        
         # use ICA outputs
        if self.denoised_outputs == 'yes' or self.denoised_outputs == 'Yes':
            self.denoised_outputs = 'YES'
        else:
            self.denoised_outputs = 'NO'
                
        if self.network_matrix_calculation == 'All':
            self.network_matrix_calculation = ['correlation','partial_correlation','dynamic_time_warping',
                                               'tangent','covariance', 'precision',
                                               'sparse_inverse_precision','sparse_inverse_covariance']
        print("Running MACCHIATO ")
        print('\t-Level that MACCHIATO will run on: %s' %(str(self.analysis_level)))
        if self.participant_label:
            print('\t-Participant ID: %s' %str(self.participant_label))
        if self.session_label:
            print('\t-Session ID: %s' %str(self.session_label))
        print('\t-Parcellation file to be used to parcellate outputs: %s' %str(self.parcellation_file))
        print('\t-Short hand parcellation name to be used: %s' %str(self.parcellation_name))
        print('\t-Network matrix metric/s to compute: %s' %(self.network_matrix_calculation))
        print("\t-Whether or not to compute Fisher's r-to-z transform to network matrices: %s" %(self.apply_Fishers_r_to_z_transform))

        print('\t-Input registration file to be used: %s' %str(self.selected_reg_name))
        print('\t-The preprocessing pipeline that the input comes from: %s' %str(self.preprocessing_type))
        print('\t-Use ICA outputs: %s' %str(self.denoised_outputs))
        print('\t-Combine matrices/timeseries from resting state pairs within the same session: %s' %str(self.combine_resting_scans))
        print('\n')
        
    def create_bold_lists(self):
        
        '''
        Generates list or list of lists of bold timeseries paths for downstream processing

        Returns
        -------
        TYPE list
            Paths to bold timeseries
        '''
        # TODO need to incorporate the instances in which participant label and session label BOTH exist
        if self.combine_resting_scans == 'NO':
            if self.preprocessing_type == 'HCP':
                # use ICA outputs
                if self.denoised_outputs == 'YES':
                    self.ICA_string="_FIXclean"
                    if self.selected_reg_name == self.msm_all_reg_name:
                        if not self.participant_label and not self.session_label:
                            self.bolds = [f.filename for f in self.layout.get(type='clean',extensions="dtseries.nii") if self.msm_all_reg_name+'_hp2000_clean' in f.filename and 'task-rest' in f.filename]
                        elif self.participant_label and not self.session_label:
                            self.bolds = [f.filename for f in self.layout.get(subject=self.participant_label,type='clean',extensions="dtseries.nii") if self.msm_all_reg_name+'_hp2000_clean' in f.filename and 'task-rest' in f.filename] 
                        elif not self.participant_label and self.session_label:
                            self.bolds = [f.filename for f in self.layout.get(session=self.session_label,type='clean',extensions="dtseries.nii") if self.msm_all_reg_name+'_hp2000_clean' in f.filename and 'task-rest' in f.filename]
                        else:
                            self.bolds = [f.filename for f in self.layout.get(subject=self.participant_label,session=self.session_label,type='clean',extensions="dtseries.nii") if self.msm_all_reg_name+'_hp2000_clean' in f.filename and 'task-rest' in f.filename]
                    else:
                        if not self.participant_label and not self.session_label:
                            self.bolds = [f.filename for f in self.layout.get(type='clean',extensions="dtseries.nii") if '_hp2000_clean' in f.filename and 'task-rest' in f.filename and not self.msm_all_reg_name in f.filename]
                        elif self.participant_label and not self.session_label:
                            self.bolds = [f.filename for f in self.layout.get(subject=self.participant_label,type='clean',extensions="dtseries.nii") if '_hp2000_clean' in f.filename and 'task-rest' in f.filename and not self.msm_all_reg_name in f.filename]
                        elif not self.participant_label and self.session_label:
                            self.bolds = [f.filename for f in self.layout.get(session=self.session_label,type='clean',extensions="dtseries.nii",) if self.msm_all_reg_name+'_hp2000_clean' in f.filename and 'task-rest' in f.filename]
                        else:
                            self.bolds = [f.filename for f in self.layout.get(subject=self.participant_label,session=self.session_label,type='clean',extensions="dtseries.nii") if self.msm_all_reg_name+'_hp2000_clean' in f.filename and 'task-rest' in f.filename]

                # do not use ICA outputs
                else:
                    self.ICA_string=""
                    if self.selected_reg_name == self.msm_all_reg_name:
                        if not self.participant_label and not self.session_label:
                            self.bolds = [f.filename for f in self.layout.get(extensions="dtseries.nii") if self.msm_all_reg_name + '_hp2000' in f.filename and 'task-rest' in f.filename and not 'clean' in f.filename]
                        elif self.participant_label and not self.session_label:
                            self.bolds = [f.filename for f in self.layout.get(subject=self.participant_label,extensions="dtseries.nii") if self.msm_all_reg_name + '_hp2000' in f.filename and 'task-rest' in f.filename and not 'clean' in f.filename]
                        elif not self.participant_label and self.session_label:
                            self.bolds = [f.filename for f in self.layout.get(session=self.session_label,extensions="dtseries.nii") if self.msm_all_reg_name + '_hp2000' in f.filename and 'task-rest' in f.filename and not 'clean' in f.filename]
                        else:
                            self.bolds = [f.filename for f in self.layout.get(subject=self.participant_label,session=self.session_label,extensions="dtseries.nii") if self.msm_all_reg_name + '_hp2000' in f.filename and 'task-rest' in f.filename and not 'clean' in f.filename]
                    else:
                        if not self.participant_label and not self.session_label:
                            self.bolds = [f.filename for f in self.layout.get(extensions="dtseries.nii") if '_hp2000' in f.filename and not 'clean' in f.filename and 'task-rest' in f.filename and not self.msm_all_reg_name in f.filename]
                        elif self.participant_label and not self.session_label:
                            self.bolds = [f.filename for f in self.layout.get(subject=self.participant_label,extensions="dtseries.nii") if '_hp2000' in f.filename and 'task-rest' in f.filename and not 'clean' in f.filename and not self.msm_all_reg_name in f.filename]
                        elif not self.participant_label and self.session_label:
                            self.bolds = [f.filename for f in self.layout.get(session=self.session_label,extensions="dtseries.nii") if '_hp2000' in f.filename and 'task-rest' in f.filename and not 'clean' in f.filename and not self.msm_all_reg_name in f.filename]
                        else:
                            self.bolds = [f.filename for f in self.layout.get(subject=self.participant_label,session=self.session_label,extensions="dtseries.nii") if '_hp2000' in f.filename and 'task-rest' in f.filename and not 'clean' in f.filename and not self.msm_all_reg_name in f.filename]
            elif self.preprocessing_type == 'fmriprep':
                #use ICA outputs
                if self.denoised_outputs == 'YES':
                    self.ICA_string="_AROMAclean"
                    if not self.participant_label and not self.session_label:
                        self.bolds = [f.filename for f in self.layout.get(type='bold') if 'smoothAROMAnonaggr' in f.filename]
                    elif self.participant_label and not self.session_label:
                        self.bolds = [f.filename for f in self.layout.get(subject=self.participant_label,type='bold') if 'smoothAROMAnonaggr' in f.filename and 'task-rest' in f.filename]
                    elif not self.participant_label and self.session_label:
                        self.bolds = [f.filename for f in self.layout.get(session=self.session_label,type='bold') if 'smoothAROMAnonaggr' in f.filename and 'task-rest' in f.filename]
                    else:
                        self.bolds = [f.filename for f in self.layout.get(subject=self.participant_label,session=self.session_label,type='bold') if 'smoothAROMAnonaggr' in f.filename and 'task-rest' in f.filename]
                # do not use ICA outputs
                else:
                    self.ICA_string=""
                    if not self.participant_label and not self.session_label:
                        self.bolds = [f.filename for f in self.layout.get(type='bold') if 'preproc' in f.filename]
                    elif self.participant_label and not self.session_label:
                        self.bolds = [f.filename for f in self.layout.get(subject=self.participant_label,type='bold') if 'preproc' in f.filename and 'task-rest' in f.filename]
                    elif not self.participant_label and self.session_label:
                        self.bolds = [f.filename for f in self.layout.get(session =self.session_label,type='bold') if 'preproc' in f.filename and 'task-rest' in f.filename]
                    else:
                        self.bolds = [f.filename for f in self.layout.get(subject=self.participant_label,session =self.session_label,type='bold') if 'preproc' in f.filename and 'task-rest' in f.filename]
                self.bolds_ref = [f.filename for f in self.layout.get(type='boldref')]
            elif self.preprocessing_type == 'ABCD':
                 # use filered outputs
                if self.denoised_outputs == 'YES':
                    self.ICA_string="_filtered"
                    if not self.participant_label and not self.session_label:
                        self.bolds = [f.filename for f in self.layout.get(extensions="dtseries.nii",type='Atlas') if 'task-rest' in f.filename and 'DCANBOLDProc' in f.filename if 'DCANBOLDProc' in os.path.basename(os.path.dirname(f.filename))] # traditional outputs
                        #if not self.bolds > 0:
                        #    self.bolds = [f.filename for f in self.layout.get(extensions="dtseries.nii",type='timeseries') if 'task-rest' and 'run' in f.filename] # truncated outputs generated for NDA sumbission
                        #    filtered = 'NO' 
                        #else:
                        #filtered = 'YES'
                    elif self.participant_label and not self.session_label:
                        self.bolds = [f.filename for f in self.layout.get(subject=self.participant_label,extensions="dtseries.nii",type='Atlas') if 'task-rest' in f.filename and 'DCANBOLDProc' in f.filename if 'DCANBOLDProc' in os.path.basename(os.path.dirname(f.filename))] # traditional outputs
                        #if not self.bolds > 0:
                        #    self.bolds = [f.filename for f in self.layout.get(subject=self.participant_label,extensions="dtseries.nii",type='timeseries') if 'task-rest' and 'run' in f.filename] # truncated outputs generated for NDA sumbission
                        #    filtered = 'NO' 
                        #else:
                        #filtered = 'YES'
                    elif not self.participant_label and self.session_label:
                        self.bolds = [f.filename for f in self.layout.get(session=self.session_label,extensions="dtseries.nii",type='Atlas') if 'task-rest' in f.filename and 'DCANBOLDProc' in f.filename if 'DCANBOLDProc' in os.path.basename(os.path.dirname(f.filename))] # traditional outputs
                        #if not self.bolds > 0:
                        #    self.bolds = [f.filename for f in self.layout.get(session=self.session_label,extensions="dtseries.nii",type='timeseries') if 'task-rest' and 'run' in f.filename] # truncated outputs generated for NDA sumbission
                        #    filtered = 'NO' 
                        #else:
                        #filtered = 'YES'
                    else:
                        self.bolds = [f.filename for f in self.layout.get(subject=self.participant_label,session=self.session_label,extensions="dtseries.nii",type='Atlas') if 'task-rest' in f.filename and 'DCANBOLDProc' in f.filename if 'DCANBOLDProc' in os.path.basename(os.path.dirname(f.filename))] # traditional outputs
                        #if not self.bolds > 0:
                        #    self.bolds = [f.filename for f in self.layout.get(subject=self.participant_label,session=self.session_label,extensions="dtseries.nii",type='timeseries') if 'task-rest' and 'run' in f.filename] # truncated outputs generated for NDA sumbission
                        #    filtered = 'NO'
                        #else:
                        #filtered = 'YES'
                    #if filtered == 'NO':
                    #    raise NotImplementedError # TODO: will need to add this, however, the use cases in which they may happen are limited
                        
                else: 
                    if not self.participant_label and not self.session_label:
                        self.bolds = [f.filename for f in self.layout.get(extensions="dtseries.nii",type='timeseries') if 'task-rest' in f.filename and 'run' in f.filename] # truncated outputs generated for NDA sumbission
                        #filtered = 'NO'
                    elif self.participant_label and not self.session_label:
                        self.bolds = [f.filename for f in self.layout.get(subject=self.participant_label,extensions="dtseries.nii",type='timeseries') if 'task-rest' in f.filename and 'run' in f.filename] # truncated outputs generated for NDA sumbission
                        #filtered = 'NO' 
                    elif not self.participant_label and self.session_label:
                        self.bolds = [f.filename for f in self.layout.get(session=self.session_label,extensions="dtseries.nii",type='timeseries') if 'task-rest' in f.filename and 'run' in f.filename] # truncated outputs generated for NDA sumbission
                        #filtered = 'NO' 
                    else:
                        self.bolds = [f.filename for f in self.layout.get(subject=self.participant_label,session=self.session_label,extensions="dtseries.nii",type='timeseries') if 'task-rest' in f.filename and 'run' in f.filename] # truncated outputs generated for NDA sumbission
                        #filtered = 'NO'
        else:
            self.combined_bolds_list = [] # will combine within one scanning session no matter how many resting state scans collected
            if len(self.layout.get_sessions()) > 0:
                # retreive sessions data based on inputted (or not) participant and session labels
                if not self.participant_label and not self.session_label:
                    scanning_sessions = self.layout.get_sessions()
                elif self.participant_label and not self.session_label:
                    scanning_sessions = self.layout.get_sessions(subject=self.participant_label)
                elif not self.participant_label and self.session_label:
                    scanning_sessions = self.layout.get_sessions(session=self.session_label)
                else:
                    scanning_sessions = self.layout.get_sessions(session=self.session_label, subject=self.participant_label)
                for subject in self.layout.get_subjects(session=scanning_sessions):
                    for scanning_session in self.layout.get_sessions(subject=subject):
                        if self.preprocessing_type == 'HCP':
                            # use ICA outputs
                            if self.denoised_outputs == 'YES':
                                self.ICA_string="_FIXclean"
                                if self.selected_reg_name == self.msm_all_reg_name:
                                    self.bolds = [f.filename for f in self.layout.get(type='clean',extensions="dtseries.nii",subject=subject,session=scanning_session) if self.msm_all_reg_name+'_hp2000_clean' in f.filename]
                                else:
                                    self.bolds = [f.filename for f in self.layout.get(type='clean',extensions="dtseries.nii",subject=subject,session=scanning_session) if '_hp2000_clean' in f.filename and not self.msm_all_reg_name in f.filename]
                            # do not use ICA outputs
                            else:
                                self.ICA_string=""
                                if self.selected_reg_name == self.msm_all_reg_name:
                                    self.bolds = [f.filename for f in self.layout.get(extensions="dtseries.nii",subject=subject,session=scanning_session) if self.msm_all_reg_name + '_hp2000' in f.filename and not 'clean' in f.filename]
                                else:
                                    self.bolds = [f.filename for f in self.layout.get(extensions="dtseries.nii",subject=subject,session=scanning_session) if '_hp2000' in f.filename and not 'clean' in f.filename and not self.msm_all_reg_name in f.filename]
                        elif self.preprocessing_type == 'fmriprep':
                            #use ICA outputs
                            if self.denoised_outputs == 'YES':
                                self.ICA_string="_AROMAclean"
                                self.bolds = [f.filename for f in self.layout.get(type='bold',subject=subject,session=scanning_session) if 'smoothAROMAnonaggr' in f.filename]
                            # do not use ICA outputs
                            else:
                                self.ICA_string=""
                                self.bolds = [f.filename for f in self.layout.get(type='bold') if 'preproc' in f.filename]
                            self.bolds_ref = [f.filename for f in self.layout.get(type='boldref')]
                        elif self.preprocessing_type == 'ABCD':    
                            if self.denoised_outputs == 'YES':
                                self.ICA_string='_filtered'
                                self.bolds = [f.filename for f in self.layout.get(subject=subject,session=scanning_session,type='timeseries',extensions="dtseries.nii") if 'desc-filtered' in f.filename and 'task-rest' in f.filename]
                                if not len(self.bolds) > 0:
                                    self.bolds = [f.filename for f in self.layout.get(subject=subject,session=scanning_session,type='Atlas',extensions="dtseries.nii") if 'DCANBOLDProc' in f.filename and 'task-rest' in f.filename if not 'DCANBOLDProc' in os.path.basename(os.path.dirname(f.filename))]
                            else:
                                self.ICA_string=''
                                self.bolds = [f.filename for f in self.layout.get(subject=subject,session=scanning_session,extensions="dtseries.nii",type='timeseries') if 'task-rest' in f.filename and 'run' in f.filename] # truncated outputs generated for NDA sumbission    
                        self.combined_bolds_list.append(self.bolds)
            else:
                # retreive subject data based on inputted (or not) participant label
                if not self.participant_label:
                    scanning_sessions = self.layout.get_subjects()
                elif self.participant_label:
                    scanning_sessions = self.layout.get_subjects(subject=self.participant_label)
                for scanning_session in scanning_sessions:
                    if self.preprocessing_type == 'HCP':
                        # use ICA outputs
                        if self.denoised_outputs == 'YES':
                            self.ICA_string="_FIXclean"
                            if self.selected_reg_name == self.msm_all_reg_name:
                                self.bolds = [f.filename for f in self.layout.get(type='clean',extensions="dtseries.nii",subject=scanning_session) if self.msm_all_reg_name+'_hp2000_clean' in f.filename]
                            else:
                                self.bolds = [f.filename for f in self.layout.get(type='clean',extensions="dtseries.nii",subject=scanning_session) if '_hp2000_clean' in f.filename and not self.msm_all_reg_name in f.filename]
                        # do not use ICA outputs
                        else:
                            self.ICA_string=""
                            if self.selected_reg_name == self.msm_all_reg_name:
                                self.bolds = [f.filename for f in self.layout.get(extensions="dtseries.nii",subject=scanning_session) if self.msm_all_reg_name + '_hp2000' in f.filename and not 'clean' in f.filename]
                            else:
                                self.bolds = [f.filename for f in self.layout.get(extensions="dtseries.nii",subject=scanning_session) if '_hp2000' in f.filename and not 'clean' in f.filename and not self.msm_all_reg_name in f.filename]
                    elif self.preprocessing_type == 'fmriprep':
                        #use ICA outputs
                        if self.denoised_outputs == 'YES':
                            self.ICA_string="_AROMAclean"
                            self.bolds = [f.filename for f in self.layout.get(type='bold',subject=scanning_session) if 'smoothAROMAnonaggr' in f.filename]
                        # do not use ICA outputs
                        else:
                            self.ICA_string=""
                            self.bolds = [f.filename for f in self.layout.get(type='bold') if 'preproc' in f.filename]
                        self.bolds_ref = [f.filename for f in self.layout.get(type='boldref')]
                    self.combined_bolds_list.append(self.bolds)
        self.bolds = self.combined_bolds_list
    
    # TODO incorporate variables into class
    def execute_MACCHIATO_instances(self):
        # now parse arguments and print to standard output (STDOUT), slight difference in call if done as batch or participant
        
        # set up multiprocessing/parallelization allocation
        comm = MPI.COMM_WORLD
        rank = comm.rank
        
        # retreive number of, height and width of matrices 
        image_count = len(self.bolds)
        try:
            read_parcel_file = cifti.read(self.parcellation_file)
        except TypeError:
            print('This does not look like a CIFTI parcel file. Must exit')
        parcel_file_label_tuple = read_parcel_file[1][0][0][1]
        parcel_labels = []
        for value in parcel_file_label_tuple:
                if not '???' in parcel_file_label_tuple[value][0]:
                    parcel_labels.append(parcel_file_label_tuple[value][0])
        height = len(parcel_labels) # height of functional connectome
        width = len(parcel_labels) # width of functional connectome
        pprint("============================================================================")
        pprint(" Running %d parallel MPI processes" % comm.size)
        pprint(" Processing %d images of size %d x %d" % (image_count, width, height))
        
        comm.Barrier()    ### Start Stopwatch ###
        today=date.today()
        pretty_today=today.strftime("%b-%d-%Y")
        with h5py.File(os.path.join(self.output_dir,'MACCHIATO_output_'+pretty_today+'.hdf5'),'w',driver='mpio',comm=comm) as hdf5:
            if type(self.network_matrix_calculation) == str:
                network_matrix_dset = hdf5.create_dataset(name=self.network_matrix_calculation, shape=(image_count,height,width), dtype='f')
                network_matrix_dset.attrs['parcel_file'] = self.parcellation_file
                network_matrix_dset.attrs['parcel_name'] = self.parcellation_name    
                for i in range(rank, image_count, comm.size):                  
                    bold = self.bolds[i]
                    network_metric_init = NetworkIO(output_dir=self.output_dir, 
                                        cifti_data=bold, 
                                        parcel_file=self.parcellation_file, 
                                        parcel_name=self.parcellation_name, 
                                        network_metric=self.network_matrix_calculation,
                                        fishers_r_to_z_transform=self.apply_Fishers_r_to_z_transform)
                    
                    metric_data = network_metric_init.create_network_matrix()
                    network_matrix_dset[i,:,:]=metric_data
            else:
                for j in self.network_matrix_calculation:
                    network_matrix_dset = hdf5.create_dataset(name=j, shape=(image_count,height,width), dtype='f')
                    network_matrix_dset.attrs['parcel_file'] = self.parcellation_file
                    network_matrix_dset.attrs['parcel_name'] = self.parcellation_name    
                    for i in range(rank, image_count, comm.size):                  
                        bold = self.bolds[i]
                        network_metric_init = NetworkIO(output_dir=self.output_dir, 
                                            cifti_data=bold, 
                                            parcel_file=self.parcellation_file, 
                                            parcel_name=self.parcellation_name, 
                                            network_metric=self.network_matrix_calculation,
                                            fishers_r_to_z_transform=self.apply_Fishers_r_to_z_transform)
                        
                        metric_data = network_metric_init.create_network_matrix()
                        network_matrix_dset[i,:,:]=metric_data
                
if __name__ == '__main__':
    arg_string = sys.argv[1:]
    args_dict = {}
    for arg in arg_string:
        key=arg.split('=')[0]
        value=arg.split('=')[1]
        args_dict[key] = value
    MACCHIATO_setup(args_dict)
