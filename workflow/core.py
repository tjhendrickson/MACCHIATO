#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  7 17:22:55 2020

@author: timothy
"""
from __future__ import print_function
import os
import nibabel
from .connectivity_metrics import NetworkIO, GraphTheoryIO
#from timeseries
import pdb
import argparse
from functools import partial
from subprocess import Popen, PIPE
import subprocess
from mpi4py import MPI
import h5py
import numpy as np
import cifti
from pprint import pprint
from bids.grabbids import BIDSLayout

# specify arguments that MACCHIATO accepts
parser = argparse.ArgumentParser(description='')
parser.add_argument('args', help='The arguments passed from run.py',required=True)
args = parser.parse_args()


class MACCHIATO_setup:
    def __init__(self,group,preprocessing_type,ICA_outputs,
                 graph_theory,network_matrix_calculation,input_dir,parcel_name,
                 parcel_file,fishers_r_to_z_transform,selected_reg_name,
                 timeseries_processing,combine_resting_scans,**kwargs):
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
        selected_reg_name : string
            What type of registration do you want to use? Choices are "MSMAll_2_d40_WRN" and "NONE"'
        timeseries_processing : string
            Whether or not to output wavelet matrices 
        combine_resting_scans : string
            If multiple of the same resting state BIDS file type exist should they be combined?

        Raises
        ------
        TypeError
            DESCRIPTION.
        ValueError
            DESCRIPTION.

        Returns a dictionary with keys and values labeled below
        -------
        layout : bids.grabbids.BIDSLayout object describing BIDS dataset
        ICA_outputs : TYPE
            Parsed ICA_outputs string
        graph_theory : string or list
            Parsed graph theory metrics to calculate downstream
        network_matrix_calculation : string or list
            Parsed network metrics to calculate downstream
        combine_resting_scans : string
            Parsed combine_resting_scan string
        wavelet : string
            Parsed wavelet string
        bolds: list or list of lists
            Path/s to bold timeseries data which will be processed downstream
        
        Example usage
        -------
        >>> parameters = MACCHIATO_setup(group = group,
            preprocessing_type = preprocessing_type,
            ICA_outputs = use_ICA_outputs,
            graph_theory = graph_theory,
            network_matrix_calculation = network_matrix_calculation,
            input_dir = input_dir,
            parcel_name = parcel_name,
            parcel_file = parcel_file,
            fishers_r_to_z_transform = apply_Fishers_r_to_z_transform,
            selected_reg_name = reg_name,
            timeseries_processing = timeseries_processing,
            combine_resting_scans = combine_resting_scans)
        >>> parameters['ICA_outputs']
        >>> 'YES'

        '''
        # place arguments internal to class
        self.group = group
        self.preprocessing_type = preprocessing_type
        self.ICA_outputs = ICA_outputs
        self.graph_theory = graph_theory
        self.network_matrix_calculation = network_matrix_calculation
        self.input_dir = input_dir
        self.parcel_name = parcel_name
        self.parcel_file = parcel_file
        self.fishers_r_to_z_transform = fishers_r_to_z_transform
        self.selected_reg_name = selected_reg_name
        self.msm_all_reg_name="MSMAll_2_d40_WRN"
        self.timeseries_processing = timeseries_processing
        self.combine_resting_scans = combine_resting_scans
        
        for key in ('participant_label', 'session_label'):
                if key in kwargs:
                    setattr(self, key, kwargs[key])
        
        # run workflow logger
        self.workflow_logger()
        
        # create bold lists
        self.create_bold_lists()
        
    def workflow_logger(self):
        '''
        Ensures arguments specified to MACCHIATO are parsed correctly
        Prints specified arguments to standard output (STDOUT)
        '''
        if self.group == 'batch':
            try:
                self.layout = BIDSLayout(os.path.join(self.input_dir))
            except: 
                raise TypeError('Input folder: {input_dir} does not look like a BIDS folder or {input_dir} is not a path to a folder. '.format(input_dir=self.input_dir))
        elif self.group == 'participant':
            assert(self.participant_label or self.session_label), 'If "--group participant" is specified, either --participant_label and/or --session_label must also be specified. Exiting.'
            if self.participant_label and self.session_label:
                try:
                    pdb.set_trace()
                    self.layout = BIDSLayout(os.path.join(self.input_dir,self.participant_label,self.session_label))
                except:
                    raise TypeError('Input folder: {input_dir} does not look like a BIDS folder or {input_dir} is not a path to a folder. '.format(input_dir=os.path.join(self.input_dir,self.participant_label,self.session_label)))
            elif self.participant_label:
                try:
                    self.layout = BIDSLayout(os.path.join(self.input_dir,self.participant_label))
                except:
                    raise TypeError('Input folder: {input_dir} does not look like a BIDS folder or {input_dir} is not a path to a folder. '.format(input_dir=os.path.join(self.input_dir,self.participant_label)))
            # cannot forsee someone just passing session_label, however if that happens this will have to be built out        
            #elif self.session_label
        if self.combine_resting_scans == 'Yes' or self.combine_resting_scans == 'yes':
            self.combine_resting_scans = 'YES'
        elif self.combine_resting_scans == 'No' or self.combine_resting_scans == 'no':
            self.combine_resting_scans = 'NO'
        if not self.timeseries_processing == 'NONE' or self.timeseries_processing == 'yes':
            self.timeseries_processing = self.timeseries_processing[0]
        
         # use ICA outputs
        if self.ICA_outputs == 'yes' or self.ICA_outputs == 'Yes':
            self.ICA_outputs = 'YES'
        else:
            self.ICA_outputs = 'NO'
        # if all graph theory metrics are requested, transform arg from string to list
        if not self.graph_theory == 'NONE' :
            if self.graph_theory == 'All':
                self.graph_theory = ['clustering_coefficient','local_efficiency','strength','node_betweenness_centrality', 'edge_betweenness_centrality', 'eigenvector_centrality']
            else:
                self.graph_theory = self.graph_theory[0]
                
        if self.network_matrix_calculation == 'All':
            self.network_matrix_calculation = ['correlation','partial_correlation','dynamic_time_warping',
                                               'tangent','covariance', 'precision',
                                               'sparse_inverse_precision','sparse_inverse_covariance']
        print("Running MACCHIATO ")
        print('\t-Level that MACCHIATO will run on: %s' %(str(self.group)))
        if self.group == 'participant':
            if self.participant_label:
                print('\t-Participant ID: %s' %str(self.participant_label))
            if self.session_label:
                print('\t-Session ID: %s' %str(self.session_label))
        print('\t-Parcellation file to be used to parcellate outputs: %s' %str(self.parcel_file))
        print('\t-Short hand parcellation name to be used: %s' %str(self.parcel_name))
        print('\t-Network matrix metric/s to compute: %s' %(self.network_matrix_calculation))
        print("\t-Whether or not to compute Fisher's r-to-z transform to network matrices: %s" %(self.fishers_r_to_z_transform))
        print('\t-Timeseries processing method to apply: %s' %(self.timeseries_processing))
        print('\t-Graph theory metric/s to compute: %s' %(self.graph_theory))
        print('\t-Input registration file to be used: %s' %str(self.selected_reg_name))
        print('\t-The preprocessing pipeline that the input comes from: %s' %str(self.preprocessing_type))
        print('\t-Use ICA outputs: %s' %str(self.ICA_outputs))
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
        
        if self.combine_resting_scans == 'NO':
            if self.preprocessing_type == 'HCP':
                # use ICA outputs
                if self.ICA_outputs == 'YES':
                    self.ICA_string="_FIXclean"
                    if self.selected_reg_name == self.msm_all_reg_name:
                        if self.group == 'batch':
                            self.bolds = [f.filename for f in self.layout.get(type='clean',extensions="dtseries.nii",task='rest') if self.msm_all_reg_name+'_hp2000_clean' in f.filename]
                        elif self.group == 'participant':
                            pass 
                    else:
                        if self.group == 'batch':
                            self.bolds = [f.filename for f in self.layout.get(type='clean',extensions="dtseries.nii", task='rest') if '_hp2000_clean' and not self.msm_all_reg_name in f.filename]
                # do not use ICA outputs
                else:
                    self.ICA_string=""
                    if self.selected_reg_name == self.msm_all_reg_name:
                        if self.group == 'batch':
                            self.bolds = [f.filename for f in self.layout.get(extensions="dtseries.nii", task='rest') if self.msm_all_reg_name + '_hp2000' in f.filename and not 'clean' in f.filename]
                    else:
                        if self.group == 'batch':
                            self.bolds = [f.filename for f in self.layout.get(extensions="dtseries.nii", task='rest') if '_hp2000' in f.filename and not 'clean' and not self.msm_all_reg_name in f.filename]
            elif self.preprocessing_type == 'fmriprep':
                #use ICA outputs
                if self.ICA_outputs == 'YES':
                    self.ICA_string="_AROMAclean"
                    if self.group == 'batch':
                        self.bolds = [f.filename for f in self.layout.get(type='bold',task='rest') if 'smoothAROMAnonaggr' in f.filename]
                # do not use ICA outputs
                else:
                    self.ICA_string=""
                    if self.group == 'batch':
                        self.bolds = [f.filename for f in self.layout.get(type='bold',task='rest') if 'preproc' in f.filename]
                self.bolds_ref = [f.filename for f in self.layout.get(type='boldref',task='rest')]
        else:
            self.combined_bolds_list = []
            if self.layout.get_sessions() > 0:
                for scanning_session in self.layout.get_sessions():
                    # retreive subject id that is associated with session id and parse data with subject and session id
                    for subject in self.layout.get_subjects(session=scanning_session):
                        if self.preprocessing_type == 'HCP':
                            # use ICA outputs
                            if self.ICA_outputs == 'YES':
                                self.ICA_string="_FIXclean"
                                if self.selected_reg_name == self.msm_all_reg_name:
                                    self.bolds = [f.filename for f in self.layout.get(type='clean',extensions="dtseries.nii",task='rest',subject=subject,session=scanning_session) if self.msm_all_reg_name+'_hp2000_clean' in f.filename]
                                else:
                                    self.bolds = [f.filename for f in self.layout.get(type='clean',extensions="dtseries.nii", task='rest',subject=subject,session=scanning_session) if '_hp2000_clean' and not self.msm_all_reg_name in f.filename]
                            # do not use ICA outputs
                            else:
                                self.ICA_string=""
                                if self.selected_reg_name == self.msm_all_reg_name:
                                    self.bolds = [f.filename for f in self.layout.get(extensions="dtseries.nii", task='rest',subject=subject,session=scanning_session) if self.msm_all_reg_name + '_hp2000' in f.filename and not 'clean' in f.filename]
                                else:
                                    self.bolds = [f.filename for f in self.layout.get(extensions="dtseries.nii", task='rest',subject=subject,session=scanning_session) if '_hp2000' in f.filename and not 'clean' and not self.msm_all_reg_name in f.filename]
                        elif self.preprocessing_type == 'fmriprep':
                            #use ICA outputs
                            if self.ICA_outputs == 'YES':
                                self.ICA_string="_AROMAclean"
                                self.bolds = [f.filename for f in self.layout.get(type='bold',task='rest',subject=subject,session=scanning_session) if 'smoothAROMAnonaggr' in f.filename]
                            # do not use ICA outputs
                            else:
                                self.ICA_string=""
                                self.bolds = [f.filename for f in self.layout.get(type='bold',task='rest') if 'preproc' in f.filename]
                            self.bolds_ref = [f.filename for f in self.layout.get(type='boldref',task='rest')]
                        if len(self.bolds) == 2:
                            self.combined_bolds_list.append(self.bolds)
            else:
                for scanning_session in self.layout.get_subjects():
                    if self.preprocessing_type == 'HCP':
                        # use ICA outputs
                        if self.ICA_outputs == 'YES':
                            self.ICA_string="_FIXclean"
                            if self.selected_reg_name == self.msm_all_reg_name:
                                self.bolds = [f.filename for f in self.layout.get(type='clean',extensions="dtseries.nii",task='rest',subject=scanning_session) if self.msm_all_reg_name+'_hp2000_clean' in f.filename]
                            else:
                                self.bolds = [f.filename for f in self.layout.get(type='clean',extensions="dtseries.nii", task='rest',subject=scanning_session) if '_hp2000_clean' and not self.msm_all_reg_name in f.filename]
                        # do not use ICA outputs
                        else:
                            self.ICA_string=""
                            if self.selected_reg_name == self.msm_all_reg_name:
                                self.bolds = [f.filename for f in self.layout.get(extensions="dtseries.nii", task='rest',subject=scanning_session) if self.msm_all_reg_name + '_hp2000' in f.filename and not 'clean' in f.filename]
                            else:
                                self.bolds = [f.filename for f in self.layout.get(extensions="dtseries.nii", task='rest',subject=scanning_session) if '_hp2000' in f.filename and not 'clean' and not self.msm_all_reg_name in f.filename]
                    elif self.preprocessing_type == 'fmriprep':
                        #use ICA outputs
                        if self.ICA_outputs == 'YES':
                            self.ICA_string="_AROMAclean"
                            self.bolds = [f.filename for f in self.layout.get(type='bold',task='rest',subject=scanning_session) if 'smoothAROMAnonaggr' in f.filename]
                        # do not use ICA outputs
                        else:
                            self.ICA_string=""
                            self.bolds = [f.filename for f in self.layout.get(type='bold',task='rest') if 'preproc' in f.filename]
                        self.bolds_ref = [f.filename for f in self.layout.get(type='boldref',task='rest')]
                    if len(self.bolds) == 2:
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
            read_parcel_file = cifti.read(args.parcellation_file)
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
        t_start = MPI.Wtime() 
        
        with h5py.File(os.path.join(args.output_dir,'data.hdf5'),'w',driver='mpio',comm=comm) as hdf5:
            network_matrix_dset = hdf5.create_dataset(name=setupParams['network_matrix_calculation'], shape=(image_count,height,width), dtype='f')
            network_matrix_dset.attrs['parcel_file'] = args.parcellation_file
            network_matrix_dset.attrs['parcel_name'] = args.parcellation_name    
            for i in range(rank, image_count, comm.size):
                bold = setupParams['bolds'][i]
                metric_data = execute_MACCHIATO_instances(
                    preprocessing_type=args.preprocessing_type,
                    parcel_file=args.parcellation_file,
                    parcel_name=args.parcellation_name,
                    selected_reg_name=args.reg_name,
                    ICA_outputs=setupParams['ICA_outputs'],
                      fishers_r_to_z_transform=args.apply_Fishers_r_to_z_transform,
                    combine_resting_scans=setupParams['combine_resting_scans'],
                    timeseries_processing=setupParams['timeseries_processing'],
                    network_metric=setupParams['network_matrix_calculation'],
                    output_dir=args.output_dir, 
                    graph_theory_metric=setupParams['graph_theory'],
                    bold=bold)
                network_matrix_dset[i,:,:] = metric_data
            """
            execute_MACCHIATO_instances(preprocessing_type=args.preprocessing_type,
                    parcel_file=args.parcellation_file,
                    parcel_name=args.parcellation_name,
                    selected_reg_name=args.reg_name,
                    ICA_outputs=setupParams['ICA_outputs'],
                    fishers_r_to_z_transform=args.apply_Fishers_r_to_z_transform,
                    combine_resting_scans=setupParams['combine_resting_scans'],
                    wavelet=setupParams['wavelet'],
                    network_metric=setupParams['network_matrix_calculation'],
                    output_dir=args.output_dir, 
                    graph_theory_metric=setupParams['graph_theory'],
                    bold=setupParams['bolds'][0])
            """

#handles all logic of inputted data to run MACCHIATO        
def execute_MACCHIATO_instances(bold,preprocessing_type,parcel_file,parcel_name,
                  selected_reg_name,ICA_outputs,
                  combine_resting_scans,timeseries_processing,network_metric,graph_theory_metric,
                  output_dir,fishers_r_to_z_transform):
       
    if timeseries_processing == 'YES':
        pass #TODO will be passed to timeseries_metrics.py
    else:
        network_metric_init = NetworkIO(output_dir=output_dir, 
                                        cifti_file=bold, 
                                        parcel_file=parcel_file, 
                                        parcel_name=parcel_name, 
                                        network_metric=network_metric,
                                        fishers_r_to_z_transform=fishers_r_to_z_transform)
        if graph_theory_metric == 'NONE':
            if type(network_metric) == str:
                metric_data = network_metric_init.create_network_matrix()
                return metric_data # TODO this logic will have to be much stronger with full implementation
            else:
                pass
        else:
            GraphTheoryIO(output_dir=output_dir, 
                          cifti_file=bold, 
                          parcel_file=parcel_file, 
                          parcel_name=parcel_name,
                          network_metric=network_metric,
                          fishers_r_to_z_transform=fishers_r_to_z_transform,
                          graph_theory_metric=graph_theory_metric)
        
