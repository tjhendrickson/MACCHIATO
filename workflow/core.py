#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  7 17:22:55 2020

@author: timothy
"""
from __future__ import print_function
import os
import nibabel
import sys
sys.path.append('../')
from tools.connectivity_metrics import NetworkIO, GraphTheoryIO
#from tools.timeseries
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
"""
parser = argparse.ArgumentParser(description='')
parser.add_argument('input_dir', help='The directory where the preprocessed derivative needed live')
parser.add_argument('output_dir', help='The directory where the output files should be stored.')
parser.add_argument('analysis_level', choices=['participant'],help='Processing stage to be run, only "participant" in this case (see BIDS-Apps specification).')
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
parser.add_argument('--apply_Fishers_r_to_z_transform', help="For correlation outputs, should Fisher's r-to-z transformation be applied? Choises are 'Yes' or 'No'.", choices = ['YES','NO'],default='YES')
parser.add_argument('--network_matrix_calculation', help="What method to employ for network matrix estimation. "
                                                    " Choices are 'All', 'correlation','partial_correlation', "
                                                    " 'dynamic_time_warping', 'tangent', 'covariance', 'sparse_inverse_covariance', "
                                                    "  'precision', 'sparse_inverse_precision'. NOTE: Specifying sparse matrices or dynamic time warping will result in increased computation time. For more information on available methods https://nilearn.github.io/connectivity/index.html#functional-connectivity-and-resting-state", 
                                                    choices =['All','correlation','partial_correlation',
                                                              'dynamic_time_warping','covariance',
                                                              'precision','sparse_inverse_precision',
                                                              'sparse_inverse_covariance'], default='correlation',nargs='+')
parser.add_argument('--graph_theory',help="Whether or not to output graph theoretical measures. Choices are 'All', 'clustering_coefficient','local_efficiency','strength','node_betweenness_centrality', 'edge_betweenness_centrality', 'eigenvector_centrality', and 'NONE'. If there are multiple measures that are of interest but not all are, separate as many choices as interested with a space.", choices = ['All','clustering_coefficient','local_efficiency','strength',
'node_betweenness_centrality', 'edge_betweenness_centrality', 'eigenvector_centrality','NONE'], default = 'NONE',nargs='+')
parser.add_argument('--timeseries_processing',help="Modify timeseries prior to generating network matrices. Choices include: 'entropy','alff','falff','wavelet', 'NONE'.",
                    choices = ['entropy','alff','falff','wavelet','NONE'],default='NONE')
parser.add_argument('--num_cpus', help='How many concurrent CPUs to use',default=1)
args = parser.parse_args()
"""

class MACCHIATO_setup:
    def __init__(self,**kwargs):

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
        ICA_outputs : [Yes, yes, No, no]
            Use ICA (whether FIX or AROMA) outputs for network matrix estimation.
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
        graph_theory : 
            Whether or not to output graph theoretical measures
        timeseries_processing : string
            Whether or not to output wavelet matrices 


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
        self.input_dir = args[0]
        self.output_dir = args[1]
        self.analysis_level = args[2]
        self.preprocessing_type = args[3]
        self.ICA_outputs = args[4]
        self.combine_resting_scans = args[5]
        self.parcellation_file = args[6]
        self.parcellation_name = args[7]
        self.reg_name = args[8]
        self.apply_Fishers_r_to_z_transform = args[9]
        self.network_matrix_calculation = args[10]
        self.graph_theory = args[11]
        self.timeseries_processing = args[12]
        self.msm_all_reg_name="MSMAll_2_d40_WRN"
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
            read_parcel_file = cifti.read(self.parcel_file)
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
        
        with h5py.File(os.path.join(self.output_dir,'data.hdf5'),'w',driver='mpio',comm=comm) as hdf5:
            network_matrix_dset = hdf5.create_dataset(name=self.network_matrix_calculation, shape=(image_count,height,width), dtype='f')
            network_matrix_dset.attrs['parcel_file'] = self.parcel_file
            network_matrix_dset.attrs['parcel_name'] = self.parcel_name    
            for i in range(rank, image_count, comm.size):
                bold = self.bolds[i]
                if self.timeseries_processing == 'YES':
                    pass #TODO will be passed to timeseries_metrics.py
                else:
                    network_metric_init = NetworkIO(output_dir=self.output_dir, 
                                        cifti_file=bold, 
                                        parcel_file=self.parcel_file, 
                                        parcel_name=self.parcel_name, 
                                        network_metric=self.network_matrix_calculation,
                                        fishers_r_to_z_transform=self.fishers_r_to_z_transform)
                if self.graph_theory == 'NONE':
                    if type(self.network_matrix_calculation) == str:
                        metric_data = network_metric_init.create_network_matrix()
                        network_matrix_dset[i,:,:]=metric_data
                    else:
                        pass
                else:
                    GraphTheoryIO(output_dir=self.output_dir, 
                          cifti_file=bold, 
                          parcel_file=self.parcel_file, 
                          parcel_name=self.parcel_name,
                          network_metric=self.network_matrix_calculation,
                          fishers_r_to_z_transform=self.fishers_r_to_z_transform,
                          graph_theory_metric=self.graph_theory)
if __name__ == '__main__':
    pdb.set_trace()
    MACCHIATO_setup(sys.argv[1:])
