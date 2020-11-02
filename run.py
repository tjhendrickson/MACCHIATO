#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 23 16:49:56 2020

@author: hendr522
"""

import argparse
from argunparse  import ArgumentUnparser #https://github.com/mbdevpl/argunparse
import os
from workflow.core import MACCHIATO_setup
import pdb
import json
import sys

# specify arguments that MACCHIATO accepts
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

# unparse arguments to feed into main workflow
kwargs = vars(parser.parse_args())
unparser = ArgumentUnparser()
arg_string = unparser.unparse(**kwargs)
pdb.set_trace()

#mpiexec -n numprocs python -m mpi4py -m mod [arg] ..

if args.participant_label:
    if len(args.participant_label) > 1:
        os.system('./start_workflow_mpi.sh {cpus} {args}'.format(cpus=str(args.num_cpus),args=arg_string))
    else:
        print('You specified running an individual participant, but set the "--num_cpus" greater than 1. Reverting back to "--num_cpus"=1...')
        args.num_cpus=1
        os.system('./start_workflow_mpi.sh {cpus} {args}'.format(cpus=str(args.num_cpus),args=arg_string))
else:
    os.system('./start_workflow_mpi.sh {cpus} {args}'.format(cpus=str(args.num_cpus),args=arg_string))