#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  7 17:22:55 2020

@author: timothy
"""

import os
import nibabel
from .connectivity_metrics import NetworkIO, graph_theory
#from timeseries
from workflow.utils import create_file_output
import pdb


#handles all logic of inputted data to run MACCHIATO        
class execute_MACCHIATO_instances:
    def __init__(self,bold,preprocessing_type,parcel_file,parcel_name,
                  selected_reg_name,motion_confounds_filename,ICA_outputs,
                  combine_resting_scans,wavelet,network_matrix_calculation,graph_theory,
                  output_dir):
        pass
    
