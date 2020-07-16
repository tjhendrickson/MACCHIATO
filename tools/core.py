#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  7 17:22:55 2020

@author: timothy
"""

import os
import nibabel
from .connectivity_metrics import NetworkIO, GraphTheoryIO
#from timeseries
from workflow.utils import create_file_output
import pdb


#handles all logic of inputted data to run MACCHIATO        
class execute_MACCHIATO_instances:
    def __init__(self,bold,preprocessing_type,parcel_file,parcel_name,
                  selected_reg_name,ICA_outputs,
                  combine_resting_scans,wavelet,network_metric,graph_theory_metric,
                  output_dir):
        pdb.set_trace()
        if wavelet == 'YES':
            pass #TODO will be passed to timeseries_metrics.py
        else:
            if graph_theory_metric == 'NONE':
                if type(network_metric) == str:
                    NetworkIO(output_dir=output_dir, cifti_file=bold, parcel_file=parcel_file, parcel_name=parcel_name, network_metric=network_metric)
                else:
                    pass
            else:
                GraphTheoryIO(output_dir=output_dir, cifti_file=bold, parcel_file=parcel_file, parcel_name=parcel_name,network_metric=network_metric,graph_theory_metric=graph_theory_metric)
        
