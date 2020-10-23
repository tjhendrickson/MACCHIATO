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
import pdb


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
        
