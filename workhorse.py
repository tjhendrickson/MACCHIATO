#!/usr/bin/env python3

from __future__ import print_function
import argparse
import os
from functools import partial
from subprocess import Popen, PIPE
import subprocess
from workflow.utils import MACCHIATO_setup
from tools.core import execute_MACCHIATO_instances
import pdb
from mpi4py import MPI
import h5py
import numpy as np
import cifti
from pprint import pprint

# function which actually launches processes to underlying system
def run(command, env={}, cwd=None):
    merged_env = os.environ
    merged_env.update(env)
    merged_env.pop("DEBUG", None)
    print(command)
    process = Popen(command, stdout=PIPE, stderr=subprocess.STDOUT,
                    shell=True, env=merged_env, cwd=cwd)
    while True:
        line = process.stdout.readline()
        print(line)
        line = str(line)[:-1]
        if line == '' and process.poll() != None:
            break
    if process.returncode != 0:
        raise Exception("Non zero return code: %d"%process.returncode)

# specify arguments that MACCHIATO accepts
parser = argparse.ArgumentParser(description='')
parser.add_argument('args', help='The arguments passed from run.py',required=True)

args = parser.parse_args()


# now parse arguments and print to standard output (STDOUT), slight difference in call if done as batch or participant

if args.group == 'batch':
    setupParams = MACCHIATO_setup(group=args.group,
                    preprocessing_type=args.preprocessing_type,
                    ICA_outputs=args.use_ICA_outputs,
                    graph_theory=args.graph_theory,
                    network_matrix_calculation=args.network_matrix_calculation,
                    input_dir=args.input_dir,
                    parcel_name=args.parcellation_name,
                    parcel_file=args.parcellation_file,
                    fishers_r_to_z_transform=args.apply_Fishers_r_to_z_transform,
                    selected_reg_name=args.reg_name,
                    timeseries_processing=args.timeseries_processing,
                    combine_resting_scans=args.combine_resting_scans)
elif args.group == 'participant':
    setupParams = MACCHIATO_setup(group=args.group,
                    preprocessing_type=args.preprocessing_type,
                    ICA_outputs=args.use_ICA_outputs,
                    graph_theory=args.graph_theory,
                    network_matrix_calculation=args.network_matrix_calculation,
                    input_dir=args.input_dir,
                    parcel_name=args.parcellation_name,
                    fishers_r_to_z_transform=args.apply_Fishers_r_to_z_transform,
                    parcel_file=args.parcellation_file,
                    selected_reg_name=args.reg_name,
                    timeseries_processing=args.timeseries_processing,
                    combine_resting_scans=args.combine_resting_scans,
                    participant_label=args.participant_label,
                    session_label=args.session_label)

#transform "self" variables from MACCHIATO_setup to dictionary
setupParams = setupParams.__dict__

# set up multiprocessing/parallelization allocation

comm = MPI.COMM_WORLD
rank = comm.rank

# retreive number of, height and width of matrices 
image_count = len(setupParams['bolds'])
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
    for i in range(comm.rank, image_count, comm.size):
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
        
