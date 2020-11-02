#!/bin/bash

num_cpus=$1
shift
args=$@

echo mpirun -n ${num_cpus} python -m mpi4py -m workflow.core ${args}
mpirun -n ${num_cpus} python -m workflow.core ${args}