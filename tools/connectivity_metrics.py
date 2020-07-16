#!/opt/Miniconda3/bin/python

import nibabel.cifti2
import pandas as pd
import os
import cifti
import pdb
import fcntl as F
import csv
import numpy as np
from nilearn.connectome import ConnectivityMeasure
from sklearn.covariance import GraphicalLassoCV,EmpiricalCovariance
import scipy.linalg as la
import bct

class NetworkIO:
    def __init__(self,output_dir,cifti_file, parcel_file, parcel_name,network_metric):
        
        # path that data will be written to
        self.output_dir = output_dir
        # inputted cifti file
        self.cifti_file = cifti_file
        # inputted atlas/parcellation file
        self.parcel_file = parcel_file
        # shorthand name chosen for parcel file
        self.parcel_name = parcel_name
        self.network_metric = network_metric
        
        # create output folder if it does not exist
        if not os.path.isdir(self.output_dir):
            os.makedirs(self.output_dir)
        # determine fmriname
        self.fmriname = os.path.basename(cifti_file).split('.')[0]
        
        try:
            read_parcel_file = cifti.read(self.parcel_file)
        except TypeError:
            print('This does not look like a CIFTI parcel file. Must exit')
        parcel_file_label_tuple = read_parcel_file[1][0][0][1]
        parcel_labels = []
        
        for value in parcel_file_label_tuple:
                if not '???' in parcel_file_label_tuple[value][0]:
                        parcel_labels.append(parcel_file_label_tuple[value][0])
        self.parcel_labels = [str(r) for r in parcel_labels]
        
        # perform tests on inputted cifti file and parcellate timeseries
        self.cifti_tests()
        
        # now generate connectivity matrix
        self.create_network_matrix()

    def cifti_tests(self):
        # does CIFTI file exist?
        try:
            read_cifti = open(self.cifti_file)
            read_cifti.close()
        except IOError:
            print("file does not exist")
            
        # is entered CIFTI file actually a CIFTI file?
        try:
            self.cifti_load = nibabel.cifti2.cifti2.load(self.cifti_file)
        except:
            print("file does not look like a cifti file")
        cifti_file_basename = os.path.basename(self.cifti_file)
        cifti_prefix = cifti_file_basename.split(".")[0]
        cifti_suffix = '.'.join(cifti_file_basename.split(".")[1:])
        if cifti_suffix == 'dtseries.nii':
            self.new_cifti_suffix = '.ptseries.nii'
        # elif cifti_suffix == 'dscalar.nii':
        #    self.new_cifti_suffix = '.pscalar.nii'
        os.system("/opt/workbench/bin_rh_linux64/wb_command -cifti-parcellate %s %s %s %s" 
                  % (self.cifti_file, 
                     self.parcel_file, 
                     "COLUMN", 
                     os.path.join(self.output_dir,cifti_prefix) + "_"+self.parcel_name + self.new_cifti_suffix))
        parcellated_cifti_file = os.path.join(self.output_dir,cifti_prefix) + "_"+self.parcel_name + self.new_cifti_suffix
        self.parcellated_cifti_file = parcellated_cifti_file
        # does CIFTI file exist?
        try:
            read_cifti = open(self.parcellated_cifti_file)
            read_cifti.close()
        except IOError:
            print("file does not exist")    
        # is entered CIFTI file actually a CIFTI file?
        try:
            self.parcellated_cifti_load = nibabel.cifti2.cifti2.load(self.parcellated_cifti_file)
        except:
            print("file does not look like a cifti file")

    def create_network_matrix(self):
        print('\n')
        print('rsfMRI_network_metrics.py: Create network matrix ')
        print('\t-Output folder: ' + self.output_dir)
        print('\t-Cifti file: ' + self.cifti_file)
        print('\t-Parcel file: ' + self.parcel_file)
        print('\t-Network matrix method/type: ' + str(self.network_metric))
        cifti_np_array = np.array(self.parcellated_cifti_load.get_fdata())
        if self.network_metric == 'correlation':
            measure = ConnectivityMeasure(kind='correlation',cov_estimator=EmpiricalCovariance())
        elif self.network_metric == 'covariance':
            measure = ConnectivityMeasure(kind='covariance',cov_estimator=EmpiricalCovariance())
        elif self.network_metric == 'partial_correlation':
            measure = ConnectivityMeasure(kind='partial correlation',cov_estimator=EmpiricalCovariance())
        elif self.network_metric == 'precision':
            measure = ConnectivityMeasure(kind='precision',cov_estimator=EmpiricalCovariance())
        elif 'sparse' in self.network_metric:
            measure = GraphicalLassoCV()
            
        if 'sparse' in self.network_metric:
            measure.fit(cifti_np_array)
            if 'covariance' in self.network_metric:
                network_matrix = measure.covariance_
            elif 'precision' in self.network_metric:
                network_matrix = measure.precision_
        else:
            network_matrix = measure.fit_transform([cifti_np_array])[0]
        self.network_matrix = network_matrix
        # currently not vectorizing outputs, could come back to this later
        #vectorized_network_matrix = network_matrix[np.triu_indices(network_matrix.shape[0],k=1)]
        #self.vectorized_network_matrix = vectorized_network_matrix
class GraphTheoryIO(NetworkIO):
    def __init__(self,output_dir,cifti_file, parcel_file, parcel_name, network_metric, graph_theory_metric):
        #local efficiency, strength, node betweeness centrality, edge betweenness centrality local density, eigenvector centrality, clustering coefficient
        # ensure that data does not have NaNs or Infs, place diagnonal to zero, and restrict floating point to five decimals for stability

        # execute parent class NetworkIO to produce connectivity matrix
        super().__init__(output_dir,cifti_file, parcel_file, parcel_name,network_metric)
        
        self.graph_theory_metric = graph_theory_metric

        self.network_matrix[np.where(np.isinf(self.network_matrix))] = 0
        self.network_matrix[np.where(np.isnan(self.network_matrix))] = 0
        self.network_matrix[np.triu_indices(self.network_matrix.shape[0],k=1)]
        self.network_matrix = np.around(self.network_matrix, decimals=5)
        
        # separate metrics based on whether they include distances
        if graph_theory_metric == 'node_betweenness_centrality' or graph_theory_metric == 'edge_betweenness_centrality':
            # convert matrix to lengths
            lengths_network_matrix = bct.weight_conversion(self.network_matrix,'lengths')
            
            if graph_theory_metric == 'node_betweenness_centrality':
                # node betweeness centrality
                self.graph_theory_network_matrix = bct.centrality.betweenness_wei(lengths_network_matrix)
            elif graph_theory_metric == 'edge_betweenness_centrality':
                # edge betweeness centrality
                self.graph_theory_network_matrix = bct.centrality.edge_betweenness_wei(lengths_network_matrix)
        else:
            if graph_theory_metric == 'eigenvector_centrality':
                # eigenvector centrality
                self.graph_theory_network_matrix = bct.eigenvector_centrality_und(self.network_matrix)
            elif graph_theory_metric == 'local_efficiency':
                #local efficiency
                self.graph_theory_network_matrix = bct.efficiency_wei(self.network_matrix,local=True)
            elif graph_theory_metric == 'strength':        
                #node strength
                self.graph_theory_network_matrix = bct.strengths_und(self.network_matrix)
            elif graph_theory_metric ==  'clustering_coefficient':
                #clustering coefficent
                self.graph_theory_network_matrix = bct.clustering_coef_wu(self.network_matrix)








