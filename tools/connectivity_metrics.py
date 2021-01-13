#!/opt/Miniconda3/bin/python
import bct
import cifti
import csv
import fcntl as F
import nibabel.cifti2 as ci
from nilearn.connectome import ConnectivityMeasure
import numpy as np
import os
import pandas as pd
import scipy.linalg as la
from sklearn.covariance import GraphicalLassoCV,EmpiricalCovariance

class NetworkIO:
    def __init__(self,output_dir,cifti_data, parcel_file, parcel_name,network_metric,fishers_r_to_z_transform):
        '''
        Parameters
        ----------
        output_dir : string
            Where to output generated data
        cifti_data : string or list
            Path to inputted cifti file that will be processed
        parcel_name : string
            Shorthand name of the CIFTI label file. 
        parcel_file : string
            The CIFTI label file to use or used to parcellate the brain.
        network_metric : [All, correlation, partial_correlation, dynamic_time_warping, covariance, precision, sparse_inverse_precision, sparse_inverse_covariance]
            What method to employ for network matrix estimation.
        fishers_r_to_z_transform : [YES, NO]
            Should Fishe4r's r-to-z transformation be applied?

        Returns
        -------
        None.

        '''

        # path that data will be written to
        self.output_dir = output_dir
        # inputted cifti file
        self.cifti_data = cifti_data
        # inputted atlas/parcellation file
        self.parcel_file = parcel_file
        # shorthand name chosen for parcel file
        self.parcel_name = parcel_name
        self.fishers_r_to_z_transform = fishers_r_to_z_transform
        self.network_metric = network_metric
        
        # create output folder if it does not exist
        if not os.path.isdir(self.output_dir):
            os.makedirs(self.output_dir)
        # determine fmriname
        self.fmriname = os.path.basename(cifti_data).split('.')[0]
        
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
        
        if type(self.cifti_data) == list:
            example_cifti_img = ci.load(self.cifti_data[0])
            self.fmri_data_np_arr = np.zeros((example_cifti_img.shape[0],example_cifti_img.shape[1],len(self.cifti_data)))
            for idx, fmritcs in self.cifti_data:
                self.cifti_file = fmritcs
                self.cifti_tests()
                normalized_data = ((self.parcellated_cifti_data - self.parcellated_cifti_data .mean())/self.parcellated_cifti_data .std())
                self.fmri_data_np_arr[:,:,idx] = normalized_data
            self.df_cifti_load = pd.DataFrame(self.fmri_data_np_arr.mean(axis=2))
            self.cifti_np_array = self.df_cifti_load.to_numpy()

        else:
            self.cifti_file = self.cifti_data
            self.cifti_tests() # perform tests on inputted cifti file and parcellate timeseries
            self.cifti_np_array = np.array(self.parcellated_cifti_data )
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
            self.cifti_load = ci.load(self.cifti_file)
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
        self.parcellated_cifti_file = os.path.join(self.output_dir,cifti_prefix) + "_"+self.parcel_name + self.new_cifti_suffix
        # does CIFTI file exist?
        try:
            read_cifti = open(self.parcellated_cifti_file)
            read_cifti.close()
        except IOError:
            print("file does not exist")    
        # is entered CIFTI file actually a CIFTI file?
        try:
            self.parcellated_cifti_data = ci.load(self.parcellated_cifti_file).get_fdata()
        except:
            print("file does not look like a cifti file")

    def create_network_matrix(self):
        print('\n')
        print('rsfMRI_network_metrics.py: Create network matrix ')
        print('\t-Output folder: ' + self.output_dir)
        print('\t-Cifti file: ' + self.cifti_file)
        print('\t-Parcel file: ' + self.parcel_file)
        print('\t-Network matrix method/type: ' + str(self.network_metric))
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
            measure.fit(self.cifti_np_array)
            if 'covariance' in self.network_metric:
                network_matrix = measure.covariance_
            elif 'precision' in self.network_metric:
                network_matrix = measure.precision_
        else:
            network_matrix = measure.fit_transform([self.cifti_np_array])[0]
        if self.network_metric == 'correlation' and self.fishers_r_to_z_transform == 'YES':
            self.network_matrix = 0.5*(np.log(1+network_matrix)-np.log(1-network_matrix))
        else:
            self.network_matrix = network_matrix
        return self.network_matrix
        
class GraphTheoryIO(NetworkIO):
    def __init__(self,output_dir,cifti_file, parcel_file, parcel_name, network_metric, fishers_r_to_z_transform, graph_theory_metric):
        '''
        Parameters
        ----------
        output_dir : TYPE
            DESCRIPTION.
        cifti_data : TYPE
            DESCRIPTION.
        parcel_file : TYPE
            DESCRIPTION.
        parcel_name : TYPE
            DESCRIPTION.
        network_metric : TYPE
            DESCRIPTION.
        fishers_r_to_z_transform : TYPE
            DESCRIPTION.
        graph_theory_metric : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        '''
        
        #local efficiency, strength, node betweeness centrality, edge betweenness centrality local density, eigenvector centrality, clustering coefficient
        # ensure that data does not have NaNs or Infs, place diagnonal to zero, and restrict floating point to five decimals for stability

        # execute parent class NetworkIO to produce connectivity matrix
        super().__init__(output_dir,cifti_data, parcel_file, parcel_name,network_metric,fishers_r_to_z_transform)
        pdb.set_trace()
        self.graph_theory_metric = graph_theory_metric

        self.network_matrix[np.where(np.isinf(self.network_matrix))] = 0
        self.network_matrix[np.where(np.isnan(self.network_matrix))] = 0
        self.network_matrix[np.triu_indices(self.network_matrix.shape[0],k=1)]
        self.network_matrix = np.around(self.network_matrix, decimals=5)
        
        # separate metrics based on whether they include distances
        if graph_theory_metric == 'node_betweenness_centrality' or graph_theory_metric == 'edge_betweenness_centrality':
            # convert matrix to lengths
            lengths_network_matrix = bct.weight_conversion(self.network_matrix,'lengths')
            num_nodes = self.network_matrix.shape[0]
            if graph_theory_metric == 'node_betweenness_centrality':
                # node betweeness centrality
                self.graph_theory_connection_matrix,self.graph_theory_vector = bct.centrality.betweenness_wei(lengths_network_matrix)
            elif graph_theory_metric == 'edge_betweenness_centrality':
                # edge betweeness centrality
                self.graph_theory_vector = bct.centrality.edge_betweenness_wei(lengths_network_matrix)
        else:
            if graph_theory_metric == 'eigenvector_centrality':
                # eigenvector centrality
                self.graph_theory_vector = bct.eigenvector_centrality_und(self.network_matrix)
            elif graph_theory_metric == 'local_efficiency':
                #local efficiency
                self.graph_theory_vector = bct.efficiency_wei(self.network_matrix,local=True)
            elif graph_theory_metric == 'strength':        
                #node strength
                self.graph_theory_vector = bct.strengths_und(self.network_matrix)
            elif graph_theory_metric ==  'clustering_coefficient':
                #clustering coefficent
                self.graph_theory_vector = bct.clustering_coef_wu(self.network_matrix)
        return self.graph_theory_vector # TODO will likely need to wrap some logic into this return








