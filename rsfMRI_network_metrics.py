#!/opt/Anaconda3/bin/python

import nibabel.cifti2
import pandas as pd
import os
import cifti
import pdb
import fcntl as F
import csv
import numpy as np
from nilearn.connectome import ConnectivityMeasure
from sklearn.covariance import GraphicalLassoCV
import bct

class NetworkIO:
    def __init__(self,output_dir,cifti_file, parcel_file, parcel_name, method):
        # path that data will be written to
        self.output_dir = output_dir
        # inputted cifti file
        self.cifti_file = cifti_file
        # inputted atlas/parcellation file
        self.parcel_file = parcel_file
        # shorthand name chosen for parcel file
        self.parcel_name = parcel_name
        # method of network matrix to be outputted (i.e. correlation, covariance, partial correlation, etc)
        self.method = method
        
        # create output folder if it does not exist
        if not os.path.isdir(self.output_dir):
            os.makedirs(self.output_dir)
        # determine fmriname
        self.fmriname = os.path.basename(cifti_file).split('.')[0]
        
        #TODO: as the network matrices will be vectorized, a vector label will have to be created.
        # read parcel labels into list to query later
        read_parcel_file = cifti.read(self.parcel_file)
        parcel_file_label_tuple = read_parcel_file[1][0][0][1]
        parcel_labels = []
        
        for value in parcel_file_label_tuple:
                if not '???' in parcel_file_label_tuple[value][0]:
                        parcel_labels.append(parcel_file_label_tuple[value][0])
        self.parcel_labels = [str(r) for r in parcel_labels]

    def create_network_matrix(self):
        print('\n')
        print('rsfMRI_network_metrics.py: Create network matrix ')
        print('\t-Output folder: ' + self.output_dir)
        print('\t-Cifti file: ' + self.cifti_file)
        print('\t-Parcel file: ' + self.parcel_file)
        print('\t-Network matrix method/type: ' + str(self.method))
        
        # does CIFTI file exist?
        try:
            read_cifti = open(self.cifti_file)
            read_cifti.close()
        except IOError:
            print("file does not exist")
            
        # is entered CIFTI file actually a CIFTI file?
        try:
            cifti_load = nibabel.cifti2.cifti2.load(self.cifti_file)
        except:
            print("file does not look like a cifti file")
        # if file is actually a CIFTI file, parcellate it
        cifti_file_basename = os.path.basename(self.cifti_file)
        cifti_prefix = cifti_file_basename.split(".dtseries.nii")[0]
        os.system("/opt/workbench/bin_rh_linux64/wb_command -cifti-parcellate %s %s %s %s" 
                  % (self.cifti_file, 
                     self.parcel_file, 
                     "COLUMN", 
                     os.path.join(self.output_dir,cifti_prefix) + "_"+self.parcel_name + ".ptseries.nii"))
        cifti_file = os.path.join(self.output_dir,cifti_prefix) + "_"+self.parcel_name +".ptseries.nii"
        # now check if new parcellated CIFTI actually exists and is a CIFTI file.
        try:
            read_cifti = open(cifti_file)
            read_cifti.close()
        except IOError:
            print("file does not exist")    
        try:
            cifti_load = nibabel.cifti2.cifti2.load(cifti_file)
            cifti_np_array = np.array(cifti_load.get_fdata())
        except:
            print("file does not look like a cifti file")
        if self.method == 'correlation':
            #Pearson correlation coefficients with LedoitWolf covariance estimator
            measure = ConnectivityMeasure(kind='correlation')
        elif self.method == 'covariance':
            #LedoitWolf estimator
            measure = ConnectivityMeasure(kind='covariance')
        elif self.method == 'partial-correlation':
            # Partial correlation with LedoitWolf covariance estimator
            measure = ConnectivityMeasure(kind='partial correlation')
        elif self.method == 'tangent':
            measure = ConnectivityMeasure(kind='tangent')
        elif self.method == 'precision':
            measure = ConnectivityMeasure(kind='precision')
        elif 'sparse' in self.method:
            measure = GraphicalLassoCV()
            
        if 'sparse' in self.method:
            measure.fit(cifti_np_array)
            if 'covariance' in self.method:
                network_matrix = measure.covariance_
            elif 'precision' in self.method:
                network_matrix = measure.precision_
        else:
            network_matrix = measure.fit_transform([cifti_np_array])[0]
        vectorized_network_matrix = network_matrix[np.triu_indices(network_matrix.shape[0],k=1)]

        self.network_matrix = network_matrix
        self.vectorized_network_matrix = vectorized_network_matrix
    
    def create_graph_theory_outputs(self,graph_theory_outputs):
        #local efficiency, strength, average path distance, betweeness centrality, local density, degree centrality, cost centrality, eigenvector centrality, clustering coefficient
        
        #betweeness centrality
        #bct.centrality.betweenness_wei(self.network_matrix) # if weighted connection matrix
        #bct.centrality.betweenness_bin(self.network_matrix) #  if binary connection matrix 
        
        
        pass
    
    def create_text_output(self,ICAstring,text_output_dir,text_output_format,level):
        print('\n')
        print('rsfMRI_network_metrics.py: Create Text Output ')
        print('\t-Text output folder: %s' %str(text_output_dir))
        print('\t-Text output format: %s'%str(text_output_format))
        print('\t-Cifti file: %s' %str(self.cifti_file))
        print('\t-Parcel file: %s' %str(self.parcel_file))
        print('\t-Parcel name: %s' %str(self.parcel_name))
        print('\t-Network matrix output method/type: %s' %str(self.method))
        print('\t-The fmri file name: %s' %str(self.fmriname))
        print('\t-ICA String to be used to find FEAT dir, if any: %s' %str(ICAstring))
        print('\t-Analysis level to output data from: %s' %str(level))
        print('\n')
        # find vectorized network matrix

        # set outputs suffixes
        if text_output_format == "CSV" or text_output_format == "csv":
            # if file exists and subject and session have yet to be added, add to file
            output_text_file = os.path.join(text_output_dir,"_".join(self.fmriname.split('_')[2:])+"_"+self.parcel_name+ICAstring+'_level'+ str(level)+'_method'+self.method+".csv")
            # prevent race condition by using "try" statement
            try:
                read_output_text_file = open(output_text_file,'r')
                read_output_text_file.close()
            except:
                # if doesn't exist create headers and add subject/session data to file
                write_output_text_file = open(output_text_file,'w')
            # file exists and is accessible, ensure that to be appended data does not yet exist on it
            fieldnames = self.parcel_labels
            # append subject and session ID to fieldname list
            if os.path.basename(self.output_dir).split('-')[0] == 'ses':
                fieldnames.insert(0,'Session ID')
            fieldnames.insert(0,'Subject ID')
            
            # if dataset is empty pandas will throw an error
            try:     
                output_text_file_df = pd.read_csv(output_text_file)
            except:
                #add header and append data to file
                write_output_text_file = open(output_text_file,'w')
                try:
                    # try holding an exclusive lock first
                    F.flock(write_output_text_file, F.LOCK_EX | F.LOCK_NB)
                except IOError:
                    raise IOError('flock() failed to hold an exclusive lock')
                writer = csv.writer(write_output_text_file)
                writer.writerow(fieldnames)
                row_data = np.squeeze(zstat_data_img.get_fdata()).tolist()
                if os.path.basename(self.output_dir).split('-')[0] == 'ses':
                    row_data.insert(0,os.path.basename(self.output_dir).split('-')[1])
                    row_data.insert(0,self.output_dir.split('sub-')[1].split('/')[0])
                else:
                    row_data.insert(0,os.path.basename(self.output_dir).split('-')[1])
                pdb.set_trace()
                writer.writerow(row_data)
                # Unlock file
                try:
                    F.flock(write_output_text_file, F.LOCK_UN)
                    write_output_text_file.close()
                except IOError:
                    raise IOError('flock() failed to unlock file.')
            # find participant if it exists
            pdb.set_trace()
            row_data = np.squeeze(zstat_data_img.get_fdata()).tolist()
            if os.path.basename(self.output_dir).split('-')[0] == 'ses':
                session_id = os.path.basename(self.output_dir).split('-')[1]
                row_data.insert(0,session_id)
                subject_id = self.output_dir.split('sub-')[1].split('/')[0]    
                row_data.insert(0,subject_id)
            else:
                subject_id = os.path.basename(self.output_dir).split('-')[1]
                row_data.insert(0,subject_id)









