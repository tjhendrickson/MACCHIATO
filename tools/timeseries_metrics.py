#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  7 17:56:52 2020

@author: timothy
"""

from scipy.signal import butter,filtfilt,periodogram
from entropy import sample_entropy
import pywt
import nibabel
import numpy as np
from subprocess import Popen, PIPE
import subprocess
# alff, entropy, wavelet,

class TimeSeriesIO:
    def __init__(self,output_dir,cifti_file,time_series_metric):
        '''
        Class initialization for timeseries processing. The primary purpose of this class is to:
        1) Performs tests on inputted arguments are in the expected format
        2) Generates new timeseries based on specified method
        '''
        self.output_dir = output_dir
        self.cifti_file = cifti_file
        self.time_series_metric = time_series_metric
        
        # perform tests on inputted cifti file
        self.cifti_tests()
        
        # process timeseries
        self.create_new_timeseries()

    def cifti_tests(self):
        # does CIFTI file exist?
        try:
            read_cifti = open(self.cifti_file)
            read_cifti.close()
        except IOError:
            print("file does not exist")
            
        # is entered CIFTI file actually a CIFTI file?
        try:
            self.cifti_load_np_array = np.array(nibabel.cifti2.cifti2.load(self.cifti_file).get_fdata())
        except:
            print("file does not look like a cifti file")
            
    def create_new_timeseries(self):
        # retreive TR
        command = '/opt/workbench/bin_rh_linux64/wb_command -file-information %s' %self.cifti_file
        process = Popen(command, stdout=PIPE, stderr=subprocess.STDOUT,
                    shell=True)
        output = process.stdout.readlines()
        for out in output:
            out = out.decode('utf-8')
            if 'Map Interval Step:' in out:
                TR = float(out.split(' ')[-1].strip())
        tps = self.cifti_load_np_array.shape[0]
        if self.time_series_metric == 'wavelet':
            #Which includes functional connectivity range scales: Scale 1 (0.23-0.45 Hz), Scale 2 (0.11-0.23 Hz), Scale 3 (0.06-0.11 Hz), Scale 4 (0.03-0.06 Hz), Scale 5 (0.01-0.03 Hz), Scale 6 (0.007-0.01 Hz).
            pass
        elif self.time_series_metric == 'entropy':
            # sample entropy
            pass
        elif self.time_series_metric == 'alff' or self.timeseries == 'falff':
            # alff and f/alff
            bf_b, bf_a = butter(N=2,Wn=np.array([0.01,0.08])/(1/TR/2),btype="bandpass") # 2nd order butterworth filter with band pass 0.01-0.08 Hz
            ar=np.array([1,0.5])
            ma=np.array([1])
            bp_filtered_cifti_timeseries = filtfilt(b=bf_b,a=bf_a,x=self.cifti_load_np_array)
            power_spectra_bp_filtered_cifti_timeseries_f, power_spectra_bp_filtered_cifti_timeseries_pxx = periodogram(x=bp_filtered_cifti_timeseries,detrend='constant',fs=TR)
            hertz=(1/TR/2)*power_spectra_bp_filtered_cifti_timeseries_f
            pass
        
        

