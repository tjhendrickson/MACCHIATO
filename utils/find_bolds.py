#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  7 16:05:44 2020

@author: timothy
"""

def create_bold_lists(layout,combine_resting_scans,
                      preprocessing_type,ICA_outputs,
                      msm_all_reg_name,selected_reg_name):
    if combine_resting_scans == 'NO':
        if preprocessing_type == 'HCP':
            # use ICA outputs
            if ICA_outputs == 'YES':
                ICA_string="_FIXclean"
                if selected_reg_name == msm_all_reg_name:
                    bolds = [f.filename for f in layout.get(type='clean',extensions="dtseries.nii",task='rest') if msm_all_reg_name+'_hp2000_clean' in f.filename]
                else:
                    bolds = [f.filename for f in layout.get(type='clean',extensions="dtseries.nii", task='rest') if '_hp2000_clean' and not msm_all_reg_name in f.filename]
            # do not use ICA outputs
            else:
                ICA_string=""
                if selected_reg_name == msm_all_reg_name:
                    bolds = [f.filename for f in layout.get(extensions="dtseries.nii", task='rest') if msm_all_reg_name + '_hp2000' in f.filename and not 'clean' in f.filename]
                else:
                    bolds = [f.filename for f in layout.get(extensions="dtseries.nii", task='rest') if '_hp2000' in f.filename and not 'clean' and not msm_all_reg_name in f.filename]
        elif preprocessing_type == 'fmriprep':
            #use ICA outputs
            if ICA_outputs == 'YES':
                ICA_string="_AROMAclean"
                bolds = [f.filename for f in layout.get(type='bold',task='rest') if 'smoothAROMAnonaggr' in f.filename]
            # do not use ICA outputs
            else:
                ICA_string=""
                bolds = [f.filename for f in layout.get(type='bold',task='rest') if 'preproc' in f.filename]
            bolds_ref = [f.filename for f in layout.get(type='boldref',task='rest')]
        return sorted(bolds)
    else:
        combined_bolds_list = []
        if layout.get_sessions() > 0:
            for scanning_session in layout.get_sessions():
                # retreive subject id that is associated with session id and parse data with subject and session id
                for subject in layout.get_subjects(session=scanning_session):
                    if preprocessing_type == 'HCP':
                        # use ICA outputs
                        if ICA_outputs == 'YES':
                            ICA_string="_FIXclean"
                            if selected_reg_name == msm_all_reg_name:
                                bolds = [f.filename for f in layout.get(type='clean',extensions="dtseries.nii",task='rest',subject=subject,session=scanning_session) if msm_all_reg_name+'_hp2000_clean' in f.filename]
                            else:
                                bolds = [f.filename for f in layout.get(type='clean',extensions="dtseries.nii", task='rest',subject=subject,session=scanning_session) if '_hp2000_clean' and not msm_all_reg_name in f.filename]
                        # do not use ICA outputs
                        else:
                            ICA_string=""
                            if selected_reg_name == msm_all_reg_name:
                                bolds = [f.filename for f in layout.get(extensions="dtseries.nii", task='rest',subject=subject,session=scanning_session) if msm_all_reg_name + '_hp2000' in f.filename and not 'clean' in f.filename]
                            else:
                                bolds = [f.filename for f in layout.get(extensions="dtseries.nii", task='rest',subject=subject,session=scanning_session) if '_hp2000' in f.filename and not 'clean' and not msm_all_reg_name in f.filename]
                    elif preprocessing_type == 'fmriprep':
                        #use ICA outputs
                        if ICA_outputs == 'YES':
                            ICA_string="_AROMAclean"
                            bolds = [f.filename for f in layout.get(type='bold',task='rest',subject=subject,session=scanning_session) if 'smoothAROMAnonaggr' in f.filename]
                        # do not use ICA outputs
                        else:
                            ICA_string=""
                            bolds = [f.filename for f in layout.get(type='bold',task='rest') if 'preproc' in f.filename]
                        bolds_ref = [f.filename for f in layout.get(type='boldref',task='rest')]
                    if len(bolds) == 2:
                        combined_bolds_list.append(bolds)
        else:
            for scanning_session in layout.get_subjects():
                if preprocessing_type == 'HCP':
                    # use ICA outputs
                    if ICA_outputs == 'YES':
                        ICA_string="_FIXclean"
                        if selected_reg_name == msm_all_reg_name:
                            bolds = [f.filename for f in layout.get(type='clean',extensions="dtseries.nii",task='rest',subject=scanning_session) if msm_all_reg_name+'_hp2000_clean' in f.filename]
                        else:
                            bolds = [f.filename for f in layout.get(type='clean',extensions="dtseries.nii", task='rest',subject=scanning_session) if '_hp2000_clean' and not msm_all_reg_name in f.filename]
                    # do not use ICA outputs
                    else:
                        ICA_string=""
                        if selected_reg_name == msm_all_reg_name:
                            bolds = [f.filename for f in layout.get(extensions="dtseries.nii", task='rest',subject=scanning_session) if msm_all_reg_name + '_hp2000' in f.filename and not 'clean' in f.filename]
                        else:
                            bolds = [f.filename for f in layout.get(extensions="dtseries.nii", task='rest',subject=scanning_session) if '_hp2000' in f.filename and not 'clean' and not msm_all_reg_name in f.filename]
                elif preprocessing_type == 'fmriprep':
                    #use ICA outputs
                    if ICA_outputs == 'YES':
                        ICA_string="_AROMAclean"
                        bolds = [f.filename for f in layout.get(type='bold',task='rest',subject=scanning_session) if 'smoothAROMAnonaggr' in f.filename]
                    # do not use ICA outputs
                    else:
                        ICA_string=""
                        bolds = [f.filename for f in layout.get(type='bold',task='rest') if 'preproc' in f.filename]
                    bolds_ref = [f.filename for f in layout.get(type='boldref',task='rest')]
                if len(bolds) == 2:
                    combined_bolds_list.append(bolds)
        return sorted(combined_bolds_list)