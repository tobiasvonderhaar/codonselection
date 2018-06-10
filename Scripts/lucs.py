"""
Created on Sun Jan 14 19:39:43 2018

@author: Tobias von der Haar
"""

import numpy as np
import pandas as pd

def read_luc(filename):
    raw_data = pd.read_csv(filename,header=None)
    #read out the three tables corresponding to factors, FLuc data and RLuc data.
    #Flatten the two-dimensonal data into linear vectors.

    Factors = raw_data.iloc[2:10,1:13].values.flatten()
    FLuc = raw_data.iloc[13:21,1:13].values.flatten()
    RLuc = raw_data.iloc[24:32,1:13].values.flatten()

    #determine all rows where factor is NOT empty
    empty_index = Factors != 'empty'
    Factors = Factors[empty_index]
    FLuc = FLuc[empty_index]
    RLuc = RLuc[empty_index]

    #check if the Factors table contains more than one factor level
    #if so, split into two factor columns at ':'
    if ':' in Factors[0]:
        luc_df = pd.DataFrame([entries.split(':') for entries in Factors])
        luc_df.columns = ['Factor 1','Factor 2']
        luc_df['FLuc'] = FLuc
        luc_df['RLuc'] = RLuc
    else:
        Factors = pd.Series(Factors)
        luc_df = pd.concat([Factors,FLuc,RLuc], axis=1)
        luc_df.columns = ['Factor','FLuc','RLuc']

    #replace the FLuc and RLuc values with the F/R ratio
    luc_df['FLuc'] = pd.to_numeric(luc_df['FLuc'])
    luc_df['RLuc'] = pd.to_numeric(luc_df['RLuc'])
    luc_df['FR_ratio'] = luc_df.FLuc / luc_df.RLuc
    del luc_df['RLuc'], luc_df['FLuc']

    return luc_df

#==============================================================================
