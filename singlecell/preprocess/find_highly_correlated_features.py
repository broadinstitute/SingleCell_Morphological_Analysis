import numpy as np
import pandas as pd

def find_correlation(data, threshold=0.9, remove_negative=False):
    
    """
    Inputs
    data : pandas DataFrame
    threshold : float
        correlation threshold, will remove one of pairs of features with a
        correlation greater than this value.
    remove_negative: Boolean
        If true then features which are highly negatively correlated will
        also be returned for removal.
   
    Output
    to_drop(list): list of column names to be removed
    """
    corr_mat = data.corr()
    if remove_negative:
        corr_mat = corr_mat.abs()
        
    upper = corr_mat.where(np.triu(np.ones(corr_mat.shape), k=1).astype(np.bool))

    to_drop = [column for column in upper.columns if any(upper[column] > threshold)]
        
    return to_drop