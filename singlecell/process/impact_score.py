"""
@author: mhaghigh

* Impact score: indicates the difference between profiles of two wells

"""

import pandas as pd
import numpy as np


# def pearson_corr_coef_score():
    
    for i in range(len(frame3.columns)):    
    x, y = frame3.iloc[ :,i].values, control['CONTROL'].values
    nas = np.logical_or(x.isnan(), y.isnan())
    corr = sp.pearsonr(x[~nas], y[~nas])
    correlation.append(corr)
    
#     return 