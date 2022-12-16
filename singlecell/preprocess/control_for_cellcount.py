"""
@author: mhaghigh
"""

import numpy as np
import pandas as pd
from sklearn import linear_model

def control_feature_y_for_variable_x(df,y_col_name,x_col_name,y_col_name_postfix):

    ransac = linear_model.RANSACRegressor(min_samples=0.9)

    x2=df[x_col_name].values[:,np.newaxis]
    y2=df[y_col_name].values[:,np.newaxis]

    lm=ransac.fit(x2, y2)  ### means that we have controled for aggregated per site features
    df[y_col_name+y_col_name_postfix]=y2 - lm.predict(x2)

    return df