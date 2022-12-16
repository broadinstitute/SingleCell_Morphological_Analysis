"""
@author: mhaghigh
"""
import pandas as pd
import numpy as np

################################################################################
def handle_nans(df_input, cp_features, thrsh_null_ratio=0.05, thrsh_std=0.0001, fill_na_method=None):
    """
    from the all df_input columns extract cell painting measurments 
    the measurments that should be used for analysis
    
    Inputs:
    df_input: dataframes with all the annotations available in the raw data
    fill_na_method: 'interpolate' or 'median' or None
                    interpolate makes sense for single cell data of arrayed experiment since it fills NA values 
                    to the nearest cell values 
    Outputs: cp_features, cp_features_analysis
    
    """
    
#     cp_features=df_input.columns[df_input.columns.str.contains("Cells_|Cytoplasm_|Nuclei_")].tolist()

    print("cp_features:",len(cp_features))
    object_type_columns=df_input[cp_features].select_dtypes([object]).columns
    df_input[object_type_columns] = df_input[object_type_columns].apply(pd.to_numeric, errors = 'coerce')

#     df_input=df_input.replace([np.inf, -np.inf,'nan'], np.nan)
#     df_input[cp_features] = df_input[cp_features].apply(pd.to_numeric, errors='coerce')
#     df_input[cp_features] = df_input[cp_features].astype(float)
    

#     thrsh_null_ratio=0.05; thrsh_std=0.0001;
    cols2remove_manyNulls=[i for i in cp_features if (df_input[i].isnull().sum(axis=0)/df_input.shape[0])\
                  >thrsh_null_ratio]   
    cols2remove_lowVars = df_input[cp_features].std()[df_input[cp_features].std() < thrsh_std].index.tolist()

    cols2removeCP = cols2remove_manyNulls + cols2remove_lowVars
    print('cols2remove_manyNulls',cols2remove_manyNulls)
    
    print('cols2remove_lowVars',cols2remove_lowVars)

    cp_features_analysis = list(set(cp_features) - set(cols2removeCP))
    print("len cp_features_analysis/nan cols/low vars:",\
          len(cp_features_analysis),len(cols2remove_manyNulls),len(cols2remove_lowVars))
#     cp_features_analysis_filt1 = list(set(cp_features) - set(cols2removeCP))
#     df_numeric_columns = df_input.select_dtypes([np.number]).columns
#     cp_features_analysis = list(set(df_numeric_columns) & set(cp_features_analysis_filt1))
     
    df_p_s=df_input.drop(cols2removeCP, axis=1);
    
#     print(cols2removeCP)
    
#     df_p_s[cp_features_analysis] = df_p_s[cp_features_analysis].interpolate()  
    if fill_na_method=='median':
        df_p_s.loc[:,cp_features_analysis] = df_p_s.loc[:,cp_features_analysis].fillna(df_p_s[cp_features_analysis].median())
    elif fill_na_method=='interpolate':
        df_p_s.loc[:,cp_features_analysis] = df_p_s.loc[:,cp_features_analysis].interpolate()
    elif fill_na_method=='drop-rows':
        print('before dropping nan rows: ',df_p_s.shape)
#         print('1',df_p_s.shape,df_p_s.dropna(subset=cp_features_analysis).reset_index(drop=True).shape)
        df_p_s=df_p_s.dropna(subset=cp_features_analysis).reset_index(drop=True)
        print('after dropping nan rows: ',df_p_s.shape)
        
    elif fill_na_method=='interpolate_sim_col': #interpolate based on the columns with highest correlation
        print('Not implemented yet! Nothing got dropped! ')

    
#     row_has_NaN = df_p_s[cp_features_analysis].isnull().any(axis=1)
#     print(row_has_NaN)
#     print(df_p_s[cp_features_analysis].dropna().shape,df_p_s[cp_features_analysis].shape)
#     df_p_s[cp_features_analysis] = df_p_s[cp_features_analysis].dropna() 
#     dataframe.fillna(0)
    
    return df_p_s, cp_features_analysis

# ################################################################################
# def handle_nans(df_input, cp_features, thrsh_null_ratio=0.05, thrsh_std=0.0001, fill_na_method=None):
#     """
#     from the all df_input columns extract cell painting measurments 
#     the measurments that should be used for analysis
    
#     Inputs:
#     df_input: dataframes with all the annotations available in the raw data
#     fill_na_method: 'interpolate' or 'median' or None
#                     interpolate makes sense for single cell data of arrayed experiment since it fills NA values 
#                     to the nearest cell values 
#     Outputs: cp_features, cp_features_analysis
    
#     """
    
# #     cp_features=df_input.columns[df_input.columns.str.contains("Cells_|Cytoplasm_|Nuclei_")].tolist()

#     df_input=df_input.replace([np.inf, -np.inf,'nan'], np.nan)
# #     df_input[cp_features] = df_input[cp_features].apply(pd.to_numeric, errors='coerce')
# #     df_input[cp_features] = df_input[cp_features].astype(float)
    

# #     thrsh_null_ratio=0.05; thrsh_std=0.0001;
#     cols2remove_manyNulls=[i for i in cp_features if (df_input[i].isnull().sum(axis=0)/df_input.shape[0])\
#                   >thrsh_null_ratio]   
#     cols2remove_lowVars = df_input[cp_features].std()[df_input[cp_features].std() < thrsh_std].index.tolist()

#     cols2removeCP = cols2remove_manyNulls + cols2remove_lowVars
# #     print(cols2removeCP)


#     cp_features_analysis_filt1 = list(set(cp_features) - set(cols2removeCP))
    
#     df_numeric_columns = df_input.select_dtypes([np.number]).columns
#     cp_features_analysis = list(set(df_numeric_columns) & set(cp_features_analysis_filt1))

#     df_p_s=df_input.drop(cols2removeCP, axis=1);
    
# #     print(cols2removeCP)
    
# #     df_p_s[cp_features_analysis] = df_p_s[cp_features_analysis].interpolate()  
#     if fill_na_method=='median':
#         df_p_s.loc[:,cp_features_analysis] = df_p_s.loc[:,cp_features_analysis].fillna(df_p_s[cp_features_analysis].median())
#     elif fill_na_method=='interpolate':
#         df_p_s.loc[:,cp_features_analysis] = df_p_s.loc[:,cp_features_analysis].interpolate()
#     elif fill_na_method=='drop-rows':
#         print('1',df_p_s.shape,df_p_s.dropna(subset=cp_features_analysis).reset_index(drop=True).shape)
#         df_p_s=df_p_s.dropna(subset=cp_features_analysis).reset_index(drop=True)
#         print(df_p_s.shape)
# #     row_has_NaN = df_p_s[cp_features_analysis].isnull().any(axis=1)
# #     print(row_has_NaN)
# #     print(df_p_s[cp_features_analysis].dropna().shape,df_p_s[cp_features_analysis].shape)
# #     df_p_s[cp_features_analysis] = df_p_s[cp_features_analysis].dropna() 
# #     dataframe.fillna(0)
    
#     return df_p_s, cp_features_analysis
