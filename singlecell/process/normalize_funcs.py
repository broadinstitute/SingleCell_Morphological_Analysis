"""
@author: mhaghigh
"""
import pandas as pd
import sklearn.preprocessing as sp

def standardize_per_catX(df,column_name,cp_features):
# column_name='Metadata_Plate'
#     cp_features=df.columns[df.columns.str.contains("Cells_|Cytoplasm_|Nuclei_")]
    df_scaled_perPlate=df.copy()
    df_scaled_perPlate[cp_features]=\
    df[cp_features+[column_name]].groupby(column_name)\
    .transform(lambda x: (x - x.mean()) / x.std()).values
    return df_scaled_perPlate





def standardize_df_columns(sc_df,cp_features_analysis,scaling_method):
    '''
    scaling method can be any method supported by sklearn package, ex. Robust, Standard, ..
    
    '''
    
#     scaling_string = "scaler = sp."+scaling_method+"Scaler()"  
#     exec(scaling_string)
#     scaler = sp.RobustScaler()
#     scaler = sp.StandardScaler()


    if scaling_method=='MinMax':
        scaler = sp.MinMaxScaler(feature_range=(-1,1))

    elif scaling_method=='Robust':
        scaler = sp.RobustScaler()
        
    elif scaling_method=='Standard':
        scaler = sp.StandardScaler()
            


    sc_df_output = sc_df.copy()
    sc_df_output[cp_features_analysis]=scaler.fit_transform(sc_df.loc[:,cp_features_analysis].values)       
    
    return sc_df_output


def zscore_df_columns_by_control(sc_control_df, sc_df,cp_features_analysis,scaling_method):
    '''
    scaling method can be any method supported by sklearn package, ex. Robust, Standard, ..
    
    '''
    
#     scaling_string = "scaler = sp."+scaling_method+"Scaler()"  
#     exec(scaling_string)

    if scaling_method=='MinMax':
        scaler = sp.MinMaxScaler(feature_range=(0,1))

    elif scaling_method=='Robust':
        scaler = sp.RobustScaler()
        
    elif scaling_method=='Standard':
        scaler = sp.StandardScaler()
        
    sc_df_output = sc_df.copy()
    scaler.fit(sc_control_df.loc[:,cp_features_analysis])
    sc_df_output.loc[:,cp_features_analysis]=scaler.transform(sc_df.loc[:,cp_features_analysis])       
    
    return sc_df_output
