"""
@author: mhaghigh
Current methods for detecting untransfected cells:

1. Based on a single cp intentensity feature (usually "Cells_Intensity_MeanIntensity_" of "DsRed" or "Protein" channel)
  - If untransfected wells exist, we define the threshold based on the feature distribution of that wells:
    - plate level threshold - this one for each input plate calculates the threshold by inspecting ditribution of control 
        wells and derving plate specific parameters

    - batch level threshold - this needs to be precomputed by an independent prior step
         * the precomputed threshold should be based on scaled or unscaled features and should be noted next to the precomputed
           numbers
2. Based on all the intensity features in a channel (usually "Cells_Intensity_MeanIntensity_" of "DsRed" or "Protein" channel)

"""

import pandas as pd
import numpy as np
import sklearn.preprocessing as sp
    
# ---------------------------------------------------------------------------------
def extract_singlecell_transfection_labels(df,transfection_params_dict):
    
    """ 
    This function performs a transfection detection method on the input dataframe
    
    Inputs: 
    ++ df (pandas df) size --> (number of single cells in a full plate)x(columns): 
    input dataframe contains single cells profiles as rows for a whole plate 
    
    ++ params (dict): a dictionary filled with transfection parameters
        keys:
            - Method:
                1.'single_intensity_feature_thrsh'
                
                    'intensity_feature_to_use': 'Cells_Intensity_MeanIntensity_DsRed'
                    'thresholding_method': 'precomputed_batch_specific_thrsh' or 'batch_plate_specific_thrsh'
                    'precomputed_params': [precomputed_bottom_thresh , precomputed_top_thresh , data_norm_used]
                
                2. 'multiple_intensity_features'
                     
    
    Returns: 
    ++ a vector with the size df.shape[0] containing transfection_labels:
    (1) transfected, (-1) uncertain gray area, (0) untransfected
  
    """      

    if transfection_params_dict['Method']=='single_intensity_feature_thrsh':
        return transfection_detection_by_single_feature(df,transfection_params_dict)
                
    elif transfection_params_dict['Method']=='multiple_intensity_features':
        return transfection_detection_by_single_feature(df,transfection_params_dict)
        
        


def transfection_detection_by_single_feature(df,params):
    
    """ 
    This function performs a single feature thresholding based transfection detection on the input dataframe
    
    Inputs: 
    ++ df (pandas df) size --> (number of single cells in a full plate)x(columns): 
    input dataframe contains single cells profiles as rows for a whole plate 
    
    ++ params (dict): a dictionary filled with transfection parameters
        keys:
            'Method':'single_intensity_feature_thrsh'   
            'intensity_feature_to_use': 'Cells_Intensity_MeanIntensity_DsRed'
            'thresholding_method': 'precomputed_batch_specific_thrsh' or 'batch_plate_specific_thrsh'
            'precomputed_params': [precomputed_bottom_thresh , precomputed_top_thresh , data_norm_used]
                                  [low_thrsh , high_thrsh]
                
        
            example inputs:
                - transfection_params_dict={'Method':'single_intensity_feature_thrsh',\
                              'intensity_feature_to_use':'Cells_Intensity_MeanIntensity_DsRed',\
                              'thresholding_method': 'precomputed_batch_specific_thrsh',\
                              'pre_detection_scaler':'MinMax'
                              'precomputed_params': [precomputed_bottom_thresh , precomputed_top_thresh]}
                              
                - transfection_params_dict={'Method':'single_intensity_feature_thrsh',\
                              'intensity_feature_to_use':'Cells_Intensity_MeanIntensity_DsRed',\
                              'thresholding_method': 'batch_plate_specific_thrsh'}                            
    
    Returns: 
    ++ a vector with the size df.shape[0]
        transfection_labels:
            - (1) transfected
            - (-1) uncertain gray area
            - (0) untransfected
  
    """      
    if params['pre_detection_scaler']:
        #         clip target feature to its .999 percentile
        qpi_up=df[params['intensity_feature_to_use']].quantile(0.999)
        qpi_low=df[params['intensity_feature_to_use']].quantile(0.01)

        df[params['intensity_feature_to_use']]=df[params['intensity_feature_to_use']].clip(qpi_low, qpi_up)  
        
        if params['pre_detection_scaler']=='MinMax':
            pre_detection_scaler = sp.MinMaxScaler(feature_range=(0,1))
            df[params['intensity_feature_to_use']]=pre_detection_scaler.fit_transform(\
                                                      df[params['intensity_feature_to_use']].values.reshape(-1, 1))
        elif params['pre_detection_scaler']=='Robust':
            pre_detection_scaler = sp.RobustScaler()
            df[params['intensity_feature_to_use']]=pre_detection_scaler.fit_transform(\
                                                      df[params['intensity_feature_to_use']].values.reshape(-1, 1))
            
        else:
            raise Exception('scaler is not among the list! please add it!')
            

    if params['thresholding_method'] =='precomputed_batch_specific_thrsh':
        
        precomputed_bottom_thresh , precomputed_top_thresh  = params['precomputed_params']
        
        if np.isnan(precomputed_bottom_thresh):
            precomputed_bottom_thresh = precomputed_top_thresh
                
        
    elif params['thresholding_method'] =='batch_plate_specific_thrsh':                
        df = df.assign(Metadata_transfection = df['pert_name'].str.contains('Control'))
        
        untransfected_feature_values=df[df['Metadata_transfection'] == True][params['intensity_feature_to_use']].values
        
        precomputed_bottom_thresh=np.percentile(untransfected_feature_values, 50);
        precomputed_top_thresh=np.percentile(untransfected_feature_values, 99);        
        
        
    
    feature_values = df[params['intensity_feature_to_use']].values.reshape(-1, 1)
    transfection_labels = [1 if feature_values[i] >= precomputed_top_thresh else\
                           (0 if feature_values[i] < precomputed_bottom_thresh else -1)\
                           for i in range(len(feature_values))]        

        
    return transfection_labels 
        
    
    
    
def transfection_detection_by_clustering(df,params):
    
    """ 
    This function performs a single feature thresholding based transfection detection on the input dataframe
    
    Inputs: 
    ++ df (pandas df) size --> (number of single cells in a full plate)x(columns): 
    input dataframe contains single cells profiles as rows for a whole plate 
    
    ++ params (dict): a dictionary filled with transfection parameters
        keys:
            'Method':'single_intensity_feature_thrsh'   
            'intensity_feature_to_use': 'Cells_Intensity_MeanIntensity_DsRed'
            'thresholding_method': 'precomputed_batch_specific_thrsh' or 'batch_plate_specific_thrsh'
            'precomputed_params': [precomputed_bottom_thresh , precomputed_top_thresh , data_norm_used]
                                  [low_thrsh , high_thrsh , scaled/raw]
                
        
            example inputs:
                - transfection_params_dict={'Method':'single_intensity_feature_thrsh',\
                              'intensity_feature_to_use':'Cells_Intensity_MeanIntensity_DsRed',\
                              'thresholding_method': 'precomputed_batch_specific_thrsh',\
                              'precomputed_params': [precomputed_bottom_thresh , precomputed_top_thresh , data_norm_used]}
                              
                - transfection_params_dict={'Method':'single_intensity_feature_thrsh',\
                              'intensity_feature_to_use':'Cells_Intensity_MeanIntensity_DsRed',\
                              'thresholding_method': 'batch_plate_specific_thrsh'}                            
    
    Returns: 
    ++ a vector with the size df.shape[0]
        transfection_labels:
            - (1) transfected
            - (-1) uncertain gray area
            - (0) untransfected
  
    """      
    

#     elif method=='wellClustering':
    perWellData_prot=perWellData.loc[:,df.columns.str.contains("_Protein")]
    perWellData_prot=perWellData_prot.fillna(perWellData_prot.median())
    nCellsPerW=df.shape[0];
#     if nCellsPerW>50;
#         nCls=10
#     else:
#         nCls=2      
#     nCls=int(perWellData.shape[0]/200)
    ncls=4;
#     clustering = KMeans(n_clusters=ncls).fit(perWellData_prot)
    clustering = SpectralClustering(n_clusters=ncls,affinity='nearest_neighbors',assign_labels="discretize").fit(perWellData_prot)
    clusterLabels=clustering.labels_#.reshape(1,preprocData.shape[0])    


    df.loc[:,'clsLabel0']=clusterLabels

    print(df[intFeatureToUse].min(),\
     df[intFeatureToUse].max())  
    wellStats0=df[['clsLabel0',intFeatureToUse]].\
    groupby(['clsLabel0']).describe().reset_index(drop=True)#.sort_values(by=['mean'])
    wellStats0.columns = wellStats0.columns.droplevel(0)

    sortedClusters=wellStats0.sort_values(by=['mean']).index;
    wellStats=wellStats0.sort_values(by=['mean']).reset_index(drop=True);
    wellStats['perc']=wellStats['count'].apply(lambda x: x/wellStats['count'].sum())



    print(wellStats)
    map4labelOrdering=dict()
    for c in range(ncls):
        map4labelOrdering[sortedClusters[c]]=c

    df['clsLabel']=df['clsLabel0'].map(map4labelOrdering)    
    print(df['Metadata_Efficiency'].unique()[0],wellStats.loc[3,'perc'])
#     perWellData['ClustersForGFP']=perWellData['clsLabel'].apply(lambda x: 1 if wellStats.iloc[x]['mean'] >= globalClusThrsh  else (0 if wellStats.iloc[x]['mean'] <= globalMedian else -1))



    wellEffAnot=df['Metadata_Efficiency'].unique()[0]

    if (wellEffAnot in ['med','high']) and wellStats.loc[3,'perc']<0.3:
        transfClusters=[2,3]

    elif  (wellEffAnot in ['low', 'x', 'very low','very low, 3, cells']) and wellStats.loc[3,'perc']>0.5:
        ncls=6;
    #     clustering = KMeans(n_clusters=ncls).fit(perWellData_prot)
        clustering = SpectralClustering(n_clusters=ncls,affinity='nearest_neighbors',assign_labels="discretize").fit(perWellData_prot)
        clusterLabels=clustering.labels_#.reshape(1,preprocData.shape[0])    


        df.loc[:,'clsLabel0']=clusterLabels
        wellStats0=df[['clsLabel0',intFeatureToUse]].groupby(['clsLabel0']).describe().reset_index(drop=True)#.sort_values(by=['mean'])
        wellStats0.columns = wellStats0.columns.droplevel(0)

        sortedClusters=wellStats0.sort_values(by=['mean']).index;
        wellStats=wellStats0.sort_values(by=['mean']).reset_index(drop=True);
        wellStats['perc']=wellStats['count'].apply(lambda x: x/wellStats['count'].sum())

        map4labelOrdering=dict()
        for c in range(ncls):
            map4labelOrdering[sortedClusters[c]]=c

        print(wellStats)
        df['clsLabel']=df['clsLabel0'].map(map4labelOrdering)    
        print('modified',df['Metadata_Efficiency'].unique()[0],wellStats.loc[ncls-1,'perc'])

    else:
        transfClusters=[3]
        df['ClustersForGFP']=df['clsLabel'].apply(lambda x: 1 if (x in transfClusters) else\
                                              (0 if (x in [1,2]) else -1))
        
#     return df

        
    return transfection_labels    
    
        