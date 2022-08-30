"""
@author: mhaghigh
"""
import pandas as pd
import numpy as np
import seaborn as sns
from sklearn import preprocessing
from sklearn.cluster import KMeans
from ..preprocess import extract_cpfeature_names
from ..preprocess import handle_nans

################################################################################

def extract_single_cell_samples(df_p_s,n_cells,cell_selection_method):
    """ 
    This function select cells based on input cell_selection_method
  
    Inputs: 
    ++ df_p_s   (pandas df) size --> (number of single cells)x(columns): 
    input dataframe contains single cells profiles as rows 
    
    ++ n_cells (dtype: int): number of cells to extract from the input dataframe
    
    ++ cell_selection_method (str): 
        - 'random' - generate n randomly selected cells
        - 'representative' - clusters the data and sample from the "closest to mean cluster"
        - 'geometric_median' - plots single sample than is the geometric median of samples 
           -> method 1 (hdmedians package): faster but can handle up to 800 rows without memory error
           -> method 2 (skfda package): slower but can handle up to 1500 rows without memory error
    
    Returns: 
    dff (pandas df): sampled input dataframe
  
    """    


    cp_features, cp_features_analysis_0 =  extract_cpfeature_names.extract_cpfeature_names(df_p_s);
    df_p_s, cp_features_analysis = handle_nans.handle_nans(df_p_s,cp_features_analysis_0);
    
    
#     print("heloo")
#     print(cp_features)
    if cell_selection_method=='random':
        dff=df_p_s.reset_index(drop=True).sample(n = n_cells, replace = False).reset_index(drop=True)

    elif cell_selection_method=='representative': 
        df_p_s[cp_features_analysis] = df_p_s[cp_features_analysis].interpolate()
        if df_p_s.shape[0]>60:
            n_cells_in_each_cluster_unif=30
        else:
            n_cells_in_each_cluster_unif=int(df_p_s.shape[0]/5) 
            
        n_clusts=int(df_p_s.shape[0]/n_cells_in_each_cluster_unif) 
        kmeans = KMeans(n_clusters=n_clusts).fit(np.nan_to_num(df_p_s[cp_features_analysis].values))
        clusterLabels=kmeans.labels_
        df_p_s['clusterLabels']=clusterLabels;
        mean_clus=kmeans.predict(df_p_s[cp_features_analysis].mean().values[np.newaxis,])
        df_ps=df_p_s[df_p_s["clusterLabels"]==mean_clus[0]]
        dff=df_ps.reset_index(drop=True).sample(n = np.min([n_cells,df_ps.shape[0]]), replace = False).reset_index(drop=True)

    elif cell_selection_method=='geometric_median':    
        import hdmedians as hd
        from skfda import FDataGrid
        from skfda.exploratory.stats import geometric_median            
        
#     #     method 1
# #         ps_arr=df_p_s[cp_features_analysis].astype(np.float32).values
#         ps_arr=df_p_s[cp_features_analysis].values
#         gms=hd.medoid(ps_arr,axis=0)
#         gm_sample_ind=np.where(np.sum((ps_arr-gms),axis=1)==0)[0]
#         df_p_s_gm=df_p_s.loc[gm_sample_ind,:]
#         dff=pd.concat([df_p_s_gm,df_p_s_gm],ignore_index=True)

    #     method 2
        ps_arr=df_p_s[cp_features_analysis].values
        X = FDataGrid(ps_arr)
        gms2 = np.squeeze(geometric_median(X).data_matrix)
        # gm2_sample_ind=np.where(np.sum((ps_arr-gms2),axis=1)==0)[0]
        gm2_sample_ind=np.array([np.argmin(np.sum(abs(ps_arr-gms2),axis=1))])
        df_p_s_gm2=df_p_s.loc[gm2_sample_ind,:]
        dff=pd.concat([df_p_s_gm2,df_p_s_gm2],ignore_index=True)
        
    return dff,cp_features_analysis