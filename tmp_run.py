import pandas as pd
import numpy as np
import seaborn as sns
sns.set(color_codes=True)
from sklearn import preprocessing
import matplotlib.pyplot as plt
# from utils import read_data, visualize_data
from utils.read_data import *
from utils.visualize_data import *
from sklearn.cluster import KMeans
import pooled_cell_painting_single_cell_visualization
# import times

import glob
rootDir='/home/ubuntu/calbucket/projects/2018_11_20_Periscope_Calico/'

# batch='20200805_A549_WG_Screen';
batch='20210422_6W_CP257';


batch_multi_name_dict={'20210422_6W_CP257':'CP257-HeLa-WG',\
                   '20200805_A549_WG_Screen':'CP186-A549-WG'}
batchName2=batch_multi_name_dict[batch]
sc_files_dir=rootDir+'workspace/software/'+batchName2+'/data/1.profiles/'+batch+'/single_cell/single_cell_by_guide/'
im_size=5500 # hardcoded for now, TODO: create a dictionary if this number is different for 257 vs 186

### metadata 
metadata_dir=rootDir+'workspace/metadata/'+batch+'/'
metadata_orig= pd.read_csv(metadata_dir+'Barcodes.csv')

# input_gene='KRT28'
box_size=100
# im_size=5500
n_cells=0
cell_selection_method='geometric_median'
channels=['DNA','Mito','Phalloidin','WGA','ER','Outline']

genes_ls=metadata_orig.gene_symbol.unique().tolist()

for igi in range(164,len(genes_ls)):
    input_gene=genes_ls[igi]
    all_guides_gms_ls=[]
    gene_guids_ls=glob.glob(sc_files_dir+'*_'+input_gene+'.csv.gz')
    for gi in gene_guids_ls:
        df_p_s=pd.read_csv(gi);
        
        for ch in channels:
            df_p_s["PathName_Corr"+ch]=rootDir+batch+'/images_corrected_cropped/'+df_p_s["Metadata_Foci_plate"]+'_'+df_p_s["Metadata_Foci_well"]+'/Corr'+ch
            df_p_s["FileName_Corr"+ch]="Corr"+ch+"_"+"Site_"+df_p_s["Metadata_Foci_site_location"].astype(str)+".tiff"

        df_p_s["Path_Outlines"]=rootDir+'workspace/analysis/'+batch+'/'+df_p_s["Metadata_Foci_plate"]+'-'+df_p_s["Metadata_Foci_well"]+'-'+df_p_s["Metadata_Foci_site_location"].astype(str)+'/'\
        +'/CorrDNA_Site_'+df_p_s["Metadata_Foci_site_location"].astype(str)+'_Overlay.png'

        df_p_s["Nuclei_Location_Center_X"]=df_p_s["Cells_AreaShape_Center_X"];
        df_p_s["Nuclei_Location_Center_Y"]=df_p_s["Cells_AreaShape_Center_Y"];
        
        df_p_s=edgeCellFilter2(df_p_s,im_size,box_size/2);

        if df_p_s.shape[0]>0:
            df_samples,cp_features_analysis = extract_single_cell_samples(df_p_s.sample(n = np.min([1500,df_p_s.shape[0]]),\
                            replace = False).reset_index(drop=True),n_cells,cell_selection_method);
            all_guides_gms_ls.append(df_samples)

    all_guides_gms_df=pd.concat(all_guides_gms_ls,ignore_index=True).drop_duplicates(ignore_index=True)
    
    if all_guides_gms_df.shape[0]==1:
        all_guides_gms_df=pd.concat([all_guides_gms_df]*2, ignore_index=True)    
    
    fig=visualize_n_SingleCell_pooled(channels,all_guides_gms_df,box_size,im_size,title=input_gene+'_'+cell_selection_method);

    
    resultsDir='/home/ubuntu/bucket/projects/2018_11_20_Periscope_Calico/workspace/visualizations/'+batch+'/geometric_median_guide_level/'
    fig.savefig(resultsDir+input_gene+'.png')  
#     plt.ioff()
