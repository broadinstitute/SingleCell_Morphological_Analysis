import pickle
import numpy as np
import pandas as pd 
import time
import sys, os
# from utils import read_data, visualize_data
from utils.read_data import *
import pandas as pd
from sqlalchemy import create_engine
from functools import reduce
import time

# from funcs.utils import readSingleCellData_sqlalch_well_subset,readSingleCellData_sqlalch_random_image_subset

drug_list_rank=pd.read_excel("/home/ubuntu/bucket/projects/2016_08_01_RadialMitochondriaDistribution_donna/\
workspace/Metadata_drugRep/drugList_20210115_uncorrForSiteAgg.xlsx")
drug_list_rank=drug_list_rank[~drug_list_rank["phenotype_abundance_pval"].isnull()]
drug_list_rank_Top=drug_list_rank.loc[0:20].append(drug_list_rank.loc[drug_list_rank.shape[0]-20:]).reset_index(drop=True)
# drug_list_rank=drug_list_rank.sort_values(by=['phenotype_abundance_t'],ascending=1)

meta_lincs=pd.read_csv("/home/ubuntu/bucket/projects/2018_04_20_Rosetta/\
workspace/results/synth_meta/meta_lincs_repLevel.csv")
meta_lincs.Metadata_mmoles_per_liter=meta_lincs.Metadata_mmoles_per_liter.values.round(2)
meta_lincs2=meta_lincs.groupby(['Metadata_broad_sample','Metadata_mmoles_per_liter','Metadata_Plate','Metadata_Well']).size().reset_index()


rootDirDrug='/home/ubuntu/bucket/projects/2015_10_05_DrugRepurposing_AravindSubramanian_GolubLab_Broad/workspace'
batchName='2016_04_01_a549_48hr_batch1_Mito_Project'

# cp_features=mergProf_treatLevel_lincs.columns[mergProf_treatLevel_lincs.columns.str.contains("Cells_|Cytoplasm_|Nuclei_")].tolist()
              
# pert_plate_well=mergProf_treatLevel_lincs.groupby(['PERT','Metadata_Plate','Metadata_Well']).size().reset_index()
# PERTS=pert_plate_well.PERT.unique().tolist()



for row in range(34,drug_list_rank_Top.shape[0]):#drug_list_rank_Top.shape[0]):#[1140:]:
    X_Metadata_broad_sample,X_Metadata_mmoles_per_liter=\
    drug_list_rank_Top.loc[row,['Metadata_broad_sample','Metadata_mmoles_per_liter']].values;
    
    
    pert=X_Metadata_broad_sample+"_"+str(X_Metadata_mmoles_per_liter)
    print(pert)
    pert_df=meta_lincs2[(meta_lincs2["Metadata_broad_sample"]==X_Metadata_broad_sample) &\
           (meta_lincs2["Metadata_mmoles_per_liter"]==X_Metadata_mmoles_per_liter)].sample(5).reset_index(drop=True)
    start_time = time.time()
    
    pert_df_allP0=[]
    pert_df_dmsos0=[]
    for j in range(pert_df.shape[0]):
        # p,wells="SQ00015195",["A13"]
        p,wells=pert_df.loc[j,"Metadata_Plate"],[pert_df.loc[j,"Metadata_Well"]]
        fileName=rootDirDrug+"/backend/"+batchName+"/"+p+"/"+p+".sqlite"
        
        print(p,wells)
        # df_p_s=readSingleCellData_sqlalch_random_image_subset(fileName,50);
        df_p_s=readSingleCellData_sqlalch_well_subset(fileName,wells);
        pert_df_allP0.append(df_p_s)
        
        ## fo DMSO
        wells=meta_lincs2[(meta_lincs2["Metadata_Plate"]==p) & (meta_lincs2["Metadata_broad_sample"]=="DMSO")].Metadata_Well.tolist()
        df_p_s0=readSingleCellData_sqlalch_well_subset(fileName,wells);
        pert_df_dmsos0.append(df_p_s0.sample(500))
    
    pert_df_allP = pd.concat(pert_df_allP0)
    pert_df_dmsos = pd.concat(pert_df_dmsos0)
    
    perWellData={}
    perWellData['cp_ss']=pert_df_allP
    perWellData['cp_ss_control']=pert_df_dmsos
    
    a_file = open("/home/ubuntu/bucket/projects/2016_08_01_RadialMitochondriaDistribution_donna/workspace/drugSCprofiles/"+pert+".pkl", "wb")
    pickle.dump(perWellData, a_file)
    a_file.close()
#     perWellData.to_pickle('/home/ubuntu/bucket/projects/2018_04_20_Rosetta/workspace/synth_l1k_ssCP_meta/'+pert);
    print("--- %s minutes ---" % ((time.time() - start_time)/60))