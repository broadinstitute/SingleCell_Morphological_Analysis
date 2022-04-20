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
from multiprocessing import Pool

# from funcs.utils import readSingleCellData_sqlalch_well_subset,readSingleCellData_sqlalch_random_image_subset

import pickle
# from funcs.utils import readSingleCellData_sqlalch_well_subset,readSingleCellData_sqlalch_random_image_subset

drug_list_rank=pd.read_excel("/home/ubuntu/bucket/projects/2016_08_01_RadialMitochondriaDistribution_donna/\
workspace/Metadata_drugRep/drugList_20210115_uncorrForSiteAgg.xlsx")
drug_list_rank=drug_list_rank[~drug_list_rank["phenotype_abundance_pval"].isnull()]
drug_list_rank_Top=drug_list_rank.loc[0:40].append(drug_list_rank.loc[drug_list_rank.shape[0]-20:]).reset_index(drop=True)
# drug_list_rank_Top=drug_list_rank.loc[0:40].reset_index(drop=True)
# drug_list_rank=drug_list_rank.sort_values(by=['phenotype_abundance_t'],ascending=1)

# meta_lincs=pd.read_csv("/home/ubuntu/bucket/projects/2018_04_20_Rosetta/\
# workspace/results/synth_meta/meta_lincs_repLevel.csv")
# meta_lincs.Metadata_mmoles_per_liter=meta_lincs.Metadata_mmoles_per_liter.values.round(2)

meta_lincs=pd.read_csv("/home/ubuntu/bucket/projects/2018_04_20_Rosetta/\
workspace/results/synth_meta/matadata_lincs_2.csv")
meta_lincs['Metadata_mmoles_per_liter']=meta_lincs.mmoles_per_liter.values.round(2)
meta_lincs=meta_lincs.rename(columns={"Assay_Plate_Barcode": "Metadata_Plate",'broad_sample':'Metadata_broad_sample','well_position':'Metadata_Well'})



meta_lincs2=meta_lincs.groupby(['Metadata_broad_sample','Metadata_mmoles_per_liter','Metadata_Plate','Metadata_Well']).size().reset_index()


rootDirDrug='/home/ubuntu/bucket/projects/2015_10_05_DrugRepurposing_AravindSubramanian_GolubLab_Broad/workspace'
batchName='2016_04_01_a549_48hr_batch1'
batchName_mito='2016_04_01_a549_48hr_batch1_Mito_Project'

# cp_features=mergProf_treatLevel_lincs.columns[mergProf_treatLevel_lincs.columns.str.contains("Cells_|Cytoplasm_|Nuclei_")].tolist()
              
# pert_plate_well=mergProf_treatLevel_lincs.groupby(['PERT','Metadata_Plate','Metadata_Well']).size().reset_index()
# PERTS=pert_plate_well.PERT.unique().tolist()


# def f(row):
for row in range(53,drug_list_rank_Top.shape[0]):#[1140:]:
    X_Metadata_broad_sample,X_Metadata_mmoles_per_liter=\
    drug_list_rank_Top.loc[row,['Metadata_broad_sample','Metadata_mmoles_per_liter']].values;
    
    
    pert=X_Metadata_broad_sample+"_"+str(X_Metadata_mmoles_per_liter)
    print(pert)
    pkl_fileName="/home/ubuntu/bucket/projects/2016_08_01_RadialMitochondriaDistribution_donna/workspace/drugSCprofiles_lincs_merged/"+pert;
#     if not os.path.exists(pkl_fileName):
#     if (os.stat(pkl_fileName).st_size/10e5)<350:
    if 1:
#     try:
        pert_df=meta_lincs2[(meta_lincs2["Metadata_broad_sample"]==X_Metadata_broad_sample) &\
               (meta_lincs2["Metadata_mmoles_per_liter"]==X_Metadata_mmoles_per_liter)].reset_index(drop=True)
    
        print(row,pert,pert_df.shape)    

        if pert_df.shape[0]>10:
            pert_df=pert_df.sample(5).reset_index(drop=True)    
    
    
        start_time = time.time()

        pert_df_allP0=[]
        pert_df_dmsos0=[]
        for j in range(pert_df.shape[0]):
            # p,wells="SQ00015195",["A13"]
            p,wells=pert_df.loc[j,"Metadata_Plate"],[pert_df.loc[j,"Metadata_Well"]]
            fileName=rootDirDrug+"/backend/"+batchName+"/"+p+"/"+p+".sqlite"

            fileName_mito=rootDirDrug+"/backend/"+batchName_mito+"/"+p+"/"+p+".sqlite"

            print(p,wells,fileName_mito)
            # df_p_s=readSingleCellData_sqlalch_random_image_subset(fileName,50);

            if os.path.exists(fileName) and os.path.exists(fileName_mito):
                meta_well_col_str="Image_Metadata_Well"
                df_p_s_orig=readSingleCellData_sqlalch_well_subset(fileName,wells,meta_well_col_str);

                meta_well_col_str_mito="Metadata_Well"
                df_p_s_mito=readSingleCellData_sqlalch_well_subset(fileName_mito,wells,meta_well_col_str_mito);


                if j==0:
                    com_feats=list(set(df_p_s_orig.columns.tolist()) & set(df_p_s_mito.columns.tolist()))
                    com_feats_map={}
    #                 mito_feats=[]
                    for c in com_feats:
    #                     mito_feats.append(c+'_2')
                        com_feats_map[c]=c+'_2'

                df_p_s_mito=df_p_s_mito.rename(columns=com_feats_map)
                mito_feats=df_p_s_mito.columns.tolist()
                temp_mito_df = pd.DataFrame(columns = mito_feats, index = df_p_s_orig.index)
                df_p_s_orig=pd.concat([df_p_s_orig,temp_mito_df],axis=1)

                for w in wells:
                    df_p_s_mito_w=df_p_s_mito[df_p_s_mito[meta_well_col_str_mito]==w]
                    df_p_s_orig_w=df_p_s_orig[df_p_s_orig[meta_well_col_str]==w]

                    sites=df_p_s_mito_w["Metadata_Site"].unique().tolist()
                    for s in sites:
                        df_p_s_mito_s=df_p_s_mito_w[df_p_s_mito_w["Metadata_Site"]==s]
                        df_p_s_orig_s=df_p_s_orig_w[df_p_s_orig_w["Image_Metadata_Site"]==s]  

                        for i in df_p_s_orig_s.index:
                        #     closest_row=foci_table_1['XY'].apply(lambda x: np.linalg.norm((x-dfInfoo2.loc[i,"XY"]), ord=1)).sort_values(ascending=True)[0:1]
                            xDiff=abs(df_p_s_mito_s['Nuclei_Location_Center_X_2'].values-df_p_s_orig_s.loc[i,"Nuclei_Location_Center_X"])
                            yDiff=abs(df_p_s_mito_s['Nuclei_Location_Center_Y_2'].values-df_p_s_orig_s.loc[i,"Nuclei_Location_Center_Y"])
                            min_dist=np.min(xDiff+yDiff)
                        #     if closest_row.iloc[0]<10:
                            if min_dist<10:
                        #             print(i)
                                ind3=np.argmin(xDiff+yDiff)
                                ind4=df_p_s_mito_s.iloc[ind3:ind3+1].index[0]
                        #         ind4=closest_row.index[0]
                                df_p_s_orig.loc[i,mito_feats]=df_p_s_mito_s.loc[ind4,mito_feats]


        #         xxxxxx


                pert_df_allP0.append(df_p_s_orig)

                ############# fo DMSO
                wells=meta_lincs2[(meta_lincs2["Metadata_Plate"]==p) & (meta_lincs2["Metadata_broad_sample"]=="DMSO")].Metadata_Well.tolist()

                print('here')
                meta_well_col_str="Image_Metadata_Well"
                df_cont_orig=readSingleCellData_sqlalch_well_subset(fileName,wells,meta_well_col_str);
                df_cont_orig=df_cont_orig.sample(100).reset_index(drop=True)

                meta_well_col_str_mito="Metadata_Well"
                df_cont_mito=readSingleCellData_sqlalch_well_subset(fileName_mito,wells,meta_well_col_str_mito);        

                df_cont_mito=df_cont_mito.rename(columns=com_feats_map)

                temp_mito_df = pd.DataFrame(columns = mito_feats, index = df_cont_orig.index)
                df_cont_orig=pd.concat([df_cont_orig,temp_mito_df],axis=1)

                for w in wells:
                    df_cont_mito_w=df_cont_mito[df_cont_mito[meta_well_col_str_mito]==w]
                    df_cont_orig_w=df_cont_orig[df_cont_orig[meta_well_col_str]==w]

                    sites=df_cont_orig_w["Image_Metadata_Site"].unique().tolist()
                    for s in sites:
                        df_cont_mito_s=df_cont_mito_w[df_cont_mito_w["Metadata_Site"]==s]
                        df_cont_orig_s=df_cont_orig_w[df_cont_orig_w["Image_Metadata_Site"]==s]  

                        for i in df_cont_orig_s.index:
                        #     closest_row=foci_table_1['XY'].apply(lambda x: np.linalg.norm((x-dfInfoo2.loc[i,"XY"]), ord=1)).sort_values(ascending=True)[0:1]
                            xDiff=abs(df_cont_mito_s['Nuclei_Location_Center_X_2'].values-df_cont_orig_s.loc[i,"Nuclei_Location_Center_X"])
                            yDiff=abs(df_cont_mito_s['Nuclei_Location_Center_Y_2'].values-df_cont_orig_s.loc[i,"Nuclei_Location_Center_Y"])
                            min_dist=np.min(xDiff+yDiff)
                        #     if closest_row.iloc[0]<10:
                            if min_dist<10:
                        #             print(i)
                                ind3=np.argmin(xDiff+yDiff)
                                ind4=df_cont_mito_s.iloc[ind3:ind3+1].index[0]
                        #         ind4=closest_row.index[0]
                                df_cont_orig.loc[i,mito_feats]=df_cont_mito_s.loc[ind4,mito_feats]

                pert_df_dmsos0.append(df_cont_orig)
    #         xxxxx
    #         df_p_s0=readSingleCellData_sqlalch_well_subset(fileName,wells);
    #         pert_df_dmsos0.append(df_p_s0.sample(500))

        pert_df_allP = pd.concat(pert_df_allP0)
        pert_df_dmsos = pd.concat(pert_df_dmsos0)

        print(pert,'  merged_df finished!!')
        perWellData={}
        perWellData['cp_ss']=pert_df_allP
        perWellData['cp_ss_control']=pert_df_dmsos

        a_file = open(pkl_fileName+".pkl", "wb")
        pickle.dump(perWellData, a_file)
        a_file.close()
    #     perWellData.to_pickle('/home/ubuntu/bucket/projects/2018_04_20_Rosetta/workspace/synth_l1k_ssCP_meta/'+pert);
        print(pert,'  merged_df saved!!')
        print("--- %s minutes ---" % ((time.time() - start_time)/60))
          
#     except Exception:
#         pass
        
#     return
    
    
# with Pool(13) as p:
#     p.map(f, range(53,60))




