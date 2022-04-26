import numpy as np
import pandas as pd 
import time
import sys, os
# from utils import read_data, visualize_data
from utils.read_data import *
from utils.visualize_data import *

# from SingleCell_Morphological_Analysis.utils.read_data import *
# from SingleCell_Morphological_Analysis.utils.visualize_data import *

import pandas as pd
import seaborn as sns
from sqlalchemy import create_engine
from functools import reduce
import time
from scipy.stats import pearsonr
from sklearn.cluster import KMeans
import pickle
import dataframe_image as dfi


controlLabel='control'

# from funcs.utils import readSingleCellData_sqlalch_well_subset,readSingleCellData_sqlalch_random_image_subset
selected_feature='Cells_RadialDistribution_MeanFrac_mito_tubeness_16of16'
images_dir="/home/ubuntu/bucket/projects/2015_10_05_DrugRepurposing_AravindSubramanian_GolubLab_Broad/2016_04_01_a549_48hr_batch1_compressed/"
res_dir="/home/ubuntu/bucket/projects/2016_08_01_RadialMitochondriaDistribution_donna/workspace/results/TopLincsCompoundsAnalysis/"

drug_list_rank=pd.read_excel("/home/ubuntu/bucket/projects/2016_08_01_RadialMitochondriaDistribution_donna/\
workspace/Metadata_drugRep/drugList_20210115_uncorrForSiteAgg.xlsx")
drug_list_rank=drug_list_rank[~drug_list_rank["phenotype_abundance_pval"].isnull()]
drug_list_rank_Top=drug_list_rank.loc[0:20].append(drug_list_rank.loc[drug_list_rank.shape[0]-20:]).reset_index(drop=True)
drug_list_rank_Top["n_replicates"]=drug_list_rank_Top[0]
# drug_list_rank=drug_list_rank.sort_values(by=['phenotype_abundance_t'],ascending=1)
drug_list_rank_Top["cc-all"]="";drug_list_rank_Top["cc-uncorr"]="";drug_list_rank_Top["cc-corr"]="";
drug_list_rank_Top["L2-all"]="";drug_list_rank_Top["L2-uncorr"]="";drug_list_rank_Top["L2-corr"]="";

repLevelLincs=pd.read_csv('/home/ubuntu/bucket/projects/2018_04_20_Rosetta/workspace/preprocessed_data/LINCS-Pilot1/CellPainting/replicate_level_cp_normalized_dmso.csv.gz')
repLevelLincs.Metadata_mmoles_per_liter=repLevelLincs.Metadata_mmoles_per_liter.values.round(2)

cp_features_m=repLevelLincs.columns[repLevelLincs.columns.str.contains("Cells_|Cytoplasm_|Nuclei_")].tolist()

meta_lincs=pd.read_csv("/home/ubuntu/bucket/projects/2018_04_20_Rosetta/\
workspace/results/synth_meta/meta_lincs_repLevel.csv")
meta_lincs.Metadata_mmoles_per_liter=meta_lincs.Metadata_mmoles_per_liter.values.round(2)
meta_lincs2=meta_lincs.groupby(['Metadata_broad_sample','Metadata_mmoles_per_liter','Metadata_Plate','Metadata_Well']).size().reset_index()


rootDirDrug='/home/ubuntu/bucket/projects/2015_10_05_DrugRepurposing_AravindSubramanian_GolubLab_Broad/workspace'
batchName='2016_04_01_a549_48hr_batch1_Mito_Project'

# cp_features=mergProf_treatLevel_lincs.columns[mergProf_treatLevel_lincs.columns.str.contains("Cells_|Cytoplasm_|Nuclei_")].tolist()

blackListFeatures = pickle.load( open("utils/blackListFeatures.pkl", "rb" ))
              
# pert_plate_well=mergProf_treatLevel_lincs.groupby(['PERT','Metadata_Plate','Metadata_Well']).size().reset_index()
# PERTS=pert_plate_well.PERT.unique().tolist()



# for row in range(1):#(drug_list_rank_Top.shape[0]):#[1140:]:
# row=0
for row in range(drug_list_rank_Top.shape[0]):
    X_Metadata_broad_sample,X_Metadata_mmoles_per_liter=\
    drug_list_rank_Top.loc[row,['Metadata_broad_sample','Metadata_mmoles_per_liter']].values;

    pert=X_Metadata_broad_sample+"_"+str(X_Metadata_mmoles_per_liter)
    print(row, pert)
    
    pert_df=meta_lincs[(meta_lincs["Metadata_broad_sample"]==X_Metadata_broad_sample) &\
           (meta_lincs["Metadata_mmoles_per_liter"]==X_Metadata_mmoles_per_liter)].reset_index(drop=True)
    plates=pert_df["Metadata_Plate"].tolist()
    
    d=X_Metadata_broad_sample;

    fileName="/home/ubuntu/bucket/projects/2016_08_01_RadialMitochondriaDistribution_donna/workspace/drugSCprofiles/"+pert+".pkl"    
    if os.path.exists(fileName):
        pair_dict = pickle.load( open(fileName, "rb" ) )
        df_ss0=pair_dict['cp_ss']
        pert_df_dmsos=pair_dict['cp_ss_control']

        cp_features=pert_df_dmsos.columns[pert_df_dmsos.columns.str.contains("Cells_|Cytoplasm_|Nuclei_")].tolist()


        cols2remove0=[i for i in df_ss0 if ((df_ss0[i]=='nan').sum(axis=0)/df_ss0.shape[0])>0.05]
        cols2remove1=df_ss0[cp_features].std()[df_ss0[cp_features].std() < 0.0001].index.tolist()
#         
#         corFeature2beremoved=list(filter(lambda x: "Correlation" in x , cp_features))         
        cols2remove_drug=cols2remove0+cols2remove1+blackListFeatures
#         print(cols2remove_drug)

        cols2remove0=[i for i in pert_df_dmsos if ((df_ss0[i]=='nan').sum(axis=0)/pert_df_dmsos.shape[0])>0.05]
        cols2remove1=pert_df_dmsos[cp_features].std()[pert_df_dmsos[cp_features].std() < 0.0001].index.tolist()
        cols2remove_dmso=cols2remove0+cols2remove1+blackListFeatures
#         print(cols2remove_dmso)

        cols2remove=list(set(cols2remove_dmso+cols2remove_drug))
        df_ss0=df_ss0.drop(cols2remove, axis=1,errors='ignore');
        pert_df_dmsos=pert_df_dmsos.drop(cols2remove, axis=1,errors='ignore');

        df_ss0=df_ss0.interpolate()
        pert_df_dmsos=pert_df_dmsos.interpolate()

        cp_features=pert_df_dmsos.columns[pert_df_dmsos.columns.str.contains("Cells_|Cytoplasm_|Nuclei_")].tolist()
        
        locFeature2beremoved=list(filter(lambda x: "_X" in x or "_Y" in x , cp_features))   
        cpFeatures4scale=list(set(cp_features)-set(locFeature2beremoved))
#         cpFeatures4scale=list((set(cp_features) & set(cp_features_m))-set(locFeature2beremoved))

        scaler = preprocessing.StandardScaler()
        df_ss0_scaled = df_ss0.copy()
        pert_df_dmsos_scaled = pert_df_dmsos.copy()
#         scaler.fit(pd.concat([df_ss0[cpFeatures4scale],pert_df_dmsos[cpFeatures4scale]], ignore_index=True))
        scaler.fit(pert_df_dmsos[cpFeatures4scale].values)
        df_ss0_scaled[cpFeatures4scale]=scaler.transform(df_ss0[cpFeatures4scale].values)     
        pert_df_dmsos_scaled[cpFeatures4scale]=scaler.transform(pert_df_dmsos[cpFeatures4scale].values) 
        
        df_ss0 = df_ss0_scaled.copy()
        pert_df_dmsos = pert_df_dmsos_scaled.copy()
        
#         mean_dmso=pert_df_dmsos[cp_features].mean().values
#         mean_drug=df_ss0[cp_features].mean().values


#         if row==0:
        columns = ['CC'];
        df_f_MI_CC= pd.DataFrame(index=cpFeatures4scale,columns=columns)
        refFeature=pert_df_dmsos[selected_feature].values
        otherFeatures=cpFeatures4scale
        for f in otherFeatures:
        #     print(f)
            comparedFeature=pert_df_dmsos[f].values;
            df_f_MI_CC.loc[f,'CC']=pearsonr(refFeature, comparedFeature)[0]    

#         uncorr_features = df_f_MI_CC[df_f_MI_CC.abs()['CC']<0.2].index.tolist()
#         corr_features = df_f_MI_CC[df_f_MI_CC.abs()['CC']>0.5].index.tolist()    

        uncorr_features = df_f_MI_CC.abs().sort_values(by=["CC"])[0:600].index.tolist()
        corr_features = df_f_MI_CC.abs().sort_values(by=["CC"])[-10:].index.tolist()


        #     from sklearn.metrics.pairwise import euclidean_distances
        #     dist2=euclidean_distances(mean_dmso.reshape(1, -1),mean_drug.reshape(1, -1))

        print("uncorr_features",len(uncorr_features))
        print("corr_features",corr_features)
        
        
        mean_dmso=pert_df_dmsos[cpFeatures4scale].mean().values
        mean_drug=df_ss0[cpFeatures4scale].mean().values
    
        
#         mean_dmso=repLevelLincs[(repLevelLincs["Metadata_Plate"].isin(plates)) &\
#                      (repLevelLincs["Metadata_broad_sample"]=="DMSO")][cpFeatures4scale].mean().values       
#         mean_drug=repLevelLincs[(repLevelLincs["Metadata_broad_sample"]==X_Metadata_broad_sample) &\
#            (repLevelLincs["Metadata_mmoles_per_liter"]==X_Metadata_mmoles_per_liter)][cpFeatures4scale].mean().values    
  
        
        drug_list_rank_Top.loc[row,["L2-all","cc-all"]]=\
        [np.linalg.norm(mean_dmso-mean_drug),pearsonr(mean_dmso, mean_drug)[0]]

        mean_dmso=pert_df_dmsos[uncorr_features].mean().values
        mean_drug=df_ss0[uncorr_features].mean().values
        drug_list_rank_Top.loc[row,["L2-uncorr","cc-uncorr"]]=\
        [np.linalg.norm(mean_dmso-mean_drug),pearsonr(mean_dmso, mean_drug)[0]]

        mean_dmso=pert_df_dmsos[corr_features].mean().values
        mean_drug=df_ss0[corr_features].mean().values
        drug_list_rank_Top.loc[row,["L2-corr","cc-corr"]]=\
        [np.linalg.norm(mean_dmso-mean_drug),pearsonr(mean_dmso, mean_drug)[0]]

        df_ss0["label"]=X_Metadata_broad_sample
        pert_df_dmsos["label"]="control"

        wtANDmtDf=pd.concat([df_ss0,pert_df_dmsos], ignore_index=True)
        channels=["Mito","AGP","DNA","RNA","ER"]
        for ch in channels:
            wtANDmtDf["PathName_Orig"+ch]=images_dir+"images/"+wtANDmtDf["Metadata_Plate"]
            wtANDmtDf["FileName_Orig"+ch]=wtANDmtDf["FileName_Orig"+ch].apply(lambda x: x.replace("tiff","png"))

        imgSize=wtANDmtDf.Metadata_ImageSizeX.values[0]
        borderLength=int(np.percentile(wtANDmtDf.Cells_AreaShape_MajorAxisLength.values, 90)/2);    
        wtANDmtDf=edgeCellFilter2(wtANDmtDf,imgSize,borderLength);  
        perc99th=np.percentile(wtANDmtDf[selected_feature].values, 99);  
        wtANDmtDf=wtANDmtDf[wtANDmtDf[selected_feature]<perc99th]

#         os.mkdir(res_dir+'/'+controlLabel+'-'+d)
        os.system("mkdir "+res_dir+'/'+controlLabel+'-'+pert)

        fh=sns.displot(wtANDmtDf, x=selected_feature, hue="label",bins=40,facet_kws={"legend_out": False});
        fh.savefig(res_dir+'/'+controlLabel+'-'+pert+'/target_hist.png')  


        df_styled=drug_list_rank_Top.loc[row:row,["Metadata_broad_sample","Metadata_mmoles_per_liter",\
            "Metadata_pert_name","Metadata_moa","Count_Cells","phenotype_abundance_pval",\
            "phenotype_abundance_t","n_replicates","L2-uncorr"]].T

        dfi.export(df_styled, res_dir+'/'+controlLabel+'-'+pert+'/df_info.png',table_conversion='matplotlib')

        if not os.path.exists(res_dir+controlLabel+'-'+pert):#+"/clusterDensity.png"):

            #### sub pop analysis
            wtANDmtDf_scaled = wtANDmtDf.copy()
            scaler = preprocessing.RobustScaler()
            locFeature2beremoved=list(filter(lambda x: "_Location_Center_X" in x or "_Location_Center_Y" in x, cp_features)) 
            cpFeatures4scale=list(set(cp_features)-set(locFeature2beremoved))
            wtANDmt=scaler.fit_transform(wtANDmtDf.loc[:,cpFeatures4scale].astype('float64'))
            wtANDmtDf_scaled[cpFeatures4scale]=wtANDmt

            nSampleSCs=6;boxSize=160;
            perc80th=np.percentile(wtANDmtDf_scaled[selected_feature].values, 80);  
            perc20th=np.percentile(wtANDmtDf_scaled[selected_feature].values, 20);  

            samples2plot_l=wtANDmtDf_scaled[wtANDmtDf_scaled[selected_feature]<perc20th].sample(nSampleSCs).reset_index(drop=True)
            samples2plot_h=wtANDmtDf_scaled[wtANDmtDf_scaled[selected_feature]>perc80th].sample(nSampleSCs).reset_index(drop=True)


            f=visualize_n_SingleCell(channels,samples2plot_l,boxSize,title="Low Target Feature Value")
            f.savefig(res_dir+'/'+controlLabel+'-'+pert+'/low_target_examplar.png')     

            f=visualize_n_SingleCell(channels,samples2plot_h,boxSize,title="high Target Feature Value")
            f.savefig(res_dir+'/'+controlLabel+'-'+pert+'/high_target_examplar.png')     


            nClus=8
            kmeans = KMeans(n_clusters=nClus, random_state=0).fit(wtANDmtDf_scaled[corr_features].values)
            clusterLabels=kmeans.labels_
            distances=kmeans.fit_transform(wtANDmtDf_scaled[corr_features].values)
            wtANDmtDf_scaled['clusterLabels']=clusterLabels;
            wtANDmtDf_scaled['dist2Mean']=np.min(distances,1);

            DirsDict={}
            DirsDict['root']=images_dir
            DirsDict['resDir']=res_dir+controlLabel+'-'+pert
            # channels=['Mito','AGP','Brightfield','ER','DNA','Outline']
            # channels=["DNA","RNA","Mito","ER","AGP"]
            disLabels=['KO']
            controlLabel='control'
            d=X_Metadata_broad_sample;
            boxSize=160

            clusteringHists(DirsDict,wtANDmtDf_scaled,controlLabel,d,nClus,corr_features,channels,boxSize)

drug_list_rank_Top.to_csv(res_dir+'/drug_list_rank_Top5.csv',index=False)