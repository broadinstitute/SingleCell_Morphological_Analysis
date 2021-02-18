# %matplotlib notebook
import pandas as pd
import numpy as np
import seaborn as sns
sns.set(color_codes=True)
from sklearn import preprocessing
# from imblearn.over_sampling import SMOTE 
from imblearn.over_sampling import SMOTE,RandomOverSampler
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score,confusion_matrix
import matplotlib.pyplot as plt
import scipy
from sklearn.decomposition import PCA
from scipy.stats import pearsonr
import os
# import boto3
import time
from scipy import stats
from bokeh.events import Tap
from bokeh.models.widgets import PreText, Select,TextInput
from sklearn.preprocessing import MinMaxScaler
import sys
sys.path.insert(0, "/home/ubuntu/2017_09_27_RareDiseases_Taipale/")
import utils.preprocessing
from universal_divergence import estimate

from bokeh.plotting import figure, output_file, show
from bokeh.models import ColumnDataSource,LinearColorMapper
from bokeh.layouts import row, column, layout, gridplot
from bokeh.models.tools import HoverTool
from bokeh.io import output_notebook
import skimage.io
from bokeh.io import curdoc
import bokeh
from skimage.color import rgb2gray
from bokeh.palettes import d3
# from bokeh.charts import BoxPlot
import holoviews as hv
from holoviews import dim
from holoviews.streams import Pipe, Buffer
from bokeh.models import HoverTool, CustomJS, ColumnDataSource, TapTool
hv.extension('bokeh')

#############################################
rootDir='/home/ubuntu/bucket/projects/2017_09_27_RareDiseases_Taipale/workspace'
ndd='6';

############################# load annotaations 
rootDir='/home/ubuntu/bucket/projects/2017_09_27_RareDiseases_Taipale/workspace'
annot_df_2 = pd.read_excel(rootDir+'/metadata/annots_set1_set2.xlsx', sheet_name=None)['Sheet1']
metaDataPlates=annot_df_2['Metadata_Plate'].unique()
# listOfPlates0=os.listdir(rootDir+'/backend/wellsSingleCells/')[1:]
# listOfPlates0=os.listdir(rootDir+'/backend/'+profType+'PerWells'+ndd+'/')
# minInAllCells,maxInAllCells=0.0003, 0.968

############################# load rank list - widget
# dfRank4bokeh1 = pd.read_excel(rootDir+'/metadata/dtFimpact4bokeh.xlsx', sheet_name=None) 
# dfRank4bokeh1 = pd.read_excel(rootDir+'/metadata/dtFimpact4bokeh.xlsx', sheet_name=None) 
dfRank = pd.read_excel(rootDir+'/metadata/impactRank-Set1and2-20200219.xlsx', sheet_name='Sheet_1') 
dfRank4bokeh=dfRank[~dfRank['CC-All'].isnull()].reset_index(drop=True)
indices = np.argsort(dfRank4bokeh['CC-All'])
dfRank4bokeh1=dfRank4bokeh.loc[indices,:].reset_index(drop=True)

dfRank4bokeh1['ticker1']=dfRank4bokeh1['Metadata_Sample'].astype(str)+'_rw:'+dfRank4bokeh1['WT-rep'].astype(int).astype(str)+'_rm:'+dfRank4bokeh1['MT-rep'].astype(int).astype(str)


DEFAULT_TICKERS=sorted(dfRank4bokeh1['ticker1'].tolist())

ticker1 = Select(title="WT-MT pairs", value=DEFAULT_TICKERS[0], options=DEFAULT_TICKERS)


########################function to remove cells on the border
def edgeCellFilter(df_1):   
    # remove cells on the border
#     imgSize=2048
    imgSize=1080
    borderLength=100
    df_1_centCells=df_1.loc[~((df_1['Nuclei_Location_Center_X']>(imgSize-borderLength)) | \
                              (df_1['Nuclei_Location_Center_X']<(borderLength))\
                            | (df_1['Nuclei_Location_Center_Y']>(imgSize-borderLength)) | \
                              (df_1['Nuclei_Location_Center_Y']<(borderLength))),:].reset_index(drop=True)

    return df_1_centCells
    
############################ load single cell data
def load_data():
    tickVal = ticker1.value
    pairWM=tickVal.split('_rw:')[0]
    
    wt=tickVal.split(' ')[0]
    wtRep=tickVal.split('_rw:')[1][0]
    mtRep=tickVal.split('_rm:')[1][0]
    
    selDF_W=dfRank[(dfRank['Metadata_Sample']==wt)&(dfRank['rep']==int(wtRep))]

    selDF_M=dfRank4bokeh1[(dfRank4bokeh1['Metadata_Sample']==pairWM)&(dfRank4bokeh1['WT-rep']==int(wtRep))&(dfRank4bokeh1['MT-rep']==int(mtRep))]
    
    plateW=selDF_W.Metadata_Plate.tolist()[0]
    wellW=selDF_W.Metadata_Well.tolist()[0]
    plateM=selDF_M.Metadata_Plate.tolist()[0]
    wellM=selDF_M.Metadata_Well.tolist()[0]
    
#     plateW,wellW='Kinase_Mutants_1_Replicate','D11'
#     plateM,wellM='Kinase_Mutants_1','E10'

    # plateW,wellW='Replicate_24','F09'
    # plateM,wellM='Replicate_24','G03'

    annotForWellW=annot_df_2.loc[(annot_df_2['Metadata_Plate']==plateW) & (annot_df_2['Metadata_Well']==wellW),:].reset_index(drop=True)
    fileNameW=rootDir+'/backend/wellsSingleCells'+ndd+'/df_'+plateW+'_'+wellW;
    fileNameUW=rootDir+'/backend/wellsSingleCells'+ndd+'/df_control_'+plateW+'_'+wellW;

    annotForWellM=annot_df_2.loc[(annot_df_2['Metadata_Plate']==plateM) & (annot_df_2['Metadata_Well']==wellM),:].reset_index(drop=True)
    fileNameM=rootDir+'/backend/wellsSingleCells'+ndd+'/df_'+plateM+'_'+wellM;
    fileNameUM=rootDir+'/backend/wellsSingleCells'+ndd+'/df_control_'+plateM+'_'+wellM;

    if (os.path.exists(fileNameW)):
        perWellDataFilteredW=edgeCellFilter(pd.read_pickle(fileNameW).reset_index(drop=True));  
        cpFeaturesW=perWellDataFilteredW.columns[(perWellDataFilteredW.columns.str.contains("Cells_|Cytoplasm_|Nuclei_"))]# & perWellDataFilteredW.columns.str.contains("_Protein")& (perWellDataFilteredW.columns.str.contains("Intensity"))]
    else:
        perWellDataFilteredW=pd.DataFrame()
        
        
    if (os.path.exists(fileNameUW)):
        perWellDataUntransW=edgeCellFilter(pd.read_pickle(fileNameUW).reset_index(drop=True));  
        cpFeaturesW=perWellDataUntransW.columns[(perWellDataUntransW.columns.str.contains("Cells_|Cytoplasm_|Nuclei_"))]# & perWellDataUntransW.columns.str.contains("_Protein")& (perWellDataUntransW.columns.str.contains("Intensity"))]
    else:
        perWellDataUntransW=pd.DataFrame()
        
        
    if (os.path.exists(fileNameM)):
        perWellDataFilteredM=edgeCellFilter(pd.read_pickle(fileNameM).reset_index(drop=True));  
        cpFeaturesM=perWellDataFilteredM.columns[(perWellDataFilteredM.columns.str.contains("Cells_|Cytoplasm_|Nuclei_"))]#  & perWellDataFilteredM.columns.str.contains("_Protein") & (perWellDataFilteredM.columns.str.contains("Intensity"))]
    else:
        perWellDataFilteredM=pd.DataFrame()
        
        
    if (os.path.exists(fileNameUM)):
        perWellDataUntransM=edgeCellFilter(pd.read_pickle(fileNameUM).reset_index(drop=True));  
        cpFeaturesM=perWellDataUntransM.columns[(perWellDataUntransM.columns.str.contains("Cells_|Cytoplasm_|Nuclei_"))]#  & perWellDataUntransM.columns.str.contains("_Protein") & (perWellDataUntransM.columns.str.contains("Intensity"))]
    else:
        perWellDataUntransM=pd.DataFrame()
   

#     print('here1')
        
#     perWellDataUntransW=pd.read_pickle(fileNameUW).reset_index(drop=True);        

#     perWellDataFilteredM=pd.read_pickle(fileNameM).reset_index(drop=True);  
#     perWellDataUntransM=pd.read_pickle(fileNameUM).reset_index(drop=True);      

#     def fun(dff):
#         wellMin,wellMax=dff['Cells_Intensity_MinIntensity_Protein'].min(),dff['Cells_Intensity_MaxIntensity_Protein'].max()
#         dff['Cells_Intensity_MinIntensity_Protein_corrected']=dff['Cells_Intensity_MinIntensity_Protein']-wellMin+minInAllCells
#         dff['Cells_Intensity_MaxIntensity_Protein_corrected']=dff['Cells_Intensity_MaxIntensity_Protein']-wellMax+maxInAllCells
#         wellMinMaxDiff=wellMax-wellMin;
#         dff['Cells_Intensity_MeanIntensity_Protein_corrected']=dff['Cells_Intensity_MeanIntensity_Protein']-(wellMinMaxDiff/2);
#         return dff
        
#     perWellDataFilteredW,perWellDataUntransW,perWellDataFilteredM,perWellDataUntransM = \
#     (df.apply(fun) for df in [perWellDataFilteredW,perWellDataUntransW,perWellDataFilteredM,perWellDataUntransM])
    
    
    def normalizeIntFeatures(dff,dffU):
        bothDf=pd.concat([dff,dffU],sort=False,ignore_index=True)
        scaler = MinMaxScaler()
        bothDf['Cells_Intensity_MaxIntensity_Protein_corrected']=scaler.fit_transform(bothDf['Cells_Intensity_MaxIntensity_Protein'].values.reshape(-1, 1))
        bothDf['Cells_Intensity_MinIntensity_Protein_corrected']=scaler.fit_transform(bothDf['Cells_Intensity_MinIntensity_Protein'].values.reshape(-1, 1))
        wellMin,wellMax=bothDf['Cells_Intensity_MinIntensity_Protein'].min(),bothDf['Cells_Intensity_MaxIntensity_Protein'].max()
        wellMinMaxDiff=wellMax-wellMin;
    #     print(wellMinMaxDiff)
        bothDf['Cells_Intensity_MeanIntensity_Protein_corrected']=bothDf['Cells_Intensity_MeanIntensity_Protein']-(wellMinMaxDiff/2);
        dff2=bothDf.iloc[0:dff.shape[0]]
        dffU2=bothDf.iloc[dff.shape[0]:]
        return dff2,dffU2
    
    

#     if perWellDataFilteredW.shape[0]>0 and perWellDataUntransW.shape[0]>0 and\
#      perWellDataFilteredM.shape[0]>0 and perWellDataUntransM.shape[0]>0:
#     perWellDataFilteredW,perWellDataUntransW=normalizeIntFeatures(perWellDataFilteredW,perWellDataUntransW)    
#     perWellDataFilteredM,perWellDataUntransM=normalizeIntFeatures(perWellDataFilteredM,perWellDataUntransM)    
    
#     print(perWellDataFilteredW.shape,perWellDataUntransW.shape,perWellDataFilteredM.shape,perWellDataUntransM.shape)
#     cpFeaturesW=perWellDataUntransW.columns[(perWellDataUntransW.columns.str.contains("_Protein")) & perWellDataUntransW.columns.str.contains("Cells_|Cytoplasm_|Nuclei_")]# & (perWellDataUntransW.columns.str.contains("Intensity"))]
#     cpFeaturesM=perWellDataUntransM.columns[(perWellDataUntransM.columns.str.contains("_Protein")) & perWellDataUntransM.columns.str.contains("Cells_|Cytoplasm_|Nuclei_")]# & (perWellDataUntransM.columns.str.contains("Intensity"))]
    
#     cpFeaturesW=perWellDataUntransW.columns[perWellDataUntransW.columns.str.contains("Cells_|Cytoplasm_|Nuclei_")]
#     cpFeaturesM=perWellDataUntransM.columns[perWellDataUntransM.columns.str.contains("Cells_|Cytoplasm_|Nuclei_")]    
    
#     cpFeatures=list(set(cpFeaturesW).intersection(cpFeaturesM).intersection(columns4scalerA))
    cpFeatures=list(set(cpFeaturesW).intersection(cpFeaturesM))
#     locFeature2beremoved=df_meanP.columns[df_meanP.columns.str.contains("_Location_Center_X|_Location_Center_Y")]
    locFeature2beremoved=list(filter(lambda x: "_Location_Center_X" in x or "_Location_Center_Y" in x , cpFeatures)) 

    corFeature2beremoved=list(filter(lambda x: "Correlation" in x , cpFeatures)) 

    allFeatures=list(set(cpFeatures)-set(locFeature2beremoved)-set(corFeature2beremoved))
#     print('here2')
    nomalizedCellsW=perWellDataFilteredW.copy()
    nomalizedCellsM=perWellDataFilteredM.copy()
    nomalizedCellsW['Label']='nomalizedCellsW'
    nomalizedCellsM['Label']='nomalizedCellsM'

    scaledCellsW=perWellDataFilteredW.copy()
    scaledCellsM=perWellDataFilteredM.copy()
    scaledCellsW['Label']='CellsW'
    scaledCellsM['Label']='CellsM'

    scaledUntransCellsW=perWellDataUntransW.copy()
    scaledUntransCellsM=perWellDataUntransM.copy()
    scaledUntransCellsW['Label']='untransCellsW'
    scaledUntransCellsM['Label']='untransCellsM'

    perWellDataFilteredW['Label']='CellsW'
    perWellDataFilteredM['Label']='CellsM'
    perWellDataUntransW['Label']='untransCellsW'
    perWellDataUntransM['Label']='untransCellsM'
    
#     print('here3')    
    scalerW = preprocessing.StandardScaler()
    scalerM = preprocessing.StandardScaler()
    scaler = preprocessing.StandardScaler()

    if perWellDataUntransW.shape[0]>1 and perWellDataUntransM.shape[0]>1:
        scalerW.fit(perWellDataUntransW.loc[:,allFeatures].astype('float64'))
        scalerM.fit(perWellDataUntransM.loc[:,allFeatures].astype('float64'))

        nomalizedCellsW.loc[:,allFeatures]=scalerW.transform(perWellDataFilteredW.loc[:,allFeatures].astype('float64'))
        nomalizedCellsM.loc[:,allFeatures]=scalerM.transform(perWellDataFilteredM.loc[:,allFeatures].astype('float64'))
        
        nomalizedCellsMeanW=nomalizedCellsW.loc[:,allFeatures].mean().to_frame().T
        avgProfileW=pd.concat([nomalizedCellsMeanW,annotForWellW],sort=False, axis=1)
        
        nomalizedCellsMeanM=nomalizedCellsM.loc[:,allFeatures].mean().to_frame().T
        avgProfileM=pd.concat([nomalizedCellsMeanM,annotForWellM],sort=False, axis=1)
        normPcc=scipy.stats.pearsonr(np.squeeze(nomalizedCellsMeanW), np.squeeze(nomalizedCellsMeanM))[0]
        avgUnPcc=scipy.stats.pearsonr(np.squeeze(perWellDataUntransW.loc[:,allFeatures].mean().to_frame()),\
                                      np.squeeze(perWellDataUntransM.loc[:,allFeatures].mean().to_frame()))[0]
#         udUnT=estimate(perWellDataUntransW.loc[:,allFeatures].values,perWellDataUntransM.loc[:,allFeatures].values)

    else:
        normPcc=np.nan;avgUnPcc=np.nan;
        
#     scaledCellsW.loc[:,allFeatures]=scaler.fit_transform(perWellDataFilteredW.loc[:,allFeatures].astype('float64'))
#     scaledCellsM.loc[:,allFeatures]=scaler.fit_transform(perWellDataFilteredM.loc[:,allFeatures].astype('float64'))

#     scaledUntransCellsW.loc[:,allFeatures]=scaler.fit_transform(perWellDataUntransW.loc[:,allFeatures].astype('float64'))
#     scaledUntransCellsM.loc[:,allFeatures]=scaler.fit_transform(perWellDataUntransM.loc[:,allFeatures].astype('float64'))

#     print('here4')    
    WellW=pd.concat([perWellDataFilteredW,perWellDataUntransW], ignore_index=True,sort=False)
    WellM=pd.concat([perWellDataFilteredM,perWellDataUntransM], ignore_index=True,sort=False)  


#     print('here5')    
    notDrop=['Cells_Intensity_IntegratedIntensity_Protein','Cells_Intensity_MeanIntensity_Protein',\
     'Cells_Intensity_StdIntensity_Protein','Cells_Intensity_MinIntensity_Protein',\
     'Cells_Intensity_MaxIntensity_Protein',\
     'Cells_Intensity_UpperQuartileIntensity_Protein','Cells_Intensity_LowerQuartileIntensity_Protein']

    to_drop=utils.preprocessing.find_correlation(WellW.loc[:,allFeatures], threshold=0.9, remove_negative=False)
    to_drop=list(set(to_drop)-set(notDrop))
    allFeatures=list(set(allFeatures)-set(to_drop))
#     print(to_drop)
    
    WellW=WellW.drop(to_drop, axis=1)
    WellM=WellM.drop(to_drop, axis=1)

    n_pc=10
    pcaT = PCA(n_components = n_pc)
    cols2add=['PC'+str(i) for i in range(n_pc)]
#     print(cols2add,WellW_PC.shape,WellW_PC.columns)
    
    WellW_PC=pd.concat([WellW.copy(),pd.DataFrame(columns=cols2add)],sort=True)
    WellM_PC=pd.concat([WellM.copy(),pd.DataFrame(columns=cols2add)],sort=True)

    
    
    WellW_PC[cols2add]=pcaT.fit_transform(WellW.loc[:,allFeatures].astype('float64'))
    WellM_PC[cols2add]=pcaT.fit_transform(WellM.loc[:,allFeatures].astype('float64'))        
        
    
    scaledWellW=WellW.copy()
    scaledWellM=WellM.copy()
    
    scaledWellW.loc[:,allFeatures]=scaler.fit_transform(WellW.loc[:,allFeatures].astype('float64'))
    scaledWellM.loc[:,allFeatures]=scaler.fit_transform(WellM.loc[:,allFeatures].astype('float64'))
    print('here6')    
#     scaledWellW_mean=scaledWellW.loc[:,allFeatures].mean().to_frame().T
#     scaledWellM_mean=scaledWellM.loc[:,allFeatures].mean().to_frame().T

    if perWellDataFilteredW.shape[0]>1 and perWellDataFilteredM.shape[0]>1:
        scaledCellsMeanW=scaledWellW.loc[scaledWellW['Label']=='CellsW',allFeatures].mean().to_frame().T
        scaledCellsMeanM=scaledWellM.loc[scaledWellM['Label']=='CellsM',allFeatures].mean().to_frame().T
        MeanW=perWellDataFilteredW.loc[:,allFeatures].astype('float64').mean().to_frame().T
        MeanM=perWellDataFilteredM.loc[:,allFeatures].astype('float64').mean().to_frame().T
#         udT=estimate(WellW_PC.loc[WellW_PC['Label']=='scaledCellsW',cols2add].values,WellM_PC.loc[WellM_PC['Label']=='scaledCellsM',cols2add].values)
        udT=np.nan;
        avgPcc=scipy.stats.pearsonr(np.squeeze(MeanW), np.squeeze(MeanM))[0]
        scalPcc=scipy.stats.pearsonr(np.squeeze(scaledCellsMeanW), np.squeeze(scaledCellsMeanM))[0]
        
    else:
        avgPcc=np.nan;scalPcc=np.nan;udT=np.nan;

    if perWellDataUntransW.shape[0]>1 and perWellDataUntransM.shape[0]>1:
        scaledUnCellsMeanW=scaledWellW.loc[scaledWellW['Label']=='untransCellsW',allFeatures].mean().to_frame().T
        scaledUnCellsMeanM=scaledWellM.loc[scaledWellM['Label']=='untransCellsM',allFeatures].mean().to_frame().T
        scalUnPcc=scipy.stats.pearsonr(np.squeeze(scaledUnCellsMeanW), np.squeeze(scaledUnCellsMeanM))[0]
#         udUnT=estimate(WellW_PC.loc[WellW_PC['Label']=='scaledUntransCellsW',cols2add].values,\
#                WellM_PC.loc[WellM_PC['Label']=='scaledUntransCellsM',cols2add].values)
        udUnT=np.nan;
    else:
        scalUnPcc=np.nan;udUnT=np.nan;
#     print(scaledWellW.loc[scaledWellW['Label']=='scaledCellsW',allFeatures].shape,\
#           scaledWellW.loc[scaledWellW['Label']=='scaledUntransCellsW',allFeatures].shape)
    
#     nomalizedCellsMeanW=nomalizedCellsW.loc[:,allFeatures].mean().to_frame().T
#     avgProfileW=pd.concat([nomalizedCellsMeanW,annotForWellW],sort=False, axis=1)


#     print(scaledWellM.loc[scaledWellM['Label']=='scaledCellsM',allFeatures].shape,\
#           scaledWellM.loc[scaledWellM['Label']=='scaledUntransCellsM',allFeatures].shape)
    
#     nomalizedCellsMeanM=nomalizedCellsM.loc[:,allFeatures].mean().to_frame().T
#     avgProfileM=pd.concat([nomalizedCellsMeanM,annotForWellM],sort=False, axis=1)


#     print('here7')



    print('Avg profiles CC: ',avgPcc)
    print('Avg untransfected profiles CC: ',avgUnPcc)
    print('Normalized profiles CC: ',normPcc)
    print('Scaled profiles CC: ',scalPcc)
    print('Scaled untransfected profiles CC: ',scalUnPcc)
    
    data = {'CC':[avgPcc, avgUnPcc,normPcc,scalPcc,scalUnPcc,udT,udUnT]} 
    df_info = pd.DataFrame(data, index =['Avg profiles', 'Avg untransfected profiles',\
                                         'Normalized profiles', 'Scaled profiles','Scaled untr profiles','UD-Trans','UD-Untrans']) 
    stats2.text = str(df_info)

    import umap
    from sklearn.manifold import TSNE
    umapT=umap.UMAP()
    tsneT = TSNE(perplexity=10)
#     data4umap=pd.concat([perWellDataFilteredW, perWellDataFilteredM,perWellDataUntransW, perWellDataUntransM], ignore_index=True,sort=False)
    data4umap=pd.concat([WellW,WellM], ignore_index=True,sort=False)
#     data4umap=pd.concat([scaledWellW,scaledWellM], ignore_index=True,sort=False)
#     data4umap=pd.concat([scaledCellsW,scaledCellsM,scaledUntransCellsW, scaledUntransCellsM], ignore_index=True,sort=False)
#     data4umap['diffMinMaxNucCyto']=abs(data4umap['Cells_Intensity_MaxIntensity_Protein_corrected'])
    data4umap['meanInt']=data4umap['Cells_Intensity_MeanIntensity_Protein']
#     data4umap['Cells_Intensity_MaxIntensity_Protein'].max()
    data4umap['diffMinMaxNucCyto']=data4umap['Cells_Intensity_UpperQuartileIntensity_Protein']-data4umap['Cells_Intensity_MinIntensity_Protein']
    data4umap['UpperQ']=data4umap['Cells_Intensity_UpperQuartileIntensity_Protein']
    data4umap['MaxInt']=data4umap['Cells_Intensity_MaxIntensity_Protein']
    data4umap['StdInt']=data4umap['Cells_Intensity_StdIntensity_Protein']
#     data4umap['meanInt']=data4umap['Cells_Intensity_MeanIntensity_Protein']
#     data4umap['diffMinMaxNucCyto']=abs(data4umap['Cells_Intensity_MaxIntensity_Protein_corrected']
#     data4umap['meanInt']=data4umap['Cells_Intensity_MeanIntensity_Protein_corrected']
    #     data4umap=pd.concat([scaledCellsW,scaledCellsM,scaledUntransCellsW, scaledUntransCellsM], ignore_index=True,sort=False)
    pcaEnabled=0;
    if pcaEnabled:
        n_pc=60
        pcaT = PCA(n_components = n_pc)
        preprocData=pcaT.fit_transform(data4umap.loc[:,allFeatures])
        preprocData4clust=np.copy(preprocData)
    else:
        preprocData4clust=data4umap.loc[:,notDrop]
        preprocData=data4umap.loc[:,allFeatures]
    
    from sklearn.cluster import AgglomerativeClustering,SpectralClustering,KMeans
    nClus=4
#     clustering = AgglomerativeClustering(n_clusters=10).fit(preprocData)
    clustering = SpectralClustering(n_clusters=nClus,affinity='nearest_neighbors',assign_labels="discretize").fit(preprocData4clust)
#     clustering = KMeans(n_clusters=nClus).fit(preprocData)
    clusterLabels=clustering.labels_#.reshape(1,preprocData.shape[0])
#     print(clusterLabels.shape,preprocData.shape)
    
    Y = umapT.fit_transform(preprocData)
#     Y = tsneT.fit_transform(preprocData)
    tsneResDF=pd.DataFrame(index=range(Y.shape[0]),columns=['one','two','Label','clsLabel','Metadata_Plate','Metadata_Well',\
                                                            'Metadata_FieldID','ObjectNumber','diffMinMaxNucCyto','meanInt','UpperQ','MaxInt','StdInt']);
    tsneResDF.loc[:,['one','two']]=Y
    tsneResDF.loc[:,'clsLabel']=clusterLabels
    tsneResDF['clsLabel'] = tsneResDF['clsLabel'].astype(str)
    tsneResDF.loc[:,['Label','Metadata_Plate','Metadata_Well','Metadata_FieldID','ObjectNumber','diffMinMaxNucCyto','meanInt','UpperQ','MaxInt','StdInt']]=\
    data4umap[['Label','Metadata_Plate','Metadata_Well','Metadata_FieldID','ObjectNumber','diffMinMaxNucCyto','meanInt','UpperQ','MaxInt','StdInt']]
#     g=sns.scatterplot(x="one", y="two", hue="Label", data=tsneResDF)
#     print(tsneResDF.shape, data4umap.shape)
#     print(tsneResDF.dtypes)
    tsneResDF['meanInt'] = tsneResDF['meanInt'].astype(float)
    tsneResDF['StdInt'] = tsneResDF['StdInt'].astype(float)
#     print(tsneResDF[['clsLabel','meanInt']].groupby(['clsLabel']).describe())
    x=tsneResDF[['clsLabel','meanInt']].groupby(['clsLabel']).describe().reset_index()#.sort_values(by=['mean'])
    x.columns = x.columns.droplevel(0)
    sortedClusters=x.sort_values(by=['mean']).index;
    map4labelOrdering=dict()
    for c in range(nClus):
        map4labelOrdering[str(sortedClusters[c])]=str(c)
        
#     print(x.sort_values(by=['mean']))
#     print(tsneResDF['clsLabel'].head())
#     print(map4labelOrdering)
    tsneResDF['clsLabel']=tsneResDF['clsLabel'].map(map4labelOrdering)
#     print(tsneResDF['clsLabel'].head())
#     x=tsneResDF[['clsLabel','StdInt']].groupby(['clsLabel']).describe().reset_index()#.sort_values(by=['mean'])
#     x.columns = x.columns.droplevel(0)
#     print(x.sort_values(by=['mean']))
    
    return tsneResDF, data4umap

############################## usefull functions 
def get_images(imInf):
    """ 
    function to output conctenated single cell images in all channles using selected
    data by hover tools
    """
    plate0,well0,field,objectN=imInf[0],imInf[1],imInf[2],imInf[3]
#     print(source2.shape)
    data4umap3=pd.DataFrame(source2.data, columns=source2.column_names)
#     data4umap3=source2.data;
#     print(data4umap3.shape)
    dfWithWTlabels=data4umap3[(data4umap3['Metadata_Plate']==plate0) & (data4umap3['Metadata_Well']==well0) &\
                 (data4umap3['Metadata_FieldID']==field) & (data4umap3['ObjectNumber']==objectN)]
    ch_p=dfWithWTlabels['FileName_OrigProtein'].values[0];
    ch_M=dfWithWTlabels['FileName_OrigMito'].values[0];
    ch_E=dfWithWTlabels['FileName_OrigER'].values[0];
    ch_D=dfWithWTlabels['FileName_OrigDNA'].values[0];
    projectPath='/home/ubuntu/bucket/projects/2017_09_27_RareDiseases_Taipale/'

    
    
    plateName=dfWithWTlabels['Metadata_Plate'].values[0]
#     wellName=dfWithWTlabels.loc[index,'Metadata_Well']
# #         print(index,wellName)
#     fieldName=dfWithWTlabels.loc[index,'Field']
#     objectNum=dfWithWTlabels.loc[index,'ObjectNumber']
#     protLoc=dfWithWTlabels.loc[index,'manual_Annot']
#         protLoc=dfWithWTlabels.loc[index,'Metadata_Location']
    batch=dfWithWTlabels['batch'].values[0]
#     ch_seg=projectPath+'/workspace/analysis/'+batch+'/'+plateName+'/analysis/'+plateName+'-'+well0+'-'+str(field)+\
#     '/binarymask/'+well0+'_s'+str(field)+'_binarymask.png'

    ch_seg=projectPath+'/workspace/analysis/'+batch+'/'+plateName+'/analysis/'+plateName+'-'+well0+'-'+str(int(field))+\
    '/Cell_outlines/'+well0+'_s'+str(int(field))+'_cell_outlines.png'
    
    
    
#         dataFileName=projectPath+batch+'/images/'+plateName+'/'   # for dataset1
    dataFileName=projectPath+batch+'/images/'+plateName+'/Images/'   # for dataset2
    
    
    
    boxSize=100;
    xCenter=int(dfWithWTlabels['Nuclei_Location_Center_X'].values[0])
    yCenter=int(dfWithWTlabels['Nuclei_Location_Center_Y'].values[0])
    print(dataFileName+ch_M,xCenter,yCenter)
    if (xCenter>boxSize) & (yCenter>boxSize) & (os.path.exists(dataFileName+ch_p)) & (os.path.exists(dataFileName+ch_M)) &\
    (os.path.exists(dataFileName+ch_E)) & (os.path.exists(dataFileName+ch_D)):
#         print(xCenter,yCenter,boxSize)
        imP=np.squeeze(skimage.io.imread(dataFileName+ch_p))[yCenter-boxSize:yCenter+boxSize,xCenter-boxSize:xCenter+boxSize]
        imM=np.squeeze(skimage.io.imread(dataFileName+ch_M))[yCenter-boxSize:yCenter+boxSize,xCenter-boxSize:xCenter+boxSize]
        imE=np.squeeze(skimage.io.imread(dataFileName+ch_E))[yCenter-boxSize:yCenter+boxSize,xCenter-boxSize:xCenter+boxSize]
        imD=np.squeeze(skimage.io.imread(dataFileName+ch_D))[yCenter-boxSize:yCenter+boxSize,xCenter-boxSize:xCenter+boxSize]
        imSeg=np.squeeze(skimage.io.imread(ch_seg))[yCenter-boxSize:yCenter+boxSize,xCenter-boxSize:xCenter+boxSize]
        imSeg = rgb2gray(imSeg)
        print(np.squeeze(skimage.io.imread(dataFileName+ch_p)).shape)
        print(imP.max(),imM.max(),imE.max(),imD.max(),imSeg.max())
#         mP=15000;mM=1500;mE=1500;mD=1500
        mP=imP.max();mM=imM.max();mE=imE.max();mD=imD.max()
    else:
        print('Cell on the edge!')
        imP=np.ones((2*boxSize,2*boxSize));imM=np.ones((2*boxSize,2*boxSize));
        imE=np.ones((2*boxSize,2*boxSize));imD=np.ones((2*boxSize,2*boxSize));
        imSeg=np.ones((2*boxSize,2*boxSize));
#         print(imP.shape,imM.shape,imE.shape,imD.shape,imSeg.shape)
        mP=1;mM=1;mE=1;mD=1
    
    
    
    return np.concatenate([imP/mP,imM/mM,imE/mE,imD/mD,imSeg/(imSeg.max())],axis=1)
#     return np.concatenate([imP,imM,imE,imD],axis=1)

def gen_data_forHover(df):
        """ function to generate data for hover tools based on input dataframe """
        d=dict(one=df['one'],
        two=df['two'],
        Label=df['Label'],
        clsLabel=df['clsLabel'],
        Metadata_Plate=df['Metadata_Plate'],
        Metadata_Well=df['Metadata_Well'],
        Metadata_FieldID=df['Metadata_FieldID'],
        #         actual_image_number=sample['actual_image_number'],
        ObjectNumber=df['ObjectNumber'],
        diff=df['diffMinMaxNucCyto'],
        meanInt=df['meanInt'],
        UpperQ=df['UpperQ'],
        MaxInt=df['MaxInt'],
        StdInt=df['StdInt'])
        return d

################################ load data and update the printed cc results 
stats2 = PreText(text='', width=600)

# sample = tsneResDF.sample(100)
tsneResDF,data4umap=load_data();

sample = tsneResDF.copy()


################################ set up plots - dataframe
d=gen_data_forHover(sample)
# colormap = {'untransCellsW': 'blue', 'untransCellsM': 'green', 'CellsW': 'red','CellsM':'orchid'}
# colors = [colormap[x] for x in sample['Label']]
from bokeh.palettes import viridis, RdBu

palette1 = RdBu(len(sample['Label'].unique()))
colors = [palette1[int(x)] for x in sample['Label']]

d["color"] = colors
source = ColumnDataSource(data=d)


##########
source2 = ColumnDataSource(data4umap)


##########
dCls=gen_data_forHover(sample)
palette2 = viridis(len(sample['clsLabel'].unique()))
colors2 = [palette2[int(x)] for x in sample['clsLabel']]
dCls["color"] = colors2
source3 = ColumnDataSource(data=dCls)
##########

hover = HoverTool()
hover.tooltips=[
    ('Label', '@Label'),
    ('clsLabel', '@clsLabel'),
    ('Metadata_Plate', '@Metadata_Plate'),    
    ('Metadata_Well', '@Metadata_Well'),    
    ('Metadata_FieldID', '@Metadata_FieldID'),    
    ('ObjectNumber', '@ObjectNumber'),    
    ('diff', '@diff'), 
    ('meanInt', '@meanInt'),
    ('UpperQ', '@UpperQ'),
    ('MaxInt', '@MaxInt'),
    ('StdInt', '@StdInt')
]

############### intensity box plot setup using holoviews in bokeh 
renderer = hv.renderer('bokeh')
def boxInts(data):
    p1hv = hv.BoxWhisker(data, ['clsLabel'], 'meanInt', label="mean intensity per cluster").sort()
    p1hv.opts(show_legend=False, width=500, cmap='Set1',ylim=(-1, 40))
#     p1hv.opts(show_legend=False, width=500)
    return p1hv

mem_stream = Buffer(sample)
p1hv_dmap = hv.DynamicMap(boxInts, kdims=[], streams=[mem_stream])
p1 = hv.render(p1hv_dmap)
# p1 = renderer(p1hv_dmap)


##################################   plot for umap and transfetection color grouping
p = figure(tools="tap,reset",tooltips=hover.tooltips)
p.circle(x='one', y='two',
         source=source,
         size=2, color='color', legend='Label')

p.title.text = 'UMAP - Applied on WT and Mutant transfected and untrasfected single cells'
p.xaxis.axis_label = 'one'
p.yaxis.axis_label = 'two'

tap_callback = CustomJS(args={'source': source})
# p.add_tools(TapTool(callback=tap_callback)) 

################################# plot for umap and clustering color grouping
p2 = figure(tools="tap,reset",tooltips=hover.tooltips)
p2.circle(x='one', y='two',
         source=source3,
         size=2, color='color', legend='clsLabel')
# p2.add_tools(TapTool(callback=tap_callback)) 

################################################  setup single cell images plots
N, M=160*2, 800
kw = dict()
kw['x_range'] = (0,M)
kw['y_range'] = (0,N)
kw['plot_width'] = M
kw['plot_height'] = N
kw['title']='               Protein                        Mito                         \
    ER                         DNA                          Segmentation        '
ts1 = figure(tools='pan,wheel_zoom,xbox_select,reset', **kw)
ts1.axis.visible = False;ts1.xgrid.visible = False;ts1.ygrid.visible = False
#        print(images[l].shape,Ms[l],Ns[l])
color = LinearColorMapper(bokeh.palettes.gray(256))
r1=ts1.image(image=[], x=0, y=0, dw=M, dh=N, color_mapper=color)


######################################## set up callbacks

def create_figure(imageInfoList):
    """ 
    Gets image info list after user tap and
    uses loaded imaged by get_images function and update the data in single cell figure
    
    """
    
    new_data1 = dict();        
    images = get_images(imageInfoList)
    
    if r1.data_source.data['image']:
        print(r1.data_source.data['image'][0].shape)
        lastImage=r1.data_source.data['image'][0][-200:,:];
        print(lastImage.shape,images.shape)
        new_data1['image'] = [np.concatenate([lastImage,images],axis=0)]
    else:
        new_data1['image'] = [images]
    
    r1.data_source.data = new_data1
#     ts.append(ts1);
    return ts1
#     return column(children=[ts1]+[ts1])




def ticker1_change(attrname, old, new):
    """ Loads WT-MT data again and starts ... """
    print(old,new)
    tsneResDF2,data4umap2=load_data();
    
    source2.data=ColumnDataSource(data4umap2).data     #### update source2
#     print('shapeHere',tsneResDF2.shape,data4umap2.shape)
    
    d2=gen_data_forHover(tsneResDF2)
    colors = [colormap[x] for x in tsneResDF2['Label']]
    d2["color"] = colors
    source.data=d2     #### update source

    dCls2=gen_data_forHover(tsneResDF2)
    palette2 = viridis(len(tsneResDF2['clsLabel'].unique()))
    colors2 = [palette2[int(x)] for x in tsneResDF2['clsLabel']]
    dCls2['color']=colors2
    source3.data=dCls2   #### update source3
    
    mem_stream.send(tsneResDF2) #### update boxplot
    series.children=[ticker1]+[row(p,p2,column(stats2,p1))]+[]
#     update()    
    

def updateP():
    if source.selected.indices != []:
        updateSingleCellFigure(source)


def updateP2():
    if source3.selected.indices != []:
        updateSingleCellFigure(source3)

    
def updateSingleCellFigure(s):   
    """ Get selected hover information and update single cell images """
    indexActive=s.selected.indices
    imageInfoList=[s.data['Metadata_Plate'][indexActive[0]],s.data['Metadata_Well'][indexActive[0]],s.data['Metadata_FieldID'][indexActive[0]],s.data['ObjectNumber'][indexActive[0]]]
#     print('imageInfoList',imageInfoList)
    fS=create_figure(imageInfoList)
#     series.children=[ticker1]+[row(p,p2,p1)]+[fS]
    series.children=[ticker1]+[row(p,p2,column(stats2,p1))]+[fS]
    return 
    
    
p.on_event(Tap, updateP)
p2.on_event(Tap, updateP2)


ticker1.on_change('value', ticker1_change)

indd=0
imageInfoList0=[d['Metadata_Plate'][indd],d['Metadata_Well'][indd],d['Metadata_FieldID'][indd],d['ObjectNumber'][indd]]

fS=create_figure(imageInfoList0)

p.add_tools(hover)
p2.add_tools(hover)

plotss = row(p,p2,column(stats2,p1))
series = column(children=[ticker1]+[plotss]+[fS])
# series2 = row(children=[widgets]+[fS])
# show(series)
# update()

curdoc().add_root(series)







