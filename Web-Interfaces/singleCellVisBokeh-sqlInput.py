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
import time
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
# rootDir='/home/ubuntu/bucket/projects/2017_10_19_Profiling_rare_ORFs/workspace'
# ndd='7';

#############################################
# rootDir='/home/ubuntu/bucket/projects/2017_09_27_RareDiseases_Taipale/workspace'
# ndd='6';

# plateName='RC4_IF_15'
# SQL_file_path=rootDir+'/backend/PILOT_1_maxproj/'+plateName+'/'+plateName+'.sqlite'
# batchName=SQL_file_path.split('backend/')[1].split('/')[0]
############################# load sql file 
# SQL_file_path='/home/ubuntu/bucket/projects/2017_10_19_Profiling_rare_ORFs/workspace/backend/20200113_96W_CP127/20X_CP_CP127_1/20X_CP_CP127_1.sqlite'
# SQL_file_path="/home/ubuntu/bucket/projects/2017_10_19_Profiling_rare_ORFs/workspace/backend/20200303_96W_CP157A/20X_CP_CP157A_2/20X_CP_CP157A_2.sqlite"
# batchName=SQL_file_path.split('backend/')[1].split('/')[0]

SQL_file_path='/home/ubuntu/bucket/projects/2017_10_19_Profiling_rare_ORFs/workspace/backend/20200615_96W_CP173/20X_CP_CP173_2/20X_CP_CP173_2.sqlite'
batchName=SQL_file_path.split('backend/')[1].split('/')[0]
########################function to remove cells on the border
# def edgeCellFilter(df_1):   
#     # remove cells on the border
#     imgSize=2048
# #     imgSize=1080
#     borderLength=100
#     df_1_centCells=df_1.loc[~((df_1['Nuclei_Location_Center_X']>(imgSize-borderLength)) | \
#                               (df_1['Nuclei_Location_Center_X']<(borderLength))\
#                             | (df_1['Nuclei_Location_Center_Y']>(imgSize-borderLength)) | \
#                               (df_1['Nuclei_Location_Center_Y']<(borderLength))),:].reset_index(drop=True)

#     return df_1_centCells

def edgeCellFilter(df_1):   
    # remove cells on the border
#     imgSize=1080
#     imgSize=2048
    print(df_1.columns[df_1.columns.str.contains("Metadata_ImageSizeX")])
    print(df_1.columns[df_1.columns.str.contains("ImageSize")])
    imgSize=df_1.Metadata_ImageSizeX.values[0]
    borderLength=int(np.percentile(df_1.Cells_AreaShape_MajorAxisLength.values, 90)/2);
    print(imgSize,borderLength)
    df_1_centCells=df_1.loc[~((df_1['Nuclei_Location_Center_X']>(imgSize-borderLength)) | \
                              (df_1['Nuclei_Location_Center_X']<(borderLength))\
                            | (df_1['Nuclei_Location_Center_Y']>(imgSize-borderLength)) | \
                              (df_1['Nuclei_Location_Center_Y']<(borderLength))),:].reset_index(drop=True)

    return df_1_centCells,borderLength

#####################################
def readSingleCellData(fileName):
    import rpy2.robjects as robjects
    from rpy2.robjects import pandas2ri
    ts=robjects.r('ts')
    from rpy2.robjects.packages import importr
    from rpy2.robjects.conversion import localconverter

    rstring="""
    function(fileName){
      
    library(tidyverse)
    library(magrittr)
    backend <- file.path(fileName)
    
    db <- src_sqlite(path = backend)
    
    image <- tbl(src = db, "Image")
    
    object <-
      tbl(src = db, "Cells") %>%
      inner_join(tbl(src = db, "Cytoplasm"),
                 by = c("TableNumber", "ImageNumber", "ObjectNumber")) %>%
      inner_join(tbl(src = db, "Nuclei"),
                 by = c("TableNumber", "ImageNumber", "ObjectNumber"))
    
    object %<>% inner_join(image, by = c("TableNumber", "ImageNumber")) 
    variables <-
      colnames(object) %>%
      stringr::str_subset("^Nuclei_|^Cells_|^Cytoplasm_")
    dt <- object %<>%
      collect()
    return(dt)
    }
    """
    
#     fileName=rootPath+"/backend/"+batchName+"/"+plateName+"/"+plateName+".sqlite"
    print(fileName)
    
    if not os.path.isfile(fileName):
        print('squlite file not exist!')
        return pd.DataFrame()
    
    rfunc=robjects.r(rstring)
    
    rdata=ts(fileName)
    r_df=rfunc(rdata)

    with localconverter(robjects.default_converter + pandas2ri.converter):
      pd_df = robjects.conversion.rpy2py(r_df)

#     pd_df = pandas2ri.ri2py_dataframe(r_df)
    #pd_df=pd_df.replace('nan', np.nan)
    cols2remove=[i for i in pd_df.columns.tolist() if ((pd_df[i]=='nan').sum(axis=0)/pd_df.shape[0])>0.05]
    print(cols2remove)
    pd_df=pd_df.drop(cols2remove, axis=1);
    #     pd_df = pd_df.fillna(pd_df.median())
#     print(1)
    pd_df=pd_df.replace('nan', np.nan)
#     print(2)
#     pd_df = pd_df.interpolate()
#     print(3)
    pd_df,borderLength=edgeCellFilter(pd_df);  
#     print(4)
    return pd_df,borderLength

start_time = time.time()
dfWithWTlabels,borderLength=readSingleCellData(SQL_file_path)
dfWithWTlabels['batch']=batchName
print("--- %s seconds ---" % (time.time() - start_time))
# pd_df=pd.read_pickle('/home/ubuntu/bucket/projects/2017_10_19_Profiling_rare_ORFs/workspace/backend/wellsSingleCells7/df_20X_CP_CP127_1_F09', compression='infer').reset_index(drop=True); 
# edgeCellFilter

cpFeatures=dfWithWTlabels.columns[dfWithWTlabels.columns.str.contains("Cells_|Cytoplasm_|Nuclei_")]# & (perWellDataUntransW.columns.str.contains("Intensity"))]

# cpFeatures=list(set(cpFeatures)-set(blackListFeatures))
locFeature2beremoved=list(filter(lambda x: "_Location_Center_X" in x or "_Location_Center_Y" in x , cpFeatures)) 
cpFeatures4scale=list(set(cpFeatures)-set(locFeature2beremoved))


scaler = preprocessing.RobustScaler()
dataScaled=scaler.fit_transform(dfWithWTlabels.loc[:,cpFeatures4scale])
dfWithWTlabels_scaled = dfWithWTlabels.copy()
dfWithWTlabels_scaled[cpFeatures4scale]=dataScaled


############################# form ticker 1
dfRank4bokeh1=dfWithWTlabels_scaled.groupby(['Metadata_Plate','Metadata_Well']).size().reset_index()
dfRank4bokeh1['ticker1']=dfRank4bokeh1['Metadata_Plate'].astype(str)+'__'+dfRank4bokeh1['Metadata_Well'].astype(str)#+'__'+dfRank4bokeh1['Metadata_Location']


DEFAULT_TICKERS=sorted(dfRank4bokeh1['ticker1'].tolist())

ticker1 = Select(title="WT-MT pairs", value=DEFAULT_TICKERS[0], options=DEFAULT_TICKERS)

    

############################ load single cell data
def load_data():  # get well data
    tickVal = ticker1.value
    well=tickVal.split('__')[1]
    
    wellData=dfWithWTlabels_scaled[dfWithWTlabels_scaled['Metadata_Well']==well]
    
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
    cpFeatures=wellData.columns[(wellData.columns.str.contains("_Protein")) & wellData.columns.str.contains("Cells_|Cytoplasm_|Nuclei_")]# & (perWellDataUntransW.columns.str.contains("Intensity"))]
#     cpFeaturesM=perWellDataUntransM.columns[(perWellDataUntransM.columns.str.contains("_Protein")) & perWellDataUntransM.columns.str.contains("Cells_|Cytoplasm_|Nuclei_")]# & (perWellDataUntransM.columns.str.contains("Intensity"))]
    
#     cpFeaturesW=perWellDataUntransW.columns[perWellDataUntransW.columns.str.contains("Cells_|Cytoplasm_|Nuclei_")]
#     cpFeaturesM=perWellDataUntransM.columns[perWellDataUntransM.columns.str.contains("Cells_|Cytoplasm_|Nuclei_")]    
    
#     cpFeatures=list(set(cpFeaturesW).intersection(cpFeaturesM).intersection(columns4scalerA))
#     cpFeatures=list(set(cpFeaturesW).intersection(cpFeaturesM))
#     locFeature2beremoved=df_meanP.columns[df_meanP.columns.str.contains("_Location_Center_X|_Location_Center_Y")]
    locFeature2beremoved=list(filter(lambda x: "_Location_Center_X" in x or "_Location_Center_Y" in x , cpFeatures)) 

    corFeature2beremoved=list(filter(lambda x: "Correlation" in x , cpFeatures)) 

    allFeatures=list(set(cpFeatures)-set(locFeature2beremoved)-set(corFeature2beremoved))

    import umap
    from sklearn.manifold import TSNE
    umapT=umap.UMAP()
    tsneT = TSNE(perplexity=10)
#     data4umap=pd.concat([perWellDataFilteredW, perWellDataFilteredM,perWellDataUntransW, perWellDataUntransM], ignore_index=True,sort=False)
    
    data4umap=pd.concat([wellData], ignore_index=True,sort=False)
#     data4umap['Label']='untransCellsW'
    data4umap['Label']=data4umap['Metadata_Site'].astype(str)
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
#         preprocData4clust=data4umap.loc[:,notDrop]
        preprocData4clust=data4umap.loc[:,allFeatures]
        preprocData=data4umap.loc[:,allFeatures]
    
    from sklearn.cluster import AgglomerativeClustering,SpectralClustering,KMeans
    nClus=8
#     clustering = AgglomerativeClustering(n_clusters=10).fit(preprocData)
    clustering = SpectralClustering(n_clusters=nClus,affinity='nearest_neighbors',assign_labels="discretize").fit(preprocData4clust)
#     clustering = KMeans(n_clusters=nClus).fit(preprocData)
    clusterLabels=clustering.labels_#.reshape(1,preprocData.shape[0])
#     print(clusterLabels.shape,preprocData.shape)
    
    Y = umapT.fit_transform(preprocData)
#     Y = tsneT.fit_transform(preprocData)
    tsneResDF=pd.DataFrame(index=range(Y.shape[0]),columns=['one','two','Label','clsLabel','Metadata_Plate','Metadata_Well',\
                                                            'Metadata_Site','ObjectNumber','diffMinMaxNucCyto','meanInt','UpperQ','MaxInt','StdInt']);
    tsneResDF.loc[:,['one','two']]=Y
    tsneResDF.loc[:,'clsLabel']=clusterLabels
    tsneResDF['clsLabel'] = tsneResDF['clsLabel'].astype(str)
    tsneResDF.loc[:,['Label','Metadata_Plate','Metadata_Well','Metadata_Site','ObjectNumber','diffMinMaxNucCyto','meanInt','UpperQ','MaxInt','StdInt']]=\
    data4umap[['Label','Metadata_Plate','Metadata_Well','Metadata_Site','ObjectNumber','diffMinMaxNucCyto','meanInt','UpperQ','MaxInt','StdInt']]
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
#     print(data4umap3.columns[data4umap3.columns.str.contains('FileName_Orig')])
    
    dfWithWTlabels=data4umap3[(data4umap3['Metadata_Plate']==plate0) & (data4umap3['Metadata_Well']==well0) &\
                 (data4umap3['Metadata_Site']==field) & (data4umap3['ObjectNumber']==objectN)]
    ch_p=dfWithWTlabels['FileName_OrigProtein'].values[0];
    ch_M=dfWithWTlabels['FileName_OrigMito'].values[0];
    ch_E=dfWithWTlabels['FileName_OrigER'].values[0];
    ch_D=dfWithWTlabels['FileName_OrigDNA'].values[0];
#     projectPath='/home/ubuntu/bucket/projects/2017_09_27_RareDiseases_Taipale/'
#     projectPath='/home/ubuntu/bucket/projects/2017_10_19_Profiling_rare_ORFs/'
    projectPath=rootDir.split('workspace')[0]

    
    
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

#     ch_seg=projectPath+'/workspace/analysis/'+batch+'/'+plateName+'/analysis/'+plateName+'-'+well0+'-'+str(int(field))+\
#     '/Cell_outlines/'+well0+'_s'+str(int(field))+'_cell_outlines.png'
    
    ch_seg=projectPath+'/workspace/analysis/'+batch+'/'+plateName+'/analysis/'+plateName+'-'+well0+\
    '/Cell_outlines/'+well0+'_s'+str(int(field))+'_cell_outlines.png'    
    
    dataFileName=projectPath+batch+'/images/'+plateName+'/'   # for dataset1
#     dataFileName=projectPath+batch+'/images/'+plateName+'/Images/'   # for dataset2
#     dataFileName=projectPath+batch+'/'+plateName+'/'   # for varsha arrayed data    
    
    
    boxSize=borderLength;
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
        imP=np.ones((boxSize,boxSize));imM=np.ones((boxSize,boxSize));
        imE=np.ones((boxSize,2*boxSize));imD=np.ones((boxSize,boxSize));
        imSeg=np.ones((boxSize,boxSize));
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
        Metadata_Site=df['Metadata_Site'],
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
from bokeh.palettes import viridis, inferno,all_palettes

uniqLabels=sample['Label'].unique().tolist()
# all_palettes['Viridis'][4]
palette1 = inferno(len(uniqLabels))
colors = [palette1[int(uniqLabels.index(x))] for x in sample['Label']]

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
    ('Metadata_Site', '@Metadata_Site'),    
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
    
    # randomely select 3 sites
    uniqLabels0=tsneResDF2['Label'].unique().tolist()
    randIntArr=np.random.randint(1,len(uniqLabels0), size=3)
    randSites=[uniqLabels0[i] for i in randIntArr]
    tsneResDF2m=tsneResDF2[tsneResDF2['Label'].isin(randSites)].reset_index(drop=True)
    
    
    source2.data=ColumnDataSource(data4umap2).data     #### update source2
#     print('shapeHere',tsneResDF2.shape,data4umap2.shape)
    
    d2=gen_data_forHover(tsneResDF2m)
    uniqLabels=tsneResDF2m['Label'].unique().tolist()
#     palette1 = inferno(len(uniqLabels))
#     colors = [palette1[int(uniqLabels.index(x))] for x in tsneResDF2m['Label']]
    palette1 = all_palettes['Colorblind'][len(uniqLabels)]
    colors = [palette1[int(uniqLabels.index(x))] for x in tsneResDF2m['Label']]        
    
        
#     colors = [colormap[x] for x in tsneResDF2['Label']]
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
    imageInfoList=[s.data['Metadata_Plate'][indexActive[0]],s.data['Metadata_Well'][indexActive[0]],s.data['Metadata_Site'][indexActive[0]],s.data['ObjectNumber'][indexActive[0]]]
#     print('imageInfoList',imageInfoList)
    fS=create_figure(imageInfoList)
#     series.children=[ticker1]+[row(p,p2,p1)]+[fS]
    series.children=[ticker1]+[row(p,p2,column(stats2,p1))]+[fS]
    return 
    
    
p.on_event(Tap, updateP)
p2.on_event(Tap, updateP2)


ticker1.on_change('value', ticker1_change)

indd=0
imageInfoList0=[d['Metadata_Plate'][indd],d['Metadata_Well'][indd],d['Metadata_Site'][indd],d['ObjectNumber'][indd]]

fS=create_figure(imageInfoList0)

p.add_tools(hover)
p2.add_tools(hover)

plotss = row(p,p2,column(stats2,p1))
series = column(children=[ticker1]+[plotss]+[fS])
# series2 = row(children=[widgets]+[fS])
# show(series)
# update()

curdoc().add_root(series)




