# %matplotlib notebook
# sns.set(color_codes=True)
import pandas as pd
import numpy as np
# import seaborn as sns
from PyQt5 import QtCore, QtGui, QtWidgets
# from sklearn import preprocessing
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


#############################################
# rootDir='/home/ubuntu/bucket/projects/2017_09_27_RareDiseases_Taipale/workspace'
# ndd='6';

############################# load annotaations 
SQL_file_path='/home/ubuntu/bucket/projects/2017_10_19_Profiling_rare_ORFs/workspace/backend/20200113_96W_CP127/20X_CP_CP127_1/20X_CP_CP127_1.sqlite'

########################function to remove cells on the border
def edgeCellFilter(df_1):   
    # remove cells on the border
    imgSize=1080
    borderLength=100
    df_1_centCells=df_1.loc[~((df_1['Nuclei_Location_Center_X']>(imgSize-borderLength)) | \
                              (df_1['Nuclei_Location_Center_X']<(borderLength))\
                            | (df_1['Nuclei_Location_Center_Y']>(imgSize-borderLength)) | \
                              (df_1['Nuclei_Location_Center_Y']<(borderLength))),:].reset_index(drop=True)

    return df_1_centCells

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
    pd_df=edgeCellFilter(pd_df);  
#     print(4)
    return pd_df
pd_df=readSingleCellData(SQL_file_path)

############################# load rank list 
dfRank4bokeh1=pd_df.groupby(['Metadata_Plate','Metadata_Well']).size().reset_index()
dfRank4bokeh1['ticker1']=dfRank4bokeh1['Metadata_Plate'].astype(str)+'__'+dfRank4bokeh1['Metadata_Well'].astype(str)


DEFAULT_TICKERS=sorted(dfRank4bokeh1['ticker1'].tolist())

ticker1 = Select(title="WT-MT pairs", value=DEFAULT_TICKERS[0], options=DEFAULT_TICKERS)

    
############################ load single cell data
def load_data():  # get well data
    tickVal = ticker1.value
    well=tickVal.split('__')[1]
    
    wellData=pd_df[pd_df['Metadata_Well']==well]
    
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


# plate0,well0=plateM,wellM
# plate0,well0=plateW,wellW
# field=21; objectN=10;
# field=5; objectN=16;
# field=3; objectN=55;
# field=15; objectN=38;
# field=11; objectN=40;


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
# output_notebook()
# output_file('columndatasource_example.html')

# output_file("color_scatter.html", title="color_scatter.py example", mode="cdn")
# df = pd.read_csv('thor_wwii.csv')
stats2 = PreText(text='', width=600)
# sample = tsneResDF.sample(100)
tsneResDF,data4umap=load_data();

sample = tsneResDF.copy()
colormap = {'untransCellsW': 'blue', 'untransCellsM': 'green', 'CellsW': 'red','CellsM':'orchid'}
colors = [colormap[x] for x in sample['Label']]

# sourceBox= ColumnDataSource(sample)
# source = ColumnDataSource(sample)
# source = ColumnDataSource(sample)
d=dict(one=sample['one'],
        two=sample['two'],
        color=colors,
        Label=sample['Label'],
        clsLabel=sample['clsLabel'],
        Metadata_Plate=sample['Metadata_Plate'],
        Metadata_Well=sample['Metadata_Well'],
        Metadata_FieldID=sample['Metadata_FieldID'],
        #         actual_image_number=sample['actual_image_number'],
        ObjectNumber=sample['ObjectNumber'],
        diff=sample['diffMinMaxNucCyto'],
        meanInt=sample['meanInt'],
        UpperQ=sample['UpperQ'],
        MaxInt=sample['MaxInt'],
        StdInt=sample['StdInt']
    )

# dCls=d.copy()
from bokeh.palettes import viridis

palette2 = viridis(len(sample['clsLabel'].unique()))
# print(sample['clsLabel'].unique())
# color_map2 = bokeh.models.CategoricalColorMapper(factors=sample['clsLabel'].unique().tolist(),palette=palette2)

colors2 = [palette2[int(x)] for x in sample['clsLabel']]
# dCls['color']=color_map

source = ColumnDataSource(data=d)

source2 = ColumnDataSource(data4umap)

# print(len(d),len(dCls))
dCls=dict(one=sample['one'],
        two=sample['two'],
        color=colors2,
        Label=sample['Label'],
        clsLabel=sample['clsLabel'],
        Metadata_Plate=sample['Metadata_Plate'],
        Metadata_Well=sample['Metadata_Well'],
        Metadata_FieldID=sample['Metadata_FieldID'],
        #         actual_image_number=sample['actual_image_number'],
        ObjectNumber=sample['ObjectNumber'],
        diff=sample['diffMinMaxNucCyto'],
        meanInt=sample['meanInt'],
        UpperQ=sample['UpperQ'],
        MaxInt=sample['MaxInt'],
        StdInt=sample['StdInt']
    )
source3 = ColumnDataSource(data=dCls)


measured1=sample['meanInt'].values;
hist1, edges1 = np.histogram(measured1, density=True, bins=50)
# kde1 = stats.gaussian_kde(measured1.astype(float))

dDist=dict(hist_1=hist1,
           edges_1_right=edges1[1:],
           edges_1_left=edges1[:-1])
#            kde_1=kde1)

sourceDis=ColumnDataSource(data=dDist)

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

# p1 = figure(title='dis - meanInt', tools='', background_fill_color="#fafafa")
# # p1.quad(top='hist_1', bottom=0, left='edges_1_left', right='edges_1_right',
# #        fill_color="navy", line_color="white", alpha=0.5,source=sourceDis)
# p1 = BoxPlot(sample, values = "meanInt", label = "clsLabel",  
#             color = "blue", title = "mean intensity per cluster", 
#              legend = "top_right") 

import holoviews as hv
from holoviews import dim
from holoviews.streams import Pipe, Buffer
hv.extension('bokeh')
renderer = hv.renderer('bokeh')

# pipe = Pipe(data=sample)


def boxInts(data):
    p1hv = hv.BoxWhisker(data, ['clsLabel'], 'meanInt', label="mean intensity per cluster").sort()
    p1hv.opts(show_legend=False, width=500, cmap='Set1',ylim=(-1, 20))
#     p1hv.opts(show_legend=False, width=500)
    return p1hv

mem_stream = Buffer(sample)
# When run live, this cell's output should match the behavior of the GIF below
p1hv_dmap = hv.DynamicMap(boxInts, kdims=[], streams=[mem_stream])
p1 = hv.render(p1hv_dmap)
# p1 = renderer(p1hv_dmap)

# p1.line(x, pdf1, line_color="#ff8888", line_width=4, alpha=0.7, legend="PDF")

#p1.y_range.start = 0;p1.legend.location = "top_left";p1.legend.background_fill_color = "#fefefe"
# p1.xaxis.axis_label = 'x';p1.yaxis.axis_label = 'Pr(x)';p1.grid.grid_line_color="white"

##################################
p = figure(tools="tap,reset",tooltips=hover.tooltips)
p.circle(x='one', y='two',
         source=source,
         size=2, color='color', legend='Label')

p.title.text = 'UMAP - Applied on WT and Mutant transfected and untrasfected single cells'
p.xaxis.axis_label = 'one'
p.yaxis.axis_label = 'two'

#################################
p2 = figure(tools="tap,reset",tooltips=hover.tooltips)
p2.circle(x='one', y='two',
         source=source3,
         size=2, color='color', legend='clsLabel')

################################################
N, M=160*2, 800
kw = dict()
kw['x_range'] = (0,M)
kw['y_range'] = (0,N)
kw['plot_width'] = M
kw['plot_height'] = N
kw['title']='               Protein                        Mito                         \
    ER                         DNA                          Segmentation        '
ts1 = figure(tools='pan,wheel_zoom,xbox_select,reset', **kw)
ts1.axis.visible = False
ts1.xgrid.visible = False
ts1.ygrid.visible = False
#        print(images[l].shape,Ms[l],Ns[l])
color = LinearColorMapper(bokeh.palettes.gray(256))
r1=ts1.image(image=[], x=0, y=0, dw=M, dh=N, color_mapper=color)

def get_images(imInf):
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
# print(hover.tooltips)
# def read_rgb(imagePath):
#     boxSize=200;
#     xCenter,yCenter=200,200
    
#     img = np.squeeze(skimage.io.imread(imagePath))[yCenter-boxSize:yCenter+boxSize,xCenter-boxSize:xCenter+boxSize]
# #     img = skimage.io.imread(imagePath)
#     print(img.shape)
# #     img=np.flip(img,0)
# #     rgb_img = img.astype(np.uint8)
#     return img

def create_figure(imageInfoList):
    
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
    print(old,new)
    tsneResDF2,data4umap2=load_data();
    source2.data=ColumnDataSource(data4umap2).data
#     print('shapeHere',tsneResDF2.shape,data4umap2.shape)
    colors = [colormap[x] for x in tsneResDF2['Label']]
    d2=dict(one=tsneResDF2['one'],
        two=tsneResDF2['two'],
        color=colors,
        Label=tsneResDF2['Label'],
        clsLabel=tsneResDF2['clsLabel'],
        Metadata_Plate=tsneResDF2['Metadata_Plate'],
        Metadata_Well=tsneResDF2['Metadata_Well'],
        Metadata_FieldID=tsneResDF2['Metadata_FieldID'],
        #         actual_image_number=sample['actual_image_number'],
        ObjectNumber=tsneResDF2['ObjectNumber'],
        diff=tsneResDF2['diffMinMaxNucCyto'],
        meanInt=tsneResDF2['meanInt'],
        UpperQ=tsneResDF2['UpperQ'],
        MaxInt=tsneResDF2['MaxInt'],
        StdInt=tsneResDF2['StdInt'])
    source.data=d2
    
    dCls2=d2.copy()
    palette2 = viridis(len(tsneResDF2['clsLabel'].unique()))
    colors2 = [palette2[int(x)] for x in tsneResDF2['clsLabel']]
#     dCls2['color']=color_map
    dCls2=dict(one=tsneResDF2['one'],
        two=tsneResDF2['two'],
        color=colors2,
        Label=tsneResDF2['Label'],
        clsLabel=tsneResDF2['clsLabel'],
        Metadata_Plate=tsneResDF2['Metadata_Plate'],
        Metadata_Well=tsneResDF2['Metadata_Well'],
        Metadata_FieldID=tsneResDF2['Metadata_FieldID'],
        #         actual_image_number=sample['actual_image_number'],
        ObjectNumber=tsneResDF2['ObjectNumber'],
        diff=tsneResDF2['diffMinMaxNucCyto'],
        meanInt=tsneResDF2['meanInt'],
        UpperQ=tsneResDF2['UpperQ'],
        MaxInt=tsneResDF2['MaxInt'],
        StdInt=tsneResDF2['StdInt'])
    
    source3.data=dCls2
    
    measured1=tsneResDF2['meanInt'].values;
    hist1, edges1 = np.histogram(measured1, density=True, bins=50)
    kde1 = stats.gaussian_kde(measured1.astype(float))

    dDist=dict(hist_1=hist1,
               edges_1_right=edges1[1:],
               edges_1_left=edges1[:-1])
#                kde_1=kde1)  

    sourceDis.data=dDist
    mem_stream.send(tsneResDF2)
#     make_plot_dist()
#     p1hv = hv.BoxWhisker(tsneResDF2, ['clsLabel'], 'meanInt', label="mean intensity per cluster").sort()
#     p1hv.opts(show_legend=False, width=500, box_fill_color=dim('origin').str(), cmap='Set1')
#     p1 = hv.render(p1hv)
#     p1.refresh()
    update()    
    
# def update_notes(attrname, old, new):
#     t1, t2 = ticker1.value, ticker2.value.split(' _')[0]
#     notesTxtFile=DATA_DIR+'/WT-MUT-'+t1+'_'+t2+'/observerNotes.txt'
#     fileW = open(notesTxtFile,"w") 
# #    print(new,inputText.value)
# #    inputText.value=str(new)
#     fileW.write(new)  
#     fileW.close()
    
# def update(selected=None):
#     print(hover.dataspecs())
#     print('hello')
#     fS=create_figure()
#     pp=make_plot()
#     t1, t2 = ticker1.value, ticker2.value.split(' _')[0]
# #    print(t1,t2)
#     stats2.text = str(load_ticker_info(t1,t2))
#     inputText.value=str(load_observerNotes(t1,t2))
#     series.children=[column(widgets, pp)]+[column(row(stats2,inputText),fS)]
    
# inputText.on_change('value', update_notes)
# hover.on_change('tooltips', ticker0_change)

def update():
#     print('source.selected.indices',source.selected.indices,len(source.data),len(source.data['Metadata_Plate']))
#     print(source.selected)
    indexActive=source.selected.indices
#     indexActive2=source3.selected.indices
#     indexActive=indexActive1+indexActive2
#     print('indexActive',indexActive,indexActive1,indexActive2)
    
    imageInfoList=[source.data['Metadata_Plate'][indexActive[0]],source.data['Metadata_Well'][indexActive[0]],source.data['Metadata_FieldID'][indexActive[0]],source.data['ObjectNumber'][indexActive[0]]]
#     print('imageInfoList',imageInfoList)
    fS=create_figure(imageInfoList)
#     series.children=[ticker1]+[row(p,p2,p1)]+[fS]
    series.children=[ticker1]+[row(p,p2,column(stats2,p1))]+[fS]
#     print(source.data['Label'][indexActive])
#     print(source.data['ObjectNumber'][indexActive])

def update2():
#     print('source.selected.indices',source.selected.indices,len(source.data),len(source.data['Metadata_Plate']))
#     print(source.selected)
    indexActive=source3.selected.indices
    imageInfoList=[source3.data['Metadata_Plate'][indexActive[0]],source3.data['Metadata_Well'][indexActive[0]],source3.data['Metadata_FieldID'][indexActive[0]],source3.data['ObjectNumber'][indexActive[0]]]
#     print('imageInfoList',imageInfoList)
    fS=create_figure(imageInfoList)
#     series.children=[ticker1]+[row(p,p2,p1)]+[fS]
    series.children=[ticker1]+[row(p,p2,column(stats2,p1))]+[fS]

    
    
    
p.on_event(Tap, update)
p2.on_event(Tap, update2)

ticker1.on_change('value', ticker1_change)

# indexActive=source.selected.indices[0]
indd=0
imageInfoList0=[d['Metadata_Plate'][indd],d['Metadata_Well'][indd],d['Metadata_FieldID'][indd],d['ObjectNumber'][indd]]

fS=create_figure(imageInfoList0)

p.add_tools(hover)
p2.add_tools(hover)

# plotss = row(p,p2,p1)
plotss = row(p,p2,column(stats2,p1))
series = column(children=[ticker1]+[plotss]+[fS])
# series2 = row(children=[widgets]+[fS])
# show(series)
# update()

curdoc().add_root(series)

# show(series)





