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
rootDir='/home/ubuntu/bucket/projects/2017_10_19_Profiling_rare_ORFs/workspace'
inputWellsDf=pd.DataFrame(os.listdir(rootDir+'/backend/plate_well_dfs'))
inputWellsDf['Metadata_Plate']=inputWellsDf[0].apply(lambda x: '_'.join(str(x).split('_')[1:-1]))
inputWellsDf['Metadata_Well']=inputWellsDf[0].apply(lambda x: str(x).split('_')[-1])

#############################################

# cpFeatures=dfWithWTlabels.columns[dfWithWTlabels.columns.str.contains("Cells_|Cytoplasm_|Nuclei_")]# & (perWellDataUntransW.columns.str.contains("Intensity"))]

# # cpFeatures=list(set(cpFeatures)-set(blackListFeatures))
# locFeature2beremoved=list(filter(lambda x: "_Location_Center_X" in x or "_Location_Center_Y" in x , cpFeatures)) 
# cpFeatures4scale=list(set(cpFeatures)-set(locFeature2beremoved))


# scaler = preprocessing.RobustScaler()
# dataScaled=scaler.fit_transform(dfWithWTlabels.loc[:,cpFeatures4scale])
# dfWithWTlabels_scaled = dfWithWTlabels.copy()
# dfWithWTlabels_scaled[cpFeatures4scale]=dataScaled


############################# form ticker 1
inputWellsDf['ticker1']=inputWellsDf['Metadata_Plate'].astype(str)+'__'+inputWellsDf['Metadata_Well'].astype(str)#+'__'+dfRank4bokeh1['Metadata_Location']


DEFAULT_TICKERS=sorted(inputWellsDf['ticker1'].tolist())[10:]

ticker1 = Select(title="WT-MT pairs", value=DEFAULT_TICKERS[0], options=DEFAULT_TICKERS)

    

############################ load single cell data
def load_data():  # get well data
    tickVal = ticker1.value
    plate=tickVal.split('__')[0]
    well=tickVal.split('__')[1]
    
    plateWellFilePath=rootDir+'/backend/plate_well_dfs/df_'+plate+'_'+well
    wellData=pd.read_pickle(plateWellFilePath, compression='infer'); 
    
    wellMaxPixelInt_P=wellData['Cells_Intensity_MaxIntensity_Protein'].values.max()
    wellMaxPixelInt_D=wellData['Cells_Intensity_MaxIntensity_DNA'].values.max()
    wellMaxPixelInt_E=wellData['Cells_Intensity_MaxIntensity_ER'].values.max()
    wellMaxPixelInt_M=wellData['Cells_Intensity_MaxIntensity_Mito'].values.max()
#     global wellMaxPix
#     wellMaxPix=[wellMaxPixelInt_P,wellMaxPixelInt_D,wellMaxPixelInt_E,wellMaxPixelInt_M]
#     wellData=dfWithWTlabels_scaled[dfWithWTlabels_scaled['Metadata_Well']==well]
    print(wellData['Metadata_Site'].unique())
    cpFeatures=wellData.columns[(wellData.columns.str.contains("_Protein")) & wellData.columns.str.contains("Cells_|Cytoplasm_|Nuclei_") & (wellData.columns.str.contains("Intensity"))]
#     cpFeaturesM=perWellDataUntransM.columns[(perWellDataUntransM.columns.str.contains("_Protein")) & perWellDataUntransM.columns.str.contains("Cells_|Cytoplasm_|Nuclei_")]# & (perWellDataUntransM.columns.str.contains("Intensity"))]
    
#     cpFeatures=wellData.columns[(wellData.columns.str.contains("_Protein"))]
    
    locFeature2beremoved=list(filter(lambda x: "_Location_Center_X" in x or "_Location_Center_Y" in x , cpFeatures)) 
    corFeature2beremoved=list(filter(lambda x: "Correlation" in x , cpFeatures)) 

    allFeatures=list(set(cpFeatures)-set(locFeature2beremoved)-set(corFeature2beremoved))

    
    scaler = preprocessing.RobustScaler()
    dataScaled=scaler.fit_transform(wellData.loc[:,allFeatures])
    wellData_scaled = wellData.copy()
    wellData_scaled[allFeatures]=dataScaled

    
    
    import umap
    from sklearn.manifold import TSNE
    umapT=umap.UMAP()
    tsneT = TSNE(perplexity=10)
#     data4umap=pd.concat([perWellDataFilteredW, perWellDataFilteredM,perWellDataUntransW, perWellDataUntransM], ignore_index=True,sort=False)
    
    data4umap=pd.concat([wellData_scaled], ignore_index=True,sort=False)
#     data4umap['Label']='untransCellsW'
    data4umap['Label']=data4umap['Metadata_Site'].astype(str)
#     data4umap=pd.concat([scaledWellW,scaledWellM], ignore_index=True,sort=False)
#     data4umap=pd.concat([scaledCellsW,scaledCellsM,scaledUntransCellsW, scaledUntransCellsM], ignore_index=True,sort=False)
#     data4umap['diffMinMaxNucCyto']=abs(data4umap['Cells_Intensity_MaxIntensity_Protein_corrected'])
    data4umap['meanInt']=wellData['Cells_Intensity_MeanIntensity_Protein']
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
#     print(clusterLabels.shape,preprocData.shape,preprocData.min(),preprocData.max())
    
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
#     print(data4umap['Label'].unique())
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

    ch_seg=projectPath+'/workspace/analysis/'+batch+'/'+plateName+'/analysis/'+plateName+'-'+well0+'-'+str(int(field))+\
    '/Cell_Outlines/'+well0+'_s'+str(int(field))+'_cell_outlines.tiff'  ## for Varsha's data
    
#     ch_seg=projectPath+'/workspace/analysis/'+batch+'/'+plateName+'/analysis/'+plateName+'-'+well0+\
#     '/Cell_outlines/'+well0+'_s'+str(int(field))+'_cell_outlines.png'    
    
    dataFileName=projectPath+batch+'/images/'+plateName+'/'   # for dataset1 and Varsha's pilots
#     dataFileName=projectPath+batch+'/images/'+plateName+'/Images/'   # for dataset2
#     dataFileName=projectPath+batch+'/'+plateName+'/'   # for varsha arrayed data  (changed this, now like dataset1)  
    
    
    boxSize=50;
    xCenter=int(dfWithWTlabels['Nuclei_Location_Center_X'].values[0])
    yCenter=int(dfWithWTlabels['Nuclei_Location_Center_Y'].values[0])
    print(dataFileName+ch_M,xCenter,yCenter)
    if (xCenter>boxSize) & (yCenter>boxSize) & (os.path.exists(dataFileName+ch_p)) & (os.path.exists(dataFileName+ch_M)) &\
    (os.path.exists(dataFileName+ch_E)) & (os.path.exists(dataFileName+ch_D)):
#         print(xCenter,yCenter,boxSize)
        imP=np.squeeze(skimage.io.imread(dataFileName+ch_p))
        mP=np.percentile(imP, 99)#imP.max()
        imP=imP[yCenter-boxSize:yCenter+boxSize,xCenter-boxSize:xCenter+boxSize]
        
        imM=np.squeeze(skimage.io.imread(dataFileName+ch_M))
        mM=np.percentile(imM, 99)#imP.max()
        imM=imM[yCenter-boxSize:yCenter+boxSize,xCenter-boxSize:xCenter+boxSize]
        
        imE=np.squeeze(skimage.io.imread(dataFileName+ch_E))        
        mE=np.percentile(imE, 99)#imP.max()
        imE=imE[yCenter-boxSize:yCenter+boxSize,xCenter-boxSize:xCenter+boxSize]
        
        imD=np.squeeze(skimage.io.imread(dataFileName+ch_D))
        mD=np.percentile(imD, 99)#imP.max()
        imD=imD[yCenter-boxSize:yCenter+boxSize,xCenter-boxSize:xCenter+boxSize]
        
        print(ch_seg)
        if (os.path.exists(ch_seg)):
            imSeg=np.squeeze(skimage.io.imread(ch_seg))[yCenter-boxSize:yCenter+boxSize,xCenter-boxSize:xCenter+boxSize]
            imSeg = rgb2gray(imSeg)
        else:
            imSeg=np.ones((2*boxSize,2*boxSize));
#         print(np.squeeze(skimage.io.imread(dataFileName+ch_p)).shape)
#         print(imP.max(),imM.max(),imE.max(),imD.max(),imSeg.max())
#         mP=15000;mM=1500;mE=1500;mD=1500
#         mP=imP.max();mM=imM.max();mE=imE.max();mD=imD.max()
        print('max',imP.max(),imM.max(),imE.max(),imD.max())
#         mP=wellMaxPix[0];mM=wellMaxPix[3];mE=wellMaxPix[2];mD=wellMaxPix[1]
#         print('max2',mP,mM,mE,mD)
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
#     print(data.max())
    p1hv = hv.BoxWhisker(data, ['clsLabel'], 'meanInt', label="mean intensity per cluster").sort()
    p1hv.opts(show_legend=False, width=500, cmap='Set1')#,ylim=(-1, 40))
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
N, M=160*3, 800
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
#     uniqLabels0=tsneResDF2['Label'].unique().tolist()
#     randIntArr=np.random.randint(1,len(uniqLabels0), size=3)
#     randSites=[uniqLabels0[i] for i in randIntArr]
#     tsneResDF2m=tsneResDF2[tsneResDF2['Label'].isin(randSites)].reset_index(drop=True)
    tsneResDF2m=tsneResDF2.copy()
    
    source2.data=dict(ColumnDataSource(data4umap2).data)    #### update source2
#     print('shapeHere',tsneResDF2.shape,data4umap2.shape)
    
    d2=gen_data_forHover(tsneResDF2m)
    uniqLabels=tsneResDF2m['Label'].unique().tolist()
#     palette1 = inferno(len(uniqLabels))
#     colors = [palette1[int(uniqLabels.index(x))] for x in tsneResDF2m['Label']]
#     palette1 = all_palettes['Category20c'][len(uniqLabels)]
    palette1 = viridis(len(uniqLabels))
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




