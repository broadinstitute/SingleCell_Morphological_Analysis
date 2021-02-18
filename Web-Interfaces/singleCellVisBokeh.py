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
from bokeh.events import Tap
# import utils.visualization




# annot_df_2['rep']
blackListFeatures=[
'Nuclei_Correlation_Manders_AGP_DNA',
'Nuclei_Correlation_Manders_AGP_ER',
'Nuclei_Correlation_Manders_AGP_Mito',
'Nuclei_Correlation_Manders_AGP_RNA',
'Nuclei_Correlation_Manders_DNA_AGP',
'Nuclei_Correlation_Manders_DNA_ER',
'Nuclei_Correlation_Manders_DNA_Mito',
'Nuclei_Correlation_Manders_DNA_RNA',
'Nuclei_Correlation_Manders_ER_AGP',
'Nuclei_Correlation_Manders_ER_DNA',
'Nuclei_Correlation_Manders_ER_Mito',
'Nuclei_Correlation_Manders_ER_RNA',
'Nuclei_Correlation_Manders_Mito_AGP',
'Nuclei_Correlation_Manders_Mito_DNA',
'Nuclei_Correlation_Manders_Mito_ER',
'Nuclei_Correlation_Manders_Mito_RNA',
'Nuclei_Correlation_Manders_RNA_AGP',
'Nuclei_Correlation_Manders_RNA_DNA',
'Nuclei_Correlation_Manders_RNA_ER',
'Nuclei_Correlation_Manders_RNA_Mito',
'Nuclei_Correlation_RWC_AGP_DNA',
'Nuclei_Correlation_RWC_AGP_ER',
'Nuclei_Correlation_RWC_AGP_Mito',
'Nuclei_Correlation_RWC_AGP_RNA',
'Nuclei_Correlation_RWC_DNA_AGP',
'Nuclei_Correlation_RWC_DNA_ER',
'Nuclei_Correlation_RWC_DNA_Mito',
'Nuclei_Correlation_RWC_DNA_RNA',
'Nuclei_Correlation_RWC_ER_AGP',
'Nuclei_Correlation_RWC_ER_DNA',
'Nuclei_Correlation_RWC_ER_Mito',
'Nuclei_Correlation_RWC_ER_RNA',
'Nuclei_Correlation_RWC_Mito_AGP',
'Nuclei_Correlation_RWC_Mito_DNA',
'Nuclei_Correlation_RWC_Mito_ER',
'Nuclei_Correlation_RWC_Mito_RNA',
'Nuclei_Correlation_RWC_RNA_AGP',
'Nuclei_Correlation_RWC_RNA_DNA',
'Nuclei_Correlation_RWC_RNA_ER',
'Nuclei_Correlation_RWC_RNA_Mito',
'Nuclei_Granularity_14_AGP',
'Nuclei_Granularity_14_DNA',
'Nuclei_Granularity_14_ER',
'Nuclei_Granularity_14_Mito',
'Nuclei_Granularity_14_RNA',
'Nuclei_Granularity_15_AGP',
'Nuclei_Granularity_15_DNA',
'Nuclei_Granularity_15_ER',
'Nuclei_Granularity_15_Mito',
'Nuclei_Granularity_15_RNA',
'Nuclei_Granularity_16_AGP',
'Nuclei_Granularity_16_DNA',
'Nuclei_Granularity_16_ER',
'Nuclei_Granularity_16_Mito',
'Nuclei_Granularity_16_RNA']





rootDir='/home/ubuntu/bucket/projects/2017_09_27_RareDiseases_Taipale/workspace'
dataset='Set2'
ndd='3'
profType='mean'
normalization=''; #'' for no normaliztion 'n_' for per untransfected cells normalization
scaleMeanProfilesForEachPlate=1;
############################# load annotaations 
AnnotSet2 = pd.read_excel(rootDir+'/metadata/Set2Annots20190801.xlsx', sheet_name=None) 
print(AnnotSet2.keys())
df_1=AnnotSet2['Replicate_Plates']
df_1['batch']='Maxproj_Replicates_Original_Screen'
df_2=AnnotSet2['Kinase_Mutants']
# df_2=df_2.drop([159])
df_2['batch']='Maxproj_Kinase_Plates'
df_3=AnnotSet2['Common_Variants']
df_3['batch']='Maxproj_Common_Variants'
df_4=AnnotSet2['Cancer_Mutants']
df_4['batch']='Maxproj_Cancer_Mutations_Screen'
annot_df_2 = pd.concat([df_1,df_2,df_3,df_4],axis=0,ignore_index=True)
metaDataPlates=annot_df_2['Metadata_Plate'].unique()
# listOfPlates0=os.listdir(rootDir+'/backend/wellsSingleCells/')[1:]
listOfPlates0=os.listdir(rootDir+'/backend/'+profType+'PerWells'+ndd+'/')
annot_df_2['Metadata_Sample']=annot_df_2['Metadata_Sample'].str.rstrip()

######################### load mean profiles 
# listOfPlatesAndWells=[p for p in listOfPlates0 if p.split('df_n_')[1:][0][0:-4] in metaDataPlates]
strProf='df_'+normalization
listOfPlates=[p for p in listOfPlates0 if p.split(strProf)[1:] in metaDataPlates]
# listOfPlates.remove('df_n_Replicate_21')
scaler_Plate= preprocessing.StandardScaler()

df = pd.DataFrame();
for p in listOfPlates: #[0:1]: ['df_RC4_IF_05']:
    fileNameToSave=rootDir+'/backend/'+profType+'PerWells'+ndd+'/'+p;
    transfectedMeanPerWell=pd.read_pickle(fileNameToSave, compression='infer'); 
    if scaleMeanProfilesForEachPlate:
        dataScaled=scaler_Plate.fit_transform(transfectedMeanPerWell.loc[:,transfectedMeanPerWell.columns.str.contains("Cells_|Cytoplasm_|Nuclei_")])
        transfectedMeanPerWell.loc[:,transfectedMeanPerWell.columns.str.contains("Cells_|Cytoplasm_|Nuclei_")]=dataScaled
    df=df.append(transfectedMeanPerWell, ignore_index=True,sort='True')  
    
    
############################ merge mean profiles and anotations
annot_df_2['Gene']=annot_df_2['Metadata_Sample'].apply(lambda x: str(x).split(' ')[0])
annot_df_2['Class']=annot_df_2['Metadata_Sample'].apply(lambda x: 'Mutant' if len(str(x).split(' '))>1 else 'Wild-type')
# annot_df_2['Primary location']=annot_df_2['Metadata_Location'].apply(lambda x: str(x).split(',')[0])


annot_df_2=annot_df_2.rename(columns={"Metadata_Sample":"UNIQUE"})
# annot_df_2=annot_df_2.rename(columns={"Metadata_Sample":"UNIQUE"})

if 0:
    df=df[df['n_transf']>10]

cols2remove=[i for i in df.columns.tolist() if df[i].isnull().sum(axis=0)>.05*df.shape[0]]
print(cols2remove)
df2=df.drop(cols2remove, axis=1);
# cols2remove=[i for i in df2.columns.tolist() if df2[i].isnull().sum(axis=0)>0]
# df2 = df2.fillna(df2.median())

df2 = df2.interpolate()
if 'Plate' in df2.columns:
    df2=df2.rename(columns={"Plate":"Metadata_Plate"})
df2_2=pd.merge(df2, annot_df_2, how='inner',on=['Metadata_Plate','Metadata_Well']);

scaler = preprocessing.StandardScaler()
df_scaled = df2.copy()

# dataScaled=scaler.fit_transform(df2.loc[:,df2.columns.str.contains("Cells_|Cytoplasm_|Nuclei_")])
# df_scaled[df2.columns[df2.columns.str.contains("Cells_|Cytoplasm_|Nuclei_")].tolist()]=dataScaled

df4=pd.merge(df_scaled, annot_df_2, how='inner',on=['Metadata_Plate','Metadata_Well']);
df4=df4.rename(columns={"Metadata_Location":"manual_Annot","rep_x":"rep"})

df4=df4.assign(Metadata_transfection=df4['manual_Annot'].str.contains('expressed'))
df4=df4[df4['Metadata_transfection']==False]

df_meanP=df4.copy()
df_meanP['Gene']=df_meanP['Gene'].str.rstrip()
df_meanP['UNIQUE']=df_meanP['UNIQUE'].str.rstrip()

locFeature2beremoved=df_meanP.columns[df_meanP.columns.str.contains("_Location_Center_X|_Location_Center_Y")]

df_meanP=df_meanP[list(set(df_meanP.columns.tolist())-set(blackListFeatures)-set(locFeature2beremoved))]




df_meanP2=df_meanP.copy().reset_index(drop=True)

# df_meanP2=df_meanP2.replace('cytoplasm', 'er')


# dfWithProtFeatures_meanP=df_meanP.loc[:,df_meanP.columns.str.contains("Protein") & df_meanP.columns.str.contains("Cells|Cytoplasm|Nuclei")];
# dfWithAllFeatures_meanP=df_meanP.loc[:,df_meanP.columns.str.contains("Cells|Cytoplasm|Nuclei")];

dfWithProtFeatures_meanP=df_meanP2.loc[:,df_meanP2.columns.str.contains("Protein") & df_meanP2.columns.str.contains("Cells_|Cytoplasm_|Nuclei_")];
dfWithOutProtFeatures_meanP=df_meanP2.loc[:,~(df_meanP2.columns.str.contains("Protein")) & df_meanP2.columns.str.contains("Cells_|Cytoplasm_|Nuclei_")];
dfWithAllFeatures_meanP=df_meanP2.loc[:,df_meanP2.columns.str.contains("Cells_|Cytoplasm_|Nuclei_")];


scaler_P = preprocessing.StandardScaler()
scaler_P.fit(dfWithProtFeatures_meanP)
columns4scalerP=dfWithProtFeatures_meanP.columns;

scaler_NP = preprocessing.StandardScaler()
scaler_NP.fit(dfWithOutProtFeatures_meanP)
columns4scalerNP=dfWithOutProtFeatures_meanP.columns;

scaler_A = preprocessing.StandardScaler()
scaler_A.fit(dfWithAllFeatures_meanP)
columns4scalerA=dfWithAllFeatures_meanP.columns;

dataScaled_p=scaler_P.transform(dfWithProtFeatures_meanP)
dataScaled_np=scaler_NP.transform(dfWithOutProtFeatures_meanP)
dataScaled_a=scaler_A.transform(dfWithAllFeatures_meanP)


listOfPlates_2=annot_df_2['Metadata_Plate'].unique()
ndd='3'
# df_trans=pd.DataFrame(columns=['Metadata_Plate','Metadata_Well','nT','nU','Problem','Metadata_Efficiency','Metadata_Location'])
# for p in listOfPlates_2:
#     print(p)
# wells=annot_df_2.loc[annot_df_2['Metadata_Plate']==p,'Metadata_Well'].tolist()
#     df_agg_plate=pd.DataFrame();#pd.DataFrame(index=range(len(wells)))
#     for wi in range(len(wells)):
#     for wi in range(1):        
plateW,wellW='Kinase_Mutants_1_Replicate','D11'
plateM,wellM='Kinase_Mutants_1','E10'

# plateW,wellW='Replicate_24','F09'
# plateM,wellM='Replicate_24','G03'

annotForWellW=annot_df_2.loc[(annot_df_2['Metadata_Plate']==plateW) & (annot_df_2['Metadata_Well']==wellW),:].reset_index(drop=True)
fileNameW=rootDir+'/backend/wellsSingleCells'+ndd+'/df_'+plateW+'_'+wellW;
fileNameUW=rootDir+'/backend/wellsSingleCells'+ndd+'/df_control_'+plateW+'_'+wellW;

annotForWellM=annot_df_2.loc[(annot_df_2['Metadata_Plate']==plateM) & (annot_df_2['Metadata_Well']==wellM),:].reset_index(drop=True)
fileNameM=rootDir+'/backend/wellsSingleCells'+ndd+'/df_'+plateM+'_'+wellM;
fileNameUM=rootDir+'/backend/wellsSingleCells'+ndd+'/df_control_'+plateM+'_'+wellM;

perWellDataFilteredW=pd.read_pickle(fileNameW).reset_index(drop=True);  
perWellDataUntransW=pd.read_pickle(fileNameUW).reset_index(drop=True);        

perWellDataFilteredM=pd.read_pickle(fileNameM).reset_index(drop=True);  
perWellDataUntransM=pd.read_pickle(fileNameUM).reset_index(drop=True);      

print(perWellDataFilteredW.shape,perWellDataUntransW.shape,perWellDataFilteredM.shape,perWellDataUntransM.shape)
cpFeaturesW=perWellDataUntransW.columns[perWellDataUntransW.columns.str.contains("Cells_|Cytoplasm_|Nuclei_")]
cpFeaturesM=perWellDataUntransM.columns[perWellDataUntransM.columns.str.contains("Cells_|Cytoplasm_|Nuclei_")]
cpFeatures=list(set(cpFeaturesW).intersection(cpFeaturesM).intersection(columns4scalerA))

allFeatures=list(set(cpFeatures)-set(blackListFeatures))

nomalizedCellsW=perWellDataFilteredW.copy()
nomalizedCellsM=perWellDataFilteredM.copy()
nomalizedCellsW['Label']='nomalizedCellsW'
nomalizedCellsM['Label']='nomalizedCellsM'

scaledCellsW=perWellDataFilteredW.copy()
scaledCellsM=perWellDataFilteredM.copy()
scaledCellsW['Label']='scaledCellsW'
scaledCellsM['Label']='scaledCellsM'

scaledUntransCellsW=perWellDataUntransW.copy()
scaledUntransCellsM=perWellDataUntransM.copy()
scaledUntransCellsW['Label']='scaledUntransCellsW'
scaledUntransCellsM['Label']='scaledUntransCellsM'


scalerW = preprocessing.StandardScaler()
scalerM = preprocessing.StandardScaler()
scaler = preprocessing.StandardScaler()

scalerW.fit(perWellDataUntransW.loc[:,allFeatures].astype('float64'))
scalerM.fit(perWellDataUntransM.loc[:,allFeatures].astype('float64'))

nomalizedCellsW.loc[:,allFeatures]=scalerW.transform(perWellDataFilteredW.loc[:,allFeatures].astype('float64'))
nomalizedCellsM.loc[:,allFeatures]=scalerM.transform(perWellDataFilteredM.loc[:,allFeatures].astype('float64'))

scaledCellsW.loc[:,allFeatures]=scaler.fit_transform(perWellDataFilteredW.loc[:,allFeatures].astype('float64'))
scaledCellsM.loc[:,allFeatures]=scaler.fit_transform(perWellDataFilteredM.loc[:,allFeatures].astype('float64'))

scaledUntransCellsW.loc[:,allFeatures]=scaler.fit_transform(perWellDataUntransW.loc[:,allFeatures].astype('float64'))
scaledUntransCellsM.loc[:,allFeatures]=scaler.fit_transform(perWellDataUntransM.loc[:,allFeatures].astype('float64'))

scaledCellsMeanW=scaledCellsW.loc[:,allFeatures].mean().to_frame().T
scaledUnCellsMeanW=scaledUntransCellsW.loc[:,allFeatures].mean().to_frame().T
nomalizedCellsMeanW=nomalizedCellsW.loc[:,allFeatures].mean().to_frame().T
avgProfileW=pd.concat([nomalizedCellsMeanW,annotForWellW], axis=1)
MeanW=perWellDataFilteredW.loc[:,allFeatures].astype('float64').mean().to_frame().T

scaledCellsMeanM=scaledCellsM.loc[:,allFeatures].mean().to_frame().T
scaledUnCellsMeanM=scaledUntransCellsM.loc[:,allFeatures].mean().to_frame().T
nomalizedCellsMeanM=nomalizedCellsM.loc[:,allFeatures].mean().to_frame().T
avgProfileM=pd.concat([nomalizedCellsMeanM,annotForWellM], axis=1)
MeanM=perWellDataFilteredM.loc[:,allFeatures].astype('float64').mean().to_frame().T

avgPcc=scipy.stats.pearsonr(MeanW.T, MeanM.T)[0][0]
avgUnPcc=scipy.stats.pearsonr(perWellDataUntransW.loc[:,allFeatures].astype('float64').mean().to_frame(),\
                              perWellDataUntransM.loc[:,allFeatures].astype('float64').mean().to_frame())[0][0]
scalPcc=scipy.stats.pearsonr(scaledCellsMeanW.T, scaledCellsMeanM.T)[0][0]
scalUnPcc=scipy.stats.pearsonr(scaledUnCellsMeanW.T, scaledUnCellsMeanM.T)[0][0]
normPcc=scipy.stats.pearsonr(nomalizedCellsMeanW.T, nomalizedCellsMeanM.T)[0][0]

print('Avg profiles CC: ',avgPcc)
print('Avg untransfected profiles CC: ',avgUnPcc)
print('Normalized profiles CC: ',normPcc)
print('Scaled profiles CC: ',scalPcc)
print('Scaled untransfected profiles CC: ',scalUnPcc)


import umap
umapT=umap.UMAP()
data4umap=pd.concat([scaledCellsW,scaledCellsM,scaledUntransCellsW, scaledUntransCellsM], ignore_index=True,sort=False)
Y = umapT.fit_transform(data4umap.loc[:,allFeatures])
tsneResDF=pd.DataFrame(index=range(Y.shape[0]),columns=['one','two','Label','Metadata_Plate','Metadata_Well',\
                                                        'Metadata_FieldID','ObjectNumber']);
tsneResDF.loc[:,['one','two']]=Y
tsneResDF.loc[:,['Label','Metadata_Plate','Metadata_Well','Metadata_FieldID','ObjectNumber']]=\
data4umap[['Label','Metadata_Plate','Metadata_Well','Metadata_FieldID','ObjectNumber']]
g=sns.scatterplot(x="one", y="two", hue="Label", data=tsneResDF)



# plate0,well0=plateM,wellM
# plate0,well0=plateW,wellW
# field=21; objectN=10;
# field=5; objectN=16;
# field=3; objectN=55;
# field=15; objectN=38;
# field=11; objectN=40;


from bokeh.plotting import figure, output_file, show
from bokeh.models import ColumnDataSource
from bokeh.layouts import row, column, layout, gridplot
from bokeh.models.tools import HoverTool
from bokeh.io import output_notebook
import skimage.io
from bokeh.io import curdoc
# output_notebook()
# output_file('columndatasource_example.html')

# output_file("color_scatter.html", title="color_scatter.py example", mode="cdn")
# df = pd.read_csv('thor_wwii.csv')

# sample = tsneResDF.sample(100)
sample = tsneResDF.copy()
colormap = {'scaledUntransCellsW': 'blue', 'scaledUntransCellsM': 'green', 'scaledCellsW': 'red','scaledCellsM':'orchid'}
colors = [colormap[x] for x in sample['Label']]


# source = ColumnDataSource(sample)
# source = ColumnDataSource(sample)
d=dict(one=sample['one'],
        two=sample['two'],
        color=colors,
        Label=sample['Label'],
        Metadata_Plate=sample['Metadata_Plate'],
        Metadata_Well=sample['Metadata_Well'],
        Metadata_FieldID=sample['Metadata_FieldID'],
        #         actual_image_number=sample['actual_image_number'],
        ObjectNumber=sample['ObjectNumber']
    )
source = ColumnDataSource(data=d)


hover = HoverTool()
hover.tooltips=[
    ('Label', '@Label'),
    ('Metadata_Plate', '@Metadata_Plate'),    
    ('Metadata_Well', '@Metadata_Well'),    
    ('Metadata_FieldID', '@Metadata_FieldID'),    
    ('ObjectNumber', '@ObjectNumber'),    
]


p = figure(tools="tap,reset",tooltips=hover.tooltips)
p.circle(x='one', y='two',
         source=source,
         size=2, color='color')

p.title.text = 'TSNE - Applied on Mito features of single cells'
p.xaxis.axis_label = 'one'
p.yaxis.axis_label = 'two'


def get_images(imInf):
    plate0,well0,field,objectN=imInf[0],imInf[1],imInf[2],imInf[3]
    dfWithWTlabels=data4umap[(data4umap['Metadata_Plate']==plate0) & (data4umap['Metadata_Well']==well0) &\
                 (data4umap['Metadata_FieldID']==field) & (data4umap['ObjectNumber']==objectN)]
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
#         dataFileName=projectPath+batch+'/images/'+plateName+'/'   # for dataset1
    dataFileName=projectPath+batch+'/images/'+plateName+'/Images/'   # for dataset2

    boxSize=100;
    xCenter=int(dfWithWTlabels['Nuclei_Location_Center_X'].values[0])
    yCenter=int(dfWithWTlabels['Nuclei_Location_Center_Y'].values[0])
    print(dataFileName+ch_M,xCenter,yCenter)
    if os.path.exists(dataFileName+ch_p) & os.path.exists(dataFileName+ch_M) &\
    os.path.exists(dataFileName+ch_E) & os.path.exists(dataFileName+ch_D):
        imP=np.squeeze(skimage.io.imread(dataFileName+ch_p))[yCenter-boxSize:yCenter+boxSize,xCenter-boxSize:xCenter+boxSize]
        imM=np.squeeze(skimage.io.imread(dataFileName+ch_M))[yCenter-boxSize:yCenter+boxSize,xCenter-boxSize:xCenter+boxSize]
        imE=np.squeeze(skimage.io.imread(dataFileName+ch_E))[yCenter-boxSize:yCenter+boxSize,xCenter-boxSize:xCenter+boxSize]
        imD=np.squeeze(skimage.io.imread(dataFileName+ch_D))[yCenter-boxSize:yCenter+boxSize,xCenter-boxSize:xCenter+boxSize]
    
    return np.concatenate([imP/(imP.max()),imM/(imM.max()),imE/(imE.max()),imD/(imD.max())],axis=1)

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
    adss='/home/ubuntu/bucket/projects/2017_09_27_RareDiseases_Taipale/Maxproj_Kinase_Plates/images/Kinase_Mutants_1/Images/r05c10f11p02-ch2sk1fk1fl1.tiff'
#     adss='/home/ubuntu/2017_09_27_RareDiseases_Taipale/cluster10_examplar.png'
    images = get_images(imageInfoList)
    N, M=200, 800
    kw = dict()
    kw['x_range'] = (0,M)
    kw['y_range'] = (0,N)
    kw['plot_width'] = M
    kw['plot_height'] = N
    ts1 = figure(tools='pan,wheel_zoom,xbox_select,reset', **kw)
    ts1.axis.visible = False
    ts1.xgrid.visible = False
    ts1.ygrid.visible = False
#        print(images[l].shape,Ms[l],Ns[l])
    r1=ts1.image(image=[], x=0, y=0, dw=M, dh=N)
    ds1 = r1.data_source;
    new_data1['image'] = [images]
    ds1.data = new_data1
#     ts.append(ts1);
    return ts1

def ticker0_change(attrname, old, new):
    print(old,new)
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
    print(source.selected.indices)
#     print(source.selected)
    indexActive=source.selected.indices[0]
    imageInfoList=[d['Metadata_Plate'][indexActive],d['Metadata_Well'][indexActive],d['Metadata_FieldID'][indexActive],d['ObjectNumber'][indexActive]]
    fS=create_figure(imageInfoList)
    series.children=[p]+[fS]
    print(d['Label'][indexActive])
    print(d['ObjectNumber'][indexActive])


p.on_event(Tap, update)


# indexActive=source.selected.indices[0]
indd=1
imageInfoList0=[d['Metadata_Plate'][indd],d['Metadata_Well'][indd],d['Metadata_FieldID'][indd],d['ObjectNumber'][indd]]

fS=create_figure(imageInfoList0)

p.add_tools(hover)
series = column(children=[p]+[fS])

# show(series)
# update()

curdoc().add_root(series)

# show(series)
