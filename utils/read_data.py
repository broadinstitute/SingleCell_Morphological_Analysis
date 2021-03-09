"""
This script contains usefull functions used in the notebooks

@author: mhaghigh
"""
import pandas as pd
import numpy as np
import seaborn as sns
from sklearn import preprocessing
import pickle
from sklearn import preprocessing
from sklearn.feature_selection import mutual_info_regression
from imblearn.over_sampling import SMOTE,RandomOverSampler
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score,confusion_matrix
import matplotlib.pyplot as plt
from imblearn.over_sampling import SMOTE,RandomOverSampler
import os
import time

import pandas as pd
from sqlalchemy import create_engine
from functools import reduce

def balanceLeaveOneCloneOutCV(model,dataScaled,labels,groupsNum):
    # 
    from sklearn.model_selection import LeaveOneGroupOut
    from sklearn.metrics import roc_auc_score
    acc=[]
    predLabels=[]; #np.zeros((labels.shape[0]), dtype=bool)
    trueLabels=[];
    predProbs=[];
    trueLabelConf=[];
    i=0;
    logo = LeaveOneGroupOut()
    acc2=[]
    featureImportanceListofLists=[]
    abundance=[]
    for train_index, test_index in logo.split(dataScaled, labels, groupsNum):
    #         print("TRAIN:", train_index, "TEST:", test_index)        
        leftOutClone=groupsNum[test_index]
        X_train, X_test = dataScaled[train_index], dataScaled[test_index]        
        y_train, y_test = labels[train_index], labels[test_index]
#         sm = SMOTE()
        sm=RandomOverSampler()
#         print(i,X_train.shape,y_train.shape)

        data_res_train, labels_res_train = sm.fit_sample(X_train, y_train)
#         print(i,data_res_train.shape,labels_res_train.shape)

        model.fit(data_res_train, labels_res_train) 
        importances = model.feature_importances_
        indicesRF = np.argsort(importances)[::-1]
        topXfs=indicesRF[0:10]
#         print(np.count_nonzero(labels_res_train==True),len(labels_res_train))
        y_model = model.predict(X_test);
        prob_model=model.predict_proba(X_test)
        accuracy=accuracy_score(y_test, y_model)
        trueLabelConf0=[prob_model[i][y_test[i]] for i in range(len(y_test))];
#         accuracy=roc_auc_score(y_test, model.predict_proba(X_test))
#         accuracy=model.score(X_test,y_test)

        if accuracy>.65:
            featureImportanceListofLists.append(topXfs)
        acc.append(accuracy);
        acc2.append(str(leftOutClone[0])+':'+str(np.round(accuracy*100,2)));
#         leftOtClone
        trueLabels+=y_test.tolist();
        predLabels+=y_model.tolist();
        trueLabelConf+=trueLabelConf0;
#         predProbs+=prob_model.tolist();
#         print(trueLabelConf0,(np.sum(np.array(trueLabelConf0)>=0.95),np.sum((np.array(trueLabelConf0)>=0.95) |(np.array(trueLabelConf0)<=0.05))))
#         abnd=(np.sum(np.array(trueLabelConf0)>=0.95)/np.sum((np.array(trueLabelConf0)>=0.95) |(np.array(trueLabelConf0)<=0.05)))
#         abundance+=[abnd];
    topXfeats = set(featureImportanceListofLists[0])
    for s in featureImportanceListofLists[1:]:
        topXfeats.intersection_update(s)
    return np.mean(acc),trueLabels,predLabels,acc2,trueLabelConf,topXfeats







def visualize_n_SingleCell(channels,dfWithWTlabels,boxSize):
    """ 
    This function plots the single cells correspoding to the input single cell dataframe
  
    Inputs: 
    ++ dfWithWTlabels   (pandas df) size --> (number of single cells)x(columns): 
    input dataframe contains single cells profiles as rows (make sure it has "Nuclei_Location_Center_X"or"Y" columns)
    
    ++ channels (dtype: list): list of channels to be displayed as columns of output image
           example: channels=['Mito','AGP','Brightfield','ER','DNA','Outline']
        * If Outline exist in the list of channels; function reads the outline image address from 
          "URL_CellOutlines" column of input dataframe, therefore, check that the addresses are correct
           before inputing them to the function, and if not, modify before input!
       
    boxSize (int): Height or Width of the square bounding box
    
    Returns: 
    f (object): handle to the figure
  
    """
    
    halfBoxSize=int(boxSize/2);
#     print(channels)
    
    import skimage.io
    f, axarr = plt.subplots(dfWithWTlabels.shape[0], len(channels),figsize=(len(channels)*2,dfWithWTlabels.shape[0]*2));
#     f.suptitle('Transfected: '+str(gfpTag))
    f.subplots_adjust(hspace=0, wspace=0)


    maxRanges={"DNA":8000,"RNA":6000,"Mito":6000,"ER":8000,"AGP":6000}
    for index in range(dfWithWTlabels.shape[0]):
        
        xCenter=int(dfWithWTlabels.loc[index,'Nuclei_Location_Center_X'])
        yCenter=int(dfWithWTlabels.loc[index,'Nuclei_Location_Center_Y'])
        
        cpi=0;
        for c in channels:
            if c=='Outline':
                imPath=dfWithWTlabels.loc[index,'URL_CellOutlines'];
            else:
                ch_D=dfWithWTlabels.loc[index,'Image_FileName_Orig'+c];
#                 print(ch_D)
    #         imageDir=imDir+subjectID+' Mito_Morphology/'
                imageDir=dfWithWTlabels.loc[index,'Image_PathName_Orig'+c]+'/'
                imPath=imageDir+ch_D
            
            imD=skimage.io.imread(imPath)[yCenter-halfBoxSize:yCenter+halfBoxSize,xCenter-halfBoxSize:xCenter+halfBoxSize]
            axarr[index,cpi].imshow(imD,cmap='gray',clim=(0, maxRanges[c]));axarr[0,cpi].set_title(c);
#             axarr[index,cpi].imshow(imD,cmap='gray');axarr[0,cpi].set_title(c);
            cpi+=1        

#         Well=dfWithWTlabels.loc[index,'Metadata_Well']
#         Site=str(dfWithWTlabels.loc[index,'Metadata_Site'])
#         imylabel=Well+'\n'+Site
#         axarr[index,0].set_ylabel(imylabel);            
            
            
#         subjectID=dfWithWTlabels.loc[index,'subject']
#         imylabel=dfWithWTlabels.label[index]+'\n'+subjectID
#         axarr[index,0].set_ylabel(imylabel);
# #     plt.tight_layout() 

    for i in range(len(channels)):
        for j in range(dfWithWTlabels.shape[0]):
            axarr[j,i].xaxis.set_major_locator(plt.NullLocator())
            axarr[j,i].yaxis.set_major_locator(plt.NullLocator())
            axarr[j,i].set_aspect('auto')
    
    return 

def readSingleCellData_r(fileName):
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
#     pd_df=edgeCellFilter(pd_df);  
#     print(4)
    return pd_df

########################function to remove cells on the border
def edgeCellFilter(df_1):   
    # remove cells on the border
#     imgSize=1080
#     imgSize=2048
    imgSize=dfWithWTlabels.Metadata_ImageSizeX.values[0]
    borderLength=int(np.percentile(dfWithWTlabels.Cells_AreaShape_MajorAxisLength.values, 90)/2);
    df_1_centCells=df_1.loc[~((df_1['Nuclei_Location_Center_X']>(imgSize-borderLength)) | \
                              (df_1['Nuclei_Location_Center_X']<(borderLength))\
                            | (df_1['Nuclei_Location_Center_Y']>(imgSize-borderLength)) | \
                              (df_1['Nuclei_Location_Center_Y']<(borderLength))),:].reset_index(drop=True)
    return df_1_centCells,borderLength



########################function to remove cells on the border
def edgeCellFilter2(df_1,imgSize,borderLength):   

    df_1_centCells=df_1.loc[~((df_1['Nuclei_Location_Center_X']>(imgSize-borderLength)) | \
                              (df_1['Nuclei_Location_Center_X']<(borderLength))\
                            | (df_1['Nuclei_Location_Center_Y']>(imgSize-borderLength)) | \
                              (df_1['Nuclei_Location_Center_Y']<(borderLength))),:].reset_index(drop=True)
    return df_1_centCells


# def edgeCellFilter(df_1):   
#     # remove cells on the border
# #     imgSize=1080
# #     imgSize=2048
#     print(df_1.columns[df_1.columns.str.contains("Metadata_ImageSizeX")])
#     print(df_1.columns[df_1.columns.str.contains("ImageSize")])
#     imgSize=df_1.Metadata_ImageSizeX.values[0]
#     borderLength=int(np.percentile(df_1.Cells_AreaShape_MajorAxisLength.values, 90)/2);
#     print(imgSize,borderLength)
#     df_1_centCells=df_1.loc[~((df_1['Nuclei_Location_Center_X']>(imgSize-borderLength)) | \
#                               (df_1['Nuclei_Location_Center_X']<(borderLength))\
#                             | (df_1['Nuclei_Location_Center_Y']>(imgSize-borderLength)) | \
#                               (df_1['Nuclei_Location_Center_Y']<(borderLength))),:].reset_index(drop=True)

#     return df_1_centCells,borderLength

def readSingleCellData_sqlalch(fileName,compartments):
    from sqlalchemy import create_engine
    sql_file="sqlite:////"+fileName
    engine = create_engine(sql_file)
    conn = engine.connect()
#     compartments=["cells", "cytoplasm", "nuclei"]
    # compartments=["Neurites","CellBodies","CellBodiesPlusNeurites","Nuclei","Cytoplasm"]
    plateDf_list=[]
    for compartment in compartments:
        compartment_query = "select * from {}".format(compartment)
        plateDf_list.append(pd.read_sql(sql=compartment_query, con=conn))

    plateDf = reduce(lambda left,right: pd.merge(left,right,on=["TableNumber", "ImageNumber", "ObjectNumber"]), plateDf_list)

    compartment_query = "select * from {}".format("Image")
    plateImageDf= pd.read_sql(sql=compartment_query, con=conn);

    plateDfwMeta = pd.merge(plateDf, plateImageDf, on=["TableNumber", "ImageNumber"])
    plateDfwMeta = plateDfwMeta.loc[:,~plateDfwMeta.columns.duplicated()]
    
    return plateDfwMeta



def readSingleCellData_sqlalch_random_image_subset(fileName,n_rand_ims):
    import pandas as pd
    from sqlalchemy import create_engine
    from functools import reduce
    # sql_file="sqlite:////"+rootDir+'backend/JUANCHI_PILOT_1/JUANCHI_PILOT_1_PLATE_1/JUANCHI_PILOT_1_PLATE_1.sqlite'
    # sql_file="sqlite:////"+rootDirDrug+"/backend/"+batchName+"/"+p+"/"+p+".sqlite"    

    sql_file="sqlite:////"+fileName
    engine = create_engine(sql_file)
    conn = engine.connect()
    compartments=["cells", "cytoplasm", "nuclei"]
    # compartments=["Neurites","CellBodies","CellBodiesPlusNeurites","Nuclei","Cytoplasm"]

    rand_img_num=np.random.choice(range(1,4000), n_rand_ims)
    list_str="("
    for i in rand_img_num:
        list_str=list_str+str(i)+',' 
    list_str=list_str[:-1]+")"

    plateDf_list=[]
    for compartment in compartments:
#         compartment_query = "select * from {}".format(compartment)
        compartment_query = "select * from {} WHERE {} IN {};".format(compartment,"ImageNumber",list_str)
        plateDf_list.append(pd.read_sql(sql=compartment_query, con=conn))

    plateDf = reduce(lambda left,right: pd.merge(left,right,on=["TableNumber", "ImageNumber", "ObjectNumber"]), plateDf_list)

#     compartment_query = "select * from {}".format("Image")
    compartment_query = "select * from {} WHERE {} IN {};".format("Image","ImageNumber",list_str)
    plateImageDf= pd.read_sql(sql=compartment_query, con=conn);

    plateDfwMeta = pd.merge(plateDf, plateImageDf, on=["TableNumber", "ImageNumber"])
    plateDfwMeta = plateDfwMeta.loc[:,~plateDfwMeta.columns.duplicated()]
    
    return plateDfwMeta   



def readSingleCellData_sqlalch_concat(fileName):
    from sqlalchemy import create_engine
    # sql_file="sqlite:////"+rootDir+'backend/JUANCHI_PILOT_1/JUANCHI_PILOT_1_PLATE_1/JUANCHI_PILOT_1_PLATE_1.sqlite'
    # sql_file="sqlite:////"+rootDirDrug+"/backend/"+batchName+"/"+p+"/"+p+".sqlite"    

    sql_file="sqlite:////"+fileName
    engine = create_engine(sql_file)
    conn = engine.connect()
    compartments=["cells", "cytoplasm", "nuclei"]
    #     compartments="cells, cytoplasm, nuclei"
    #     compartments="cytoplasm"
    # compartments=["Neurites","CellBodies","CellBodiesPlusNeurites","Nuclei","Cytoplasm"]
    #     plateDf_list=[]
    #     for compartment in compartments:
    #         compartment_query = "select * from {}".format(compartment)
    #         plateDf_list.append(pd.read_sql(sql=compartment_query, con=conn))

    #     plateDf = reduce(lambda left,right: pd.merge(left,right,on=["TableNumber", "ImageNumber", "ObjectNumber"]), plateDf_list)
    #     plateDf=pd.concat(plateDf_list, axis=1)
    #     compartment_query = "select * from {}".format(compartments[0])
    # query_cols = "TableNumber, ImageNumber, "+selected_feature

    ########################################## find common TableNumber+ImageNumber
    for c in compartments:
        query_cols = "TableNumber, ImageNumber, ObjectNumber"
        compartment_query = "select {} from {}".format(query_cols,c)
        compartTableDf0=pd.read_sql(sql=compartment_query, con=conn)
        listOfTandI0=(compartTableDf0["TableNumber"] +'_'+compartTableDf0["ImageNumber"].astype(str)+\
                     '_'+compartTableDf0["ObjectNumber"].astype(str)).tolist()
        if c==compartments[0]:
            listOfTandI=listOfTandI0.copy()
        else:
            listOfTandI=list(set(listOfTandI) & set(listOfTandI0)) 
    #         compartTableDf = compartTableDf.join(compartTableDf0,on=["TableNumber", "ImageNumber"], how='left')
    #         compartTableDf=pd.merge(compartTableDf0, compartTableDf,how='inner', on=["TableNumber", "ImageNumber"]);

#     print(len(listOfTandI))

    intersectionOfTableNs=pd.DataFrame(columns=["TableNumber","ImageNumber","ObjectNumber"])
    intersectionOfTableNs["TableNumber"]=[s.split('_')[0] for s in listOfTandI0]
    intersectionOfTableNs["ImageNumber"]=[int(s.split('_')[1]) for s in listOfTandI0]


    plateDfList=[]
    for c in compartments:
        compartment_query = "select * from {}".format(c)
        plateDf0=pd.read_sql(sql=compartment_query, con=conn)
        plateDf0["TableImage"]=plateDf0["TableNumber"] +'_'+plateDf0["ImageNumber"].astype(str)+'_'+plateDf0["ObjectNumber"].astype(str)
        plateDf0=plateDf0[plateDf0['TableImage'].isin(listOfTandI)].sort_values(["TableNumber", "ImageNumber","ObjectNumber"]).reset_index(drop=True)

        if c==compartments[0]:
            plateDfList.append(plateDf0)
            checkTanIdf=plateDf0[["TableNumber", "ImageNumber","ObjectNumber"]].sort_values(["TableNumber", "ImageNumber","ObjectNumber"]).reset_index(drop=True)
        else:        
            if plateDf0[["TableNumber", "ImageNumber","ObjectNumber"]].equals(checkTanIdf):
                plateDfList.append(plateDf0[list(set(plateDf0.columns)-set(["TableNumber", "ImageNumber","ObjectNumber"]))])
            else:
                print("change the code")

    plateDf=pd.concat(plateDfList, axis=1)            
#     print(plateDf.shape)         
    # del plateDf_list
    #     plateDf = plateDf.loc[:,["TableNumber", "ImageNumber",selected_feature]]
    # plateDf=plateDf.dropna()

    compartment_query = "select * from {}".format("Image")
    plateImageDf= pd.read_sql(sql=compartment_query, con=conn);

    plateDfwMeta = pd.merge(plateDf, plateImageDf, on=["TableNumber", "ImageNumber"])
    plateDfwMeta = plateDfwMeta.loc[:,~plateDfwMeta.columns.duplicated()]

    return plateDfwMeta



def readSingleCellData_sqlalch_well_subset(fileName,wells):
    import pandas as pd
    from sqlalchemy import create_engine
    from functools import reduce
    # sql_file="sqlite:////"+rootDir+'backend/JUANCHI_PILOT_1/JUANCHI_PILOT_1_PLATE_1/JUANCHI_PILOT_1_PLATE_1.sqlite'
    # sql_file="sqlite:////"+rootDirDrug+"/backend/"+batchName+"/"+p+"/"+p+".sqlite"    

#     wells=['A01', 'A02', 'A03', 'A04', 'A05', 'A06', 'A07', 'A08', 'A09',\
#        'A10', 'A11', 'A12', 'B01', 'B02', 'B03', 'B04', 'B05', 'B06',\
#        'B07', 'B08', 'B09', 'B10', 'B11', 'B12', 'C01', 'C02', 'C03',\
#        'C04', 'C05', 'C06', 'C07', 'C08', 'C09', 'C10', 'C11', 'C12',\
#        'D01', 'D02', 'D03', 'D04', 'D05', 'D06', 'D07', 'D08', 'D09',\
#        'D10', 'D11', 'D12', 'E01', 'E02', 'E03', 'E04', 'E05', 'E06',\
#        'E07', 'E08', 'E09', 'E10', 'E11', 'E12', 'F01', 'F02', 'F03',\
#        'F04', 'F05', 'F06', 'F07', 'F08', 'F09', 'F10', 'F11', 'F12',\
#        'G01', 'G02', 'G03', 'G04', 'G05', 'G06', 'G07', 'G08', 'G09',\
#        'G10', 'G11', 'G12', 'H01', 'H02', 'H03', 'H04', 'H05', 'H06',\
#        'H07', 'H08', 'H09', 'H10', 'H11', 'H12'];
#     meta_well_col_str="Image_Metadata_Well"
    meta_well_col_str="Metadata_Well"

    sql_file="sqlite:////"+fileName
    engine = create_engine(sql_file)
    conn = engine.connect()
    compartments=["cells", "cytoplasm", "nuclei"]
    # compartments=["Neurites","CellBodies","CellBodiesPlusNeurites","Nuclei","Cytoplasm"]

#     rand_img_num=np.random.choice(range(1,4000), n_rand_ims)
#     rand_img_num=np.array(range(50))
#     list_str="("
#     for i in rand_img_num:
#         list_str=list_str+str(i)+',' 
#     list_str=list_str[:-1]+")"
    
#     rand_img_num=wells[:40]
    rand_img_num=wells.copy()
    list_str="('"
    for i in rand_img_num:
        list_str=list_str+str(i)+"','" 
    list_str=list_str[:-2]+")"

#     compartment_query = "select * from {}".format("Image")
#     compartment_query = "select * from {} WHERE {} IN {};".format("Image","ImageNumber",list_str)
#     compartment_query = "select * from {} WHERE {} IN {};".\
    compartment_query = "select * from {} WHERE {} IN {};".\
    format("Image",meta_well_col_str,list_str)
    start1 = time.time()
    plateImageDf= pd.read_sql(sql=compartment_query, con=conn);
#     print(plateImageDf.columns[plateImageDf.columns.str.contains("Metadata_")])
#     print(plateImageDf['Metadata_Well'].unique())
    end1 = time.time()
    print('time elapsed:',(end1 - start1)/60)
    img_nums=plateImageDf.ImageNumber.unique().tolist()
    print(plateImageDf.shape,img_nums)
    list_str2="("
    for i in img_nums:
        list_str2=list_str2+str(i)+',' 
    list_str2=list_str2[:-1]+")"
    start2 = time.time()
    plateDf_list=[]
    for compartment in compartments:
#         compartment_query = "select * from {}".format(compartment)
        compartment_query = "select * from {} WHERE {} IN {};".format(compartment,"ImageNumber",list_str2)
#         compartment_query = "select * from {} WHERE {} IN {};".format(compartment,"Metadata_Well",list_str)
        plateDf_list.append(pd.read_sql(sql=compartment_query, con=conn))

    plateDf = reduce(lambda left,right: pd.merge(left,right,on=["TableNumber", "ImageNumber", "ObjectNumber"]), plateDf_list)
    end = time.time()
    print('time elapsed:',(end - start2)/60)
    
#     print(plateDf.columns[plateDf.columns.str.contains("Metadata_")])
    plateDfwMeta = pd.merge(plateDf, plateImageDf, on=["TableNumber", "ImageNumber"])
    
    del plateDf
#     gc.collect()
    
    plateDfwMeta = plateDfwMeta.loc[:,~plateDfwMeta.columns.duplicated()]
#     print(plateDfwMeta.shape)
#     print(plateDfwMeta.Image_Metadata_ImageSizeX.values[0])
#     plateDfwMeta=edgeCellFilter2(plateDfwMeta);  
    
    return plateDfwMeta   




def readSingleCellData_sqlalch_wellAndObject_subset(fileName,wells,n_rand_objs):
    import pandas as pd
    from sqlalchemy import create_engine
    from functools import reduce
    # sql_file="sqlite:////"+rootDir+'backend/JUANCHI_PILOT_1/JUANCHI_PILOT_1_PLATE_1/JUANCHI_PILOT_1_PLATE_1.sqlite'
    # sql_file="sqlite:////"+rootDirDrug+"/backend/"+batchName+"/"+p+"/"+p+".sqlite"    

#     wells=['A01', 'A02', 'A03', 'A04', 'A05', 'A06', 'A07', 'A08', 'A09',\
#        'H07', 'H08', 'H09', 'H10', 'H11', 'H12'];
    
    sql_file="sqlite:////"+fileName
    engine = create_engine(sql_file)
    conn = engine.connect()
    compartments=["cells", "cytoplasm", "nuclei"]
    # compartments=["Neurites","CellBodies","CellBodiesPlusNeurites","Nuclei","Cytoplasm"]

#     rand_img_num=np.random.choice(range(1,4000), n_rand_ims)
#     rand_img_num=np.array(range(50))
#     list_str="("
#     for i in rand_img_num:
#         list_str=list_str+str(i)+',' 
#     list_str=list_str[:-1]+")"
    
#     rand_img_num=wells[:40]
    rand_img_num=wells.copy()
    list_str="('"
    for i in rand_img_num:
        list_str=list_str+str(i)+"','" 
    list_str=list_str[:-2]+")"

#     compartment_query = "select * from {}".format("Image")
#     compartment_query = "select * from {} WHERE {} IN {};".format("Image","ImageNumber",list_str)
#     compartment_query = "select * from {} WHERE {} IN {};".\
    compartment_query = "select * from {} WHERE {} IN {};".\
    format("Image","Metadata_Well",list_str)
#     format("Image","Image_Metadata_Well",list_str)
    start1 = time.time()
    plateImageDf= pd.read_sql(sql=compartment_query, con=conn);
#     print(plateImageDf.columns[plateImageDf.columns.str.contains("Metadata_")])
#     print(plateImageDf['Metadata_Well'].unique())
    end1 = time.time()
    print('time elapsed:',(end1 - start1)/60)
    img_nums=plateImageDf.ImageNumber.unique().tolist()
    print(plateImageDf.shape,img_nums)
    
    
    list_str2="("
    for i in img_nums:
        list_str2=list_str2+str(i)+',' 
    list_str2=list_str2[:-1]+")"
    start2 = time.time()


    n_rand_objs=100
    rand_obj_num=np.random.choice(range(1,300), n_rand_objs)
    list_obj="("
    for i in rand_obj_num:
        list_obj=list_obj+str(i)+',' 
    list_obj=list_obj[:-1]+")"

    start2 = time.time()
    plateDf_list=[]
    for compartment in compartments:
    #         compartment_query = "select * from {}".format(compartment)
        compartment_query = "select * from {} WHERE {} IN {} AND {} IN {};".\
        format(compartment,"ImageNumber",list_str2,"ObjectNumber",list_obj)
    #         compartment_query = "select * from {} WHERE {} IN {};".format(compartment,"Metadata_Well",list_str)
        plateDf_list.append(pd.read_sql(sql=compartment_query, con=conn))    
        
    
    plateDf = reduce(lambda left,right: pd.merge(left,right,on=["TableNumber", "ImageNumber", "ObjectNumber"]), plateDf_list)
    end = time.time()
    print('time elapsed:',(end - start2)/60)
    
#     print(plateDf.columns[plateDf.columns.str.contains("Metadata_")])
    plateDfwMeta = pd.merge(plateDf, plateImageDf, on=["TableNumber", "ImageNumber"])
    
    del plateDf
#     gc.collect()
    
    plateDfwMeta = plateDfwMeta.loc[:,~plateDfwMeta.columns.duplicated()]
#     print(plateDfwMeta.shape)
#     print(plateDfwMeta.Image_Metadata_ImageSizeX.values[0])
#     plateDfwMeta=edgeCellFilter2(plateDfwMeta);  
    
    return plateDfwMeta   


def readSingleCellData_sqlalch_features_subset(fileName,selected_feature):

    start1 = time.time()
    # selected_feature='Cells_RadialDistribution_MeanFrac_mito_tubeness_16of16'
#     selected_feature='Cells_Intensity_IntegratedIntensity_DNA'
    # f2='Cells_Intensity_IntegratedIntensity_DNA'
    sql_file="sqlite:////"+fileName
    engine = create_engine(sql_file)
    conn = engine.connect()
    compartments=selected_feature.split("_")[0]
    query_cols = "TableNumber, ImageNumber, "+selected_feature#+", "+f2
    compartment_query = "select {} from {}".format(query_cols,compartments)
    plateDf=pd.read_sql(sql=compartment_query, con=conn)
    # del plateDf_list
    #     plateDf = plateDf.loc[:,["TableNumber", "ImageNumber",selected_feature]]
    plateDf=plateDf.dropna()

    img_query = "select * from {}".format("Image")
    plateImageDf= pd.read_sql(sql=img_query, con=conn);

    # plateImageDf2=pd.merge(plateMap2, plateImageDf, on=["Metadata_Plate", "Metadata_Well"])
    #     dmso_wells=plateMap2[plateMap2['Metadata_broad_sample']=='DMSO'];

    plateDfwMeta = pd.merge(plateDf, plateImageDf, on=["TableNumber", "ImageNumber"])

    end1 = time.time()
    print('time elapsed:',(end1 - start1)/60, " mins")
    return plateDfwMeta

def readSingleCellData_sqlalch_FeatureAndWell_subset(fileName,selected_feature,wells):

    start1 = time.time()
    # selected_feature='Cells_RadialDistribution_MeanFrac_mito_tubeness_16of16'
#     selected_feature='Cells_Intensity_IntegratedIntensity_DNA'
    # f2='Cells_Intensity_IntegratedIntensity_DNA'

    sql_file="sqlite:////"+fileName
    engine = create_engine(sql_file)
    conn = engine.connect()

    ######## Query wells from Image table
    ls_wells=wells.copy()
    list_str="('"
    for i in ls_wells:
        list_str=list_str+str(i)+"','" 
    list_str=list_str[:-2]+")"

    img_query = "select * from {} WHERE {} IN {};".\
    format("Image","Image_Metadata_Well",list_str)

    plateImageDf= pd.read_sql(sql=img_query, con=conn);
    img_nums=plateImageDf.ImageNumber.unique().tolist()

    list_str2="("
    for i in img_nums:
        list_str2=list_str2+str(i)+',' 
    list_str2=list_str2[:-1]+")"
    ###########################

    compartments=selected_feature.split("_")[0]
    query_cols = "TableNumber, ImageNumber, "+selected_feature#+", "+f2
    compartment_query = "select {} from {} WHERE {} IN {};".format(query_cols,compartments,"ImageNumber",list_str2)
    plateDf=pd.read_sql(sql=compartment_query, con=conn)

    plateDfwMeta = pd.merge(plateDf, plateImageDf, on=["TableNumber", "ImageNumber"])

    end1 = time.time()
    print('time elapsed:',(end1 - start1)/60, " mins")
    
    return plateDfwMeta
    
#### read samples of different replicates of the drug rep data 

# def read_single_cell_samples_drugX(Metadata_pert_id_dose,n_sample_per_replicate):
#     """ 
#     This function take drug+dose name and randomly sample single cells from all the replicates (sql files) of the input
#     and output single cell profiles into a single cell dataframe
  
#     Inputs: 
#     ++ Metadata_pert_id_dose   (dtype: str) --> str(Metadata_pert_id)+'_'+str(dose)
#     ++ n_sample_per_replicate   (dtype: int) --> number of random samples taken from each replicate
    
    
#     Returns: 
#     df (pandas dataframe): 
  
#     """ 
#     meta_lincs_repLevel=pd.read_csv("/home/ubuntu/bucket/projects/2018_04_20_Rosetta/\
#     workspace/results/synth_meta/meta_lincs_repLevel.csv")
    
    
#     readSingleCellData_sqlalch_well_subset(fileName,wells):
    
    
    
    
    
    
    
