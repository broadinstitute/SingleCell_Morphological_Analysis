"""
@author: mhaghigh
"""
import pandas as pd
import numpy as np
import seaborn as sns
from sklearn import preprocessing
import pickle
import matplotlib.pyplot as plt
import os
import time

import pandas as pd
from sqlalchemy import create_engine
from functools import reduce
# from filter_edge_single_cells import edgeCellFilter
# import filter_edge_single_cells

from ..preprocess import filter_out_edge_single_cells
################################################################################



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



def readSingleCellData_sqlalch_well_subset(fileName,wells,meta_well_col_str):
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
#     meta_well_col_str="Metadata_Well"

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
#     print('time elapsed:',(end1 - start1)/60)
    img_nums=plateImageDf.ImageNumber.unique().tolist()
#     print(plateImageDf.shape,img_nums)
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
#     print('time elapsed:',(end - start2)/60)
    
#     print(plateDf.columns[plateDf.columns.str.contains("Metadata_")])
    plateDfwMeta = pd.merge(plateDf, plateImageDf, on=["TableNumber", "ImageNumber"])
    
    del plateDf
#     gc.collect()
    
    plateDfwMeta = plateDfwMeta.loc[:,~plateDfwMeta.columns.duplicated()]
#     print(plateDfwMeta.shape)
#     print(plateDfwMeta.Image_Metadata_ImageSizeX.values[0])
#     plateDfwMeta=filter_out_edge_single_cells.edgeCellFilter2(plateDfwMeta);  
    
    return plateDfwMeta   




def readSingleCellData_sqlalch_wellAndObject_subset(fileName,wells,meta_well_col_str,n_rand_objs):
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
    format("Image",meta_well_col_str,list_str)
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
    format("Image","Metadata_Well",list_str)

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
    
    
def readSingleCellData_sqlalch_FeaturesAndWells(fileName,selected_features,wells):

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
    format("Image","Metadata_Well",list_str)

    plateImageDf= pd.read_sql(sql=img_query, con=conn);
    img_nums=plateImageDf.ImageNumber.unique().tolist()

    list_str2="("
    for i in img_nums:
        list_str2=list_str2+str(i)+',' 
    list_str2=list_str2[:-1]+")"
    ###########################

    plateDf_list=[]
    for selected_feature in selected_features:
        compartments=selected_feature.split("_")[0]
        query_cols = "TableNumber, ImageNumber, ObjectNumber, "+selected_feature#+", "+f2
        compartment_query = "select {} from {} WHERE {} IN {};".format(query_cols,compartments,"ImageNumber",list_str2)
        plateDf_list.append(pd.read_sql(sql=compartment_query, con=conn))
        
    plateDf = reduce(lambda left,right: pd.merge(left,right,on=["TableNumber", "ImageNumber", "ObjectNumber"]), plateDf_list)        

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
    
    