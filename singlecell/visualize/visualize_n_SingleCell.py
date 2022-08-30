"""
This script contains usefull functions used in the notebooks

@author: mhaghigh
"""
import pandas as pd
import numpy as np
import seaborn as sns
# from sklearn import preprocessing
# import pickle
import matplotlib.pyplot as plt
# from imblearn.over_sampling import SMOTE,RandomOverSampler
import os
from functools import reduce
from sklearn.cluster import KMeans
# import skimage
from skimage import exposure
from skimage.transform import resize
import skimage.io
# from utils.read_data import *
# from read_data import *


# import skimage.io




def visualize_n_SingleCell_pooled(channels,sc_df,boxSize,im_size,title=""):
    """ 
    This function plots the single cells correspoding to the input single cell dataframe
  
    Inputs: 
    ++ sc_df   (pandas df) size --> (number of single cells)x(columns): 
    input dataframe contains single cells profiles as rows (make sure it has "Nuclei_Location_Center_X"or"Y" columns)
    
    ++ channels (dtype: list): list of channels to be displayed as columns of output image
           example: channels=['Mito','AGP','Brightfield','ER','DNA','Outline']
        * If Outline exist in the list of channels; function reads the outline image address from 
          "URL_CellOutlines" column of input dataframe, therefore, check that the addresses are correct
           before inputing them to the function, and if not, modify before input!
       
    ++ boxSize (int): Height or Width of the square bounding box
    
    ++ title
    
    Returns: 
    f (object): handle to the figure
  
    """

    halfBoxSize=int(boxSize/2);
#     print(channels)
    
    
    f, axarr = plt.subplots(sc_df.shape[0], len(channels),figsize=(len(channels)*2,sc_df.shape[0]*2));
    if len(title)>0:
        print(title)
        f.suptitle(title);
    
    f.subplots_adjust(hspace=0, wspace=0)

    align_column_ch_name_map={'DNA':'DAPI_Painting','ER':'ConA','Mito':'Mito','Phalloidin':'Phalloidin','WGA':'WGA'}
    
    
    for index in range(sc_df.shape[0]):
               
        xCenter=int(sc_df.loc[index,'Nuclei_Location_Center_X'])
        yCenter=int(sc_df.loc[index,'Nuclei_Location_Center_Y'])   
            
        
        cpi=0;
        for c in channels:
            if c=='Outline':
                imPath=sc_df.loc[index,'Path_Outlines'];
#                 im_size=sc_df["Width_CorrDNA"].values[0]   #cp220
#                 im_size= 5500
#                 print(imPath)
                
                if os.path.exists(imPath):
                    imD2 = resize(skimage.io.imread(imPath),(im_size,im_size),\
                                 mode='constant',preserve_range=True,order=0).\
                astype('uint8')[yCenter-halfBoxSize:yCenter+halfBoxSize,xCenter-halfBoxSize:xCenter+halfBoxSize]
                    
                else:
                    imD2=np.zeros((boxSize,boxSize))
                clim_max=imD2.max()
            else:
#                 ch_D=sc_df.loc[index,'Image_FileName_Orig'+c];
                ch_D=sc_df.loc[index,'FileName_Corr'+c];
#                 print(ch_D)
    #         imageDir=imDir+subjectID+' Mito_Morphology/'
#                 imageDir=sc_df.loc[index,'Image_PathName_Orig'+c]+'/'
                imageDir=sc_df.loc[index,'PathName_Corr'+c]+'/'
                imPath=imageDir+ch_D
#                 print(imPath)
#                 print(yCenter,xCenter)

                xCenterC=int(sc_df.loc[index,'Nuclei_Location_Center_X'])+\
    int(sc_df.loc[index,'Align_Xshift_'+align_column_ch_name_map[c]])
                yCenterC=int(sc_df.loc[index,'Nuclei_Location_Center_Y'])+\
    int(sc_df.loc[index,'Align_Yshift_'+align_column_ch_name_map[c]]) 
            
#                 print(yCenterC,xCenterC)
#                 print(imPath)
#                 print('im_size',np.squeeze(skimage.io.imread(imPath)).shape)
                imD=np.squeeze(skimage.io.imread(imPath))
                imD1=exposure.rescale_intensity(imD,in_range=(imD.min(),np.percentile(imD, 99.95)))#output is  unit16
#                 and scaled to range 0,65535
                clim_max=imD1.max()
        
                imD2=imD1[yCenterC-halfBoxSize:yCenterC+halfBoxSize,xCenterC-halfBoxSize:xCenterC+halfBoxSize]
#                 print(imD1.min(),imD1.max())
#                 print(np.squeeze(skimage.io.imread(imPath)).shape)
#             axarr[index,cpi].imshow(imD,cmap='gray',clim=(0, maxRanges[c]));axarr[0,cpi].set_title(c);
#             print(imD.shape,'h')
            axarr[index,cpi].imshow(imD2,cmap='gray',clim=(0, clim_max));axarr[0,cpi].set_title(c);
            axarr[index,cpi].axes.xaxis.set_ticks([])
            axarr[index,cpi].axes.yaxis.set_ticks([])
            cpi+=1        

#         Well=sc_df.loc[index,'Metadata_Well']
#         Site=str(sc_df.loc[index,'Metadata_Site'])
#         imylabel=Well+'\n'+Site

        if 'label' in sc_df.columns:
            imylabel=sc_df.loc[index,'label']+'\n'+sc_df.loc[index,'Metadata_Foci_Barcode_MatchedTo_Barcode'][0:9]
        else:
            imylabel=sc_df.loc[index,'Metadata_Foci_Barcode_MatchedTo_Barcode'][0:12]
#         print(imylabel)
        axarr[index,0].set_ylabel(imylabel);            
            

    return f


def form_image_y_label_string(sc_df,info_columns):

    image_y_label='\n'.join([str(sc_df[ic]) for ic in info_columns])

    return image_y_label




def visualize_n_SingleCell(channels , sc_df , boxSize , outline = False, color=False,\
                           info_columns=[], title="",compressed=False,compressed_im_size=None):
    """ 
    This function plots the single cells correspoding to the input single cell dataframe
  
    Inputs: 
    ++ sc_df   (pandas df) size --> (number of single cells)x(columns): 
    input dataframe contains single cells profiles as rows (make sure it has "Nuclei_Location_Center_X"or"Y" columns)
    
    ++ channels (dtype: list): list of channels to be displayed as columns of output image
           example: channels=['Mito','AGP','Brightfield','ER','DNA','Outline']
        * If Outline exist in the list of channels; function reads the outline image address from 
          "URL_CellOutlines" column of input dataframe, therefore, check that the addresses are correct
           before inputing them to the function, and if not, modify before input!
       
    ++ boxSize (int): Height or Width of the square bounding box
    ++ info_columns (list): are the columns containing info we want to add next to each single cell
    Optional Inputs:
    ++ title (str)
    ++ compressed (bool) default is False,if set to True the next parameter is not optional anymore and should be provided
    ++ compressed_im_size (int), for eaxample for lincs compressed is 1080
    
    Returns: 
    f (object): handle to the figure
  
    """
    
    if compressed:
        
        original_im_size=sc_df['Image_Width_OrigDNA'].values[0]
        #         compressed_im_size=1080;
        compRatio=(compressed_im_size/original_im_size);
        
        sc_df['Nuclei_Location_Center_X']=sc_df['Nuclei_Location_Center_X']*compRatio
        sc_df['Nuclei_Location_Center_Y']=sc_df['Nuclei_Location_Center_Y']*compRatio          

    
    halfBoxSize=int(boxSize/2);
#     print(channels)
    
    columns_count=len(channels)+int(outline)+int(color)
    rows_count=sc_df.shape[0]
    
    
    
    f, axarr = plt.subplots(rows_count, columns_count,figsize=(columns_count*2,rows_count*2));
    if len(title)>0:
#         print(title)
        f.suptitle(title);
    
    f.subplots_adjust(hspace=0, wspace=0)
    for index in range(rows_count):
               
        xCenter=int(sc_df.loc[index,'Nuclei_Location_Center_X'])
        yCenter=int(sc_df.loc[index,'Nuclei_Location_Center_Y'])            
        
        sc_collage_row=np.zeros((boxSize,boxSize,columns_count))
        
        for ci in range(len(channels)):
            if index==0:
                axarr[index,ci].set_title(channels[ci]);
            
            ch_fName=sc_df.loc[index,'FileName_Orig'+channels[ci]];
            ch_pName=sc_df.loc[index,'PathName_Orig'+channels[ci]];
            
            image=skimage.io.imread(ch_pName+'/'+ch_fName)
            image_cropped = crop_single_cell_image(image,xCenter,yCenter,halfBoxSize)
            
           
            if 1:
                image_cropped= exposure.rescale_intensity(image_cropped,in_range=(image.min(),np.percentile(image, 99.95)))
            
            sc_collage_row[:,:,ci]=image_cropped
           
        
        if outline:
            imPath=sc_df.loc[index,'Path_Outlines'];
            imD1=skimage.io.imread(imPath)
            
            if compressed:
                imD1=skimage.transform.resize(imD1,[compressed_im_size,compressed_im_size],\
                                              mode='constant',preserve_range=True,order=0) 
                
                sc_collage_row[:,:,ci+1] = crop_single_cell_image(imD1,xCenter,yCenter,halfBoxSize)
                
                
        if color: 
            print('not available here yet')
     
        for c in range(columns_count):    
            axarr[index,c].imshow(sc_collage_row[:,:,c],cmap='gray',clim=(0, sc_collage_row[:,:,c].max()));
            axarr[index,c].axes.xaxis.set_ticks([]);axarr[index,c].axes.yaxis.set_ticks([])
        
        ## add info in the columns specified in the input next to to each single cell
        if info_columns:
            image_y_label = form_image_y_label_string(sc_df.loc[index,info_columns],info_columns);
            axarr[index,0].set_ylabel(image_y_label);
            
    return 


def crop_single_cell_image(image, xCenter,yCenter,halfBoxSize):


    image_cropped=image[yCenter-halfBoxSize:yCenter+halfBoxSize,\
                        xCenter-halfBoxSize:xCenter+halfBoxSize]
    return image_cropped
    
    
    
    
# def cp_image_to_rgb():
    
#         col1 = np.array([0, 0, 255], dtype=np.float) #blue
#     col2 = np.array([255, 0, 0], dtype=np.float) # red
#     col3 = np.array([0, 255, 0], dtype=np.float) #lime
#     col4 = np.array([255,255,0], dtype=np.float) # yellow
#     col5 = np.array([255,0,255], dtype=np.float)#magneta
    
#     channel_colors=[col1,col2,col3,col4,col5]
# #     comb_pars=[3,2,3,2,2]
#     comb_pars=[.1,.1,.1,.1,.1]
#     colorImagesList=[]
#     for i in range(im_cp.shape[2]):
#         image_gray=im_cp[:,:,i]
# #         print(skimage.color.gray2rgb(image_gray).shape)
# #         image_gray_normalized=normalize(image_gray)
#         image_color=(skimage.color.gray2rgb(image_gray).astype(float)/255) *channel_colors[i]
# #         print('max',image_color.max(),image_color.shape)
# #         image_color=colorize_image(image_gray_normalized, channel_colors[i])
#         colorImagesList.append(image_color)
# #         colorImagesList2 = [a * b.astype(np.uint16) for a, b in zip(comb_pars, colorImagesList)]

#     colorImage0=sum(colorImagesList).astype(np.uint8);