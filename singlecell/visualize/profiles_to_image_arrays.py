"""
@author: mhaghigh
"""

import os
import pandas as pd
import numpy as np

from skimage import exposure
from skimage.transform import resize
import skimage.io
from skimage.segmentation import flood_fill
from . import crop_cell


def scCPs_to_scIMs(channels , sc_df , boxSize , center_by='Nuclei_Location_Center', mask_cells = False ,\
                   compressed=False,compressed_im_size=None):
    """ 
    This function takes a pandas dataframe of single cell profiles and forms/outputs an array of single cells
  
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
    ++ compressed (bool) default is False,if set to True the next parameter is not optional anymore and should be provided
    ++ compressed_im_size (int), for eaxample for lincs compressed is 1080
    
    Returns: 
    im_arr (np array): dims-> n SCs , height , width , n CHs 
  
    """
    
    sc_df = sc_df.reset_index(drop=True)
    
    if compressed:
        metadata_cols_4size = sc_df.columns[sc_df.columns.str.contains('Width|ImageSize')]
        if len(metadata_cols_4size):
            original_im_size=sc_df[metadata_cols_4size[0]].values[0]
        else:
            raise Exception("No metadata columns for inferring image size are detected! Please enter original_im_size here!")
            
#         original_im_size=sc_df['Image_Width_OrigDNA'].values[0]
        #         compressed_im_size=1080;
        compRatio=(compressed_im_size/original_im_size);
        
        sc_df[center_by + '_X']=sc_df[center_by + '_X']*compRatio
        sc_df[center_by + '_Y']=sc_df[center_by + '_Y']*compRatio          

    
    halfBoxSize=int(boxSize/2);
#     print(channels)
    
    rows_count=sc_df.shape[0]
    
    im_arr=np.zeros((rows_count,boxSize,boxSize,len(channels)))
    
    for index in range(rows_count):
               
        xCenter=int(sc_df.loc[index, center_by + '_X'])
        yCenter=int(sc_df.loc[index, center_by + '_Y'])            
        
        for ci in range(len(channels)):
       
            ch_fName=sc_df.loc[index,'FileName_Orig'+channels[ci]];
            ch_pName=sc_df.loc[index,'PathName_Orig'+channels[ci]];
            
#             print(ch_pName+'/'+ch_fName)
            image=skimage.io.imread(ch_pName+'/'+ch_fName)
            image_cropped = crop_cell.crop_single_cell_image(image,xCenter,yCenter,halfBoxSize)
#             print(image.shape,image_cropped.shape,xCenter,yCenter)
           
            if 0:
                image_cropped= exposure.rescale_intensity(image_cropped,in_range=(image.min(),np.percentile(image, 99.95)))
           
                
            im_arr[index,:,:,ci]=image_cropped
            
        if mask_cells:
            imPath=sc_df.loc[index,'Path_Outlines'];
#             print(imPath,yCenter,xCenter)
            imD1=skimage.io.imread(imPath)

            cellMask1 = skimage.segmentation.flood_fill(imD1, (yCenter,xCenter), 1,connectivity=1)
            cellMask1[cellMask1 != 1] = 0  
            
            if compressed: 
                cellMask1 = skimage.transform.resize(cellMask1,\
                            [compressed_im_size,compressed_im_size],\
                            mode='constant',preserve_range=True,order=0).astype('uint8')
                cellMask1[cellMask1 != 1] = 0 
            

            cellMask2=cellMask1[yCenter-halfBoxSize:yCenter+halfBoxSize,\
                                xCenter-halfBoxSize:xCenter+halfBoxSize]                
#             print(len(channels),cellMask2.shape,cellMask1.shape) #1 (400, 400) (1040, 1388)
            
            cellMask4=np.repeat(cellMask2[:,:,np.newaxis],len(channels),axis=2)
#             print(cellMask4.shape,im_arr.shape,im_arr[index,:,:,:].shape) #(400, 400, 1) (100, 400, 400, 1)
            im_arr[index,:,:,:]=np.multiply(cellMask4,im_arr[index,:,:,:])
#             im_arr[index,:,:,:]=np.multiply(cellMask4,np.squeeze(im_arr[index,:,:,:], axis=0))           

    return im_arr