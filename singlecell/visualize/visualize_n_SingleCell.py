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
from skimage.segmentation import flood_fill
# from utils.read_data import *
# from read_data import *


# import skimage.io



# def visualize_n_SingleCell(channels , sc_df , boxSize , outline = False, color=False,\
#                            info_columns=[], title="",compressed=False,compressed_im_size=None):
    
    
def visualize_n_SingleCell_pooled(channels,sc_df, boxSize, im_size, outline = False, color=False,\
                           info_columns=[], title=""):
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
    
    plt.ioff()
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
                imD2 = crop_single_cell_image(imD1,xCenterC,yCenterC,halfBoxSize)
        
#                 imD2=imD1[yCenterC-halfBoxSize:yCenterC+halfBoxSize,xCenterC-halfBoxSize:xCenterC+halfBoxSize]
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
            
    plt.ion()
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
        metadata_cols_4size = sc_df.columns[sc_df.columns.str.contains('Width|ImageSize')]
        if len(metadata_cols_4size):
            original_im_size=sc_df[metadata_cols_4size[0]].values[0]
        else:
            raise Exception("No metadata columns for inferring image size are detected! Please enter original_im_size here!")
            
#         original_im_size=sc_df['Image_Width_OrigDNA'].values[0]
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
            
#             print(ch_pName+'/'+ch_fName)
            image=skimage.io.imread(ch_pName+'/'+ch_fName)
            image_cropped = crop_single_cell_image(image,xCenter,yCenter,halfBoxSize)
#             print(image_cropped.shape,xCenter,yCenter)
           
            if 1:
                image_cropped= exposure.rescale_intensity(image_cropped,in_range=(image.min(),np.percentile(image, 99.95)))
            
            sc_collage_row[:,:,ci]=image_cropped
           
                
                
        for c in range(len(channels)):    
            axarr[index,c].imshow(sc_collage_row[:,:,c],interpolation=None,cmap='gray',clim=(0, sc_collage_row[:,:,c].max()));
#             axarr[index,c].axes.xaxis.set_ticks([]);axarr[index,c].axes.yaxis.set_ticks([])                
                
        if color: 
#             print('not available here yet')
            color_im = CP_to_RGB_single(sc_collage_row[:,:,:len(channels)],channels)
#             print(color_im.max(),color_im.min())
            axarr[index,c+1].imshow(color_im);
#             axarr[index,c+1].axes.xaxis.set_ticks([]);axarr[index,c+1].axes.yaxis.set_ticks([]);
            
            if index==0:
                axarr[index,c+1].set_title('composite'); 
            
        if outline:
            imPath=sc_df.loc[index,'Path_Outlines'];
            imD1=skimage.io.imread(imPath)
            
            if compressed:
                imD1=skimage.transform.resize(imD1,[compressed_im_size,compressed_im_size],\
                                              mode='constant',preserve_range=True,order=0) 
                
            outline_im = crop_single_cell_image(imD1,xCenter,yCenter,halfBoxSize)
            axarr[index,columns_count-1].imshow(outline_im,cmap='gray');
#             axarr[index,columns_count-1].axes.xaxis.set_ticks([]);axarr[index,columns_count-1].axes.yaxis.set_ticks([]);

        for c in range(columns_count):  
            axarr[index,c].axes.xaxis.set_ticks([]);axarr[index,c].axes.yaxis.set_ticks([])    
        ## add info in the columns specified in the input next to to each single cell
        if info_columns:
            image_y_label = form_image_y_label_string(sc_df.loc[index,info_columns],info_columns);
            axarr[index,0].set_ylabel(image_y_label);
            
    return f


def crop_single_cell_image(image, xCenter,yCenter,halfBoxSize):

    im_h,im_w=image.shape;
#     print(im_w,im_h)
    before_y_pad=0
    after_y_pad=0
    before_x_pad=0
    after_x_pad=0

    if xCenter-halfBoxSize<0:
        before_x_pad=abs(xCenter-halfBoxSize)

    if yCenter-halfBoxSize<0:
        before_y_pad=abs(yCenter-halfBoxSize)

    if xCenter+halfBoxSize>im_w:
        after_x_pad=abs(im_w-xCenter-halfBoxSize)

    if yCenter+halfBoxSize>im_h:
        after_y_pad=abs(im_h-yCenter-halfBoxSize)

    image_cropped=image[np.maximum(yCenter-halfBoxSize,0):np.minimum(yCenter+halfBoxSize,im_h),\
                        np.maximum(xCenter-halfBoxSize,0):np.minimum(xCenter+halfBoxSize,im_w)]
#     print('image_cropped',image_cropped.shape)
    if np.max([before_y_pad, after_y_pad,before_x_pad, after_x_pad])>0:
        image_cropped=np.pad(image_cropped, ((before_y_pad, after_y_pad), (before_x_pad, after_x_pad)), 'minimum')
#         print('image_cropped',image_cropped.shape)

    return image_cropped
    
    

    
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
            image_cropped = crop_single_cell_image(image,xCenter,yCenter,halfBoxSize)
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
    
    

def CP_to_RGB_single(im_cp, channels):
    """ 
    This function takes a cell paiting image (as channels last array) and converts it to RGB 
  
    Inputs: 
    ++ im_cp   (np array) size --> (width)x(height)x(channels): 
    input dataframe contains single cells profiles as rows (make sure it has "Nuclei_Location_Center_X"or"Y" columns)
    
    ++ channels (dtype: list): list of channels to be displayed as columns of output image
           example: channels=['Mito','AGP','Brightfield','ER','DNA','Outline']
    
    Returns: 
    colorImage0 (np array): dims-> width , height , 3 (RGB channels)
  
    """    
    import matplotlib.colors as mcolors
#     print(im_cp.max
    # change channels first to channels last format
#     im_cp = np.moveaxis(im_cp, 0, 2)
#     col1 = np.array([0, 0, 255], dtype=np.float) #blue
#     col2 = np.array([255, 0, 0], dtype=np.float) # red
#     col3 = np.array([0, 255, 0], dtype=np.float) #lime
#     col4 = np.array([255,255,0], dtype=np.float) # yellow
#     col5 = np.array([255,0,255], dtype=np.float)#magneta
#     channel_colors=[col1, col2 ,col3 ,col4 ,col5]
    
    depth=255
    # channels_colormap= {'DAPI':'Blue', 'ER'=Green, 'RNA'=Yellow, 'AGP':Red (or orange) 'Mito' = Magenta (or red)
    channels_colormap= {'DNA':'Blue', 'ER':'Green', 'RNA':'Yellow', 'AGP':'Red', 'Mito':'Magenta',\
                       'DAPI':'Blue'}
    channel_colors = [np.array(mcolors.to_rgb(channels_colormap[c]))*depth for c in channels]
    

#     comb_pars=[3,2,3,2,2]
#     comb_pars=[.1,.1,.1,.1,.1]
    comb_pars = [1/im_cp.shape[2]]*im_cp.shape[2]
    colorImagesList=[]
    for i in range(im_cp.shape[2]):
        image_gray=im_cp[:,:,i]
        image_color=(skimage.color.gray2rgb(image_gray).astype(float)/depth) *channel_colors[i]*comb_pars[i]
#         print('max',image_color.max(),image_gray.max(),image_color.shape)
        colorImagesList.append(image_color)

    colorImage0=sum(colorImagesList).astype(np.uint8);
#     print(colorImage0.shape)
#     print(colorImage0.min(),colorImage0.max())
#     for j in range(3):
#         colorImage0[:,:,j]=normalize(colorImage0[:,:,j])

    colorImage0 = skimage.exposure.rescale_intensity(colorImage0,out_range=(0,depth)).astype(np.uint8)

    return colorImage0    
    
    
def normalize(img):
    # normalize to [0,1]
    percentile = 99.95
    high = np.percentile(img, percentile)
    low = np.percentile(img, 100-percentile)

#     img = np.minimum(high, img)
#     img = np.maximum(low, img)
    
#     high=np.max(img)
#     low=np.min(img)

#     img = (img - low) / (high - low) # gives float64, thus cast to 8 bit later
#     vmin, vmax = scipy.stats.scoreatpercentile(image, (0.05, 99.95))
#     vmax = min(vmax, pmax)
#     image = skimage.exposure.rescale_intensity(img, in_range=(low, high))
    image = skimage.exposure.rescale_intensity(img, in_range=(low, high),out_range=(0,255)).astype(np.uint8)
#     image = skimage.exposure.rescale_intensity(img, in_range=(-1, 1))
#     image[image>1]=1
#     image[image<-1]=-1
#     print(image.max(),image.min())
#     img = skimage.img_as_ubyte(image)
#     print(image.max(),image.min())
    return image    

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