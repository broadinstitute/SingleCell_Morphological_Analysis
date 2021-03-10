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
from functools import reduce
import skimage

def clusteringHists(DirsDict,wtANDmtDf_scaled,contLabel,d,nClus,feats2use,compartments,boxSize):
    rootDir=DirsDict['root']
    resultsDir=DirsDict['resDir']
#     DirsDict['imDir']=rootDir+"Mito_Morphology_input/images/"
    
    
#     d1=d.split(" ")[0]
    saveFormat='.png';#'.png'
#     plt.ioff()
    fig, axes = plt.subplots(1,2)
    # wtANDmtDf['clusterLabels']
    data2plotMut=wtANDmtDf_scaled[(wtANDmtDf_scaled['label'] == d)]['clusterLabels'].values
    data2plotWT=wtANDmtDf_scaled[wtANDmtDf_scaled['label'] == contLabel]['clusterLabels'].values

    histMut, bin_edges = np.histogram(data2plotMut,range(nClus+1), density=True)
    histWT, bin_edges = np.histogram(data2plotWT,range(nClus+1), density=True)

    histDiff=histMut-histWT;
    sortedDiff=np.sort(histDiff)
    # ind=[np.where(histDiff[i]==sortedDiff)[0][0] for i in range(len(histDiff))]
    ind=[]
    for i in range(len(histDiff)):
        iinndd = np.where(histDiff[i]==sortedDiff)[0].tolist()
        for j in range(len(iinndd)):
            if iinndd[j] not in ind:
                ind=ind+[iinndd[j]]
                break

    wtANDmtDf_scaled['clusterLabels2']=wtANDmtDf_scaled['clusterLabels'].replace(range(nClus),ind)
    data2plotMut=wtANDmtDf_scaled[~(wtANDmtDf_scaled['label'] == contLabel)]['clusterLabels2'].values
    data2plotWT=wtANDmtDf_scaled[wtANDmtDf_scaled['label'] == contLabel]['clusterLabels2'].values

    histMut, bin_edges = np.histogram(data2plotMut,range(nClus+1), density=True)
    histWT, bin_edges = np.histogram(data2plotWT,range(nClus+1), density=True)
    # def mapToNewLabelCats(histMut,histWT):


    sns.distplot(data2plotMut,kde=False,norm_hist=True,bins=bin_edges,label=d,ax=axes[0],color="r",hist_kws=dict(edgecolor="k"));
    sns.distplot(data2plotWT,kde=False,norm_hist=True,bins=bin_edges,label=contLabel,ax=axes[0],hist_kws=dict(edgecolor="k"))
    sns.distplot(data2plotMut,kde=False,hist=True,norm_hist=False,bins=bin_edges,label=d,ax=axes[1],color="r",hist_kws=dict(edgecolor="k"));
    sns.distplot(data2plotWT,kde=False,hist=True,norm_hist=False,bins=bin_edges,label=contLabel,ax=axes[1],hist_kws=dict(edgecolor="k"))

#   axes.xaxis.set_ticklabels(range(0,20,2)); 
    axes[0].set_ylabel('Density');axes[0].set_xlabel('cell category index');
    axes[1].set_ylabel('Histogram');axes[1].set_xlabel('cell category index');
    axes[0].legend();axes[1].legend();
    plt.tight_layout()
    os.system("mkdir -p "+resultsDir);
    fig.savefig(resultsDir+'/clusterDensity'+saveFormat)  

    meanWT=wtANDmtDf_scaled.loc[wtANDmtDf_scaled['label'] == contLabel,feats2use].mean()

    for c in range(len(histMut)):
        if histMut[c] > 0.001 or histWT[c] > 0.001:
    #         c=3;
            clusterDF=wtANDmtDf_scaled[wtANDmtDf_scaled['clusterLabels2'] == c].reset_index(drop=True)
            meanMutCluster=clusterDF.loc[~(clusterDF['label'] == contLabel),feats2use].mean();        
            diffOfMutMeanAndWTMean=pd.DataFrame(data=meanMutCluster.values-meanWT.values,columns=['Diff'],index=meanMutCluster.index);
            diffOfMutMeanAndWTMean.loc[:,'Diff2']=diffOfMutMeanAndWTMean.loc[:,'Diff'].abs()
            absFeatureImportanceSS=diffOfMutMeanAndWTMean.sort_values('Diff2',ascending=False)[:10];
            fig, axes = plt.subplots()
            sns.barplot(x='Diff', y=absFeatureImportanceSS.index, data=absFeatureImportanceSS,ax=axes)
            sns.despine()
            plt.tight_layout()   
            fig.savefig(resultsDir+'/cluster'+str(c)+'_barImpFeatures'+saveFormat)  
            plt.close('all')
            nSampleSCs=6
            if clusterDF.shape[0]> nSampleSCs:
                samples2plot=clusterDF.sort_values('dist2Mean',ascending=True).sample(nSampleSCs).reset_index(drop=True)
                title_str="Cluster "+str(c)
                f=visualize_n_SingleCell(compartments,samples2plot,boxSize,title=title_str)
                f.savefig(resultsDir+'/cluster'+str(c)+'_examplar'+saveFormat)     
                plt.close('all')    
                
    return




def visualize_n_SingleCell(channels,dfWithWTlabels,boxSize,title=""):
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
       
    ++ boxSize (int): Height or Width of the square bounding box
    
    ++ title
    
    Returns: 
    f (object): handle to the figure
  
    """
    
    halfBoxSize=int(boxSize/2);
#     print(channels)
    
    import skimage.io
    f, axarr = plt.subplots(dfWithWTlabels.shape[0], len(channels),figsize=(len(channels)*2,dfWithWTlabels.shape[0]*2));
    if len(title)>0:
        print(title)
        f.suptitle(title);
    
    f.subplots_adjust(hspace=0, wspace=0)


#     maxRanges={"DNA":8000,"RNA":6000,"Mito":6000,"ER":8000,"AGP":6000}
#     print(dfWithWTlabels.shape[0])
    for index in range(dfWithWTlabels.shape[0]):
               
        compressed=True
        if compressed:
            xCenter0=int(dfWithWTlabels.loc[index,'Nuclei_Location_Center_X'])
            yCenter0=int(dfWithWTlabels.loc[index,'Nuclei_Location_Center_Y'])            
            compressed_im_size=1080;
            original_im_size=2160;
            compRatio=(compressed_im_size/original_im_size);
            xCenter,yCenter=int(xCenter0*compRatio),\
            int(yCenter0*compRatio)  
        else:
            xCenter=int(dfWithWTlabels.loc[index,'Nuclei_Location_Center_X'])
            yCenter=int(dfWithWTlabels.loc[index,'Nuclei_Location_Center_Y'])            
        
#         print(xCenter,yCenter)
        
        cpi=0;
        for c in channels:
            if c=='Outline':
                imPath=dfWithWTlabels.loc[index,'URL_CellOutlines'];
            else:
#                 ch_D=dfWithWTlabels.loc[index,'Image_FileName_Orig'+c];
                ch_D=dfWithWTlabels.loc[index,'FileName_Orig'+c];
#                 print(ch_D)
    #         imageDir=imDir+subjectID+' Mito_Morphology/'
#                 imageDir=dfWithWTlabels.loc[index,'Image_PathName_Orig'+c]+'/'
                imageDir=dfWithWTlabels.loc[index,'PathName_Orig'+c]+'/'
                imPath=imageDir+ch_D
        
            
            imD=skimage.io.imread(imPath)[yCenter-halfBoxSize:yCenter+halfBoxSize,xCenter-halfBoxSize:xCenter+halfBoxSize]
#             axarr[index,cpi].imshow(imD,cmap='gray',clim=(0, maxRanges[c]));axarr[0,cpi].set_title(c);
            axarr[index,cpi].imshow(imD,cmap='gray');axarr[0,cpi].set_title(c);
            cpi+=1        

#         Well=dfWithWTlabels.loc[index,'Metadata_Well']
#         Site=str(dfWithWTlabels.loc[index,'Metadata_Site'])
#         imylabel=Well+'\n'+Site
#         axarr[index,0].set_ylabel(imylabel);            
    #         subjectID=dfWithWTlabels.loc[index,'subject']
    
        if "label" in dfWithWTlabels.columns.tolist():
            imylabel=dfWithWTlabels.label[index]#+'\n'+subjectID
            imylabel2="-".join(imylabel.split('-')[0:2])
            axarr[index,0].set_ylabel(imylabel2);
    # #     plt.tight_layout() 

    for i in range(len(channels)):
        for j in range(dfWithWTlabels.shape[0]):
            axarr[j,i].xaxis.set_major_locator(plt.NullLocator())
            axarr[j,i].yaxis.set_major_locator(plt.NullLocator())
            axarr[j,i].set_aspect('auto')
    
    return f




def CP_to_RGB_single(im_cp):
    # change channels first to channels last format
    channel_first=False
    if im_cp.shape[0]<10:
        channel_first=True
        im_cp = np.moveaxis(im_cp, 0, 2)
    col1 = np.array([0, 0, 255], dtype=np.uint8)
    col2 = np.array([0, 255, 0], dtype=np.uint8)
    col3 = np.array([255, 255, 0], dtype=np.uint8)
    col4 = np.array([255, 150, 0], dtype=np.uint8)
    col5 = np.array([255, 0, 0], dtype=np.uint8)
    channel_colors=[col1,col2,col3,col4,col5]
    comb_pars=[3,2,3,2,2]
    colorImagesList=[]
#     print(im_cp.shape[2])
    for i in range(im_cp.shape[2]):
        image_gray=im_cp[:,:,i]
        image_gray_normalized,_=normalize(image_gray)
        image_color=colorize_image(image_gray_normalized, channel_colors[i])
        colorImagesList.append(image_color)
        colorImagesList2 = [a * b.astype(np.uint16) for a, b in zip(comb_pars, colorImagesList)]
    colorImage0,_=normalize(sum(colorImagesList2));
    colorImage0=skimage.img_as_float64(colorImage0)
#         print(image_gray.shape,image_gray_normalized.shape,image_color.shape,colorImage0.shape)
    if channel_first:
        colorImage = np.moveaxis(colorImage0, 2, 0)
    else:
        colorImage=colorImage0.copy()
    return colorImage

def colorize_image(img, col):

    # rescale image
    img_float = img.astype(np.float)
    img_float = img_float / 255

    # colorize
    img_col_float = np.reshape(img_float, img_float.shape + (1,)) * col
    img_col_byte = img_col_float.astype(np.uint8)

    return img_col_byte
#         [64, 5, 128, 128]
#         return im_RGB

def normalize(img):

    # normalize to [0,1]
    img=abs(img.min())+img
    percentile = 99.95
    high = np.percentile(img, percentile)
    low = np.percentile(img, 100-percentile)

    img = np.minimum(high, img)
    img = np.maximum(low, img)

#     img = (img - low) / (high - low) # gives float64, thus cast to 8 bit later
#     vmin, vmax = scipy.stats.scoreatpercentile(image, (0.05, 99.95))
#     vmax = min(vmax, pmax)
    image_01 = skimage.exposure.rescale_intensity(img, in_range=(low, high))
    
#     image = skimage.exposure.rescale_intensity(img, in_range=(-1, 1))
    image_01[image_01>1]=1
    image_01[image_01<0]=0
#     image[image<-1]=-1
#     print(image.min(),image.max())
    img_255 = skimage.img_as_ubyte(image_01)
#     print(img_255.min(),img_255.max())
#     print(image_01.min(),image_01.max())
    return img_255, image_01  