#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Dec 2020

@author: mhaghigh
"""

'''
Currently running on this EC2 instance:
    Marzieh_?

Use the following command on your machine by executing:

    $ cd ~/workspace_rosetta/workspace/software/2018_04_20_Rosetta
    $ bokeh serve sample_field_mito_images.py --allow-websocket-origin=ec2-54-90-165-206.compute-1.amazonaws.com:8181 --port 8181
# bokeh serve sample_field_mito_images.py --allow-websocket-origin=ec2-54-164-178-213.compute-1.amazonaws.com:8181 --port 8181


Then navigate to the URL:

    http://ec2-54-90-165-206.compute-1.amazonaws.com:8181/sample_field_mito_images
    http://ec2-54-164-178-213.compute-1.amazonaws.com:8181/sample_field_mito_images

.. 

'''
print('hello')
try:
    from functools import lru_cache
except ImportError:
    # Python 2 does stdlib does not have lru_cache so let's just
    # create a dummy decorator to avoid crashing
    print ("WARNING: Cache for this example is available on Python 3 only.")
    def lru_cache():
        def dec(f):
            def _(*args, **kws):
                return f(*args, **kws)
            return _
        return dec

from os.path import dirname, join
from PIL import Image
import pandas as pd
import os
from bokeh.io import curdoc
from bokeh.layouts import row, column, layout, gridplot
from bokeh.models import ColumnDataSource,Button
from bokeh.models.widgets import PreText, Select,TextInput
from bokeh.plotting import figure
import skimage.io
import math
import numpy as np
from scipy import stats
import matplotlib as plt
import matplotlib.cm as cm
import copy
import sys
sys.path.insert(0, "../")
# import utils.preprocessing
from utils import visualize_data,read_data,analyze_data
import skimage

colormap =cm.get_cmap("OrRd") #choose any matplotlib colormap here
bokehpalette = [plt.colors.rgb2hex(m) for m in colormap(np.arange(0,colormap.N,50))]



############################################# Drug Rep data sets rootDir
rootDirDrug='/home/ubuntu/bucket/projects/2015_10_05_DrugRepurposing_AravindSubramanian_GolubLab_Broad/'
ImFilesRoot=rootDirDrug+'2016_04_01_a549_48hr_batch1_compressed/images'

drug_list_rank=pd.read_excel("/home/ubuntu/bucket/projects/2016_08_01_RadialMitochondriaDistribution_donna/\
workspace/Metadata_drugRep/drugList_20210115_uncorrForSiteAgg.xlsx",index_col=0).reset_index(drop=True)
meta_lincs=pd.read_csv("/home/ubuntu/bucket/projects/2018_04_20_Rosetta/\
workspace/results/synth_meta/meta_lincs_repLevel.csv")
meta_lincs.Metadata_mmoles_per_liter=meta_lincs.Metadata_mmoles_per_liter.values.round(2)
meta_lincs2=meta_lincs.groupby(['Metadata_broad_sample','Metadata_mmoles_per_liter','Metadata_Plate','Metadata_Well']).size().reset_index()


rootDirDrug='/home/ubuntu/bucket/projects/2015_10_05_DrugRepurposing_AravindSubramanian_GolubLab_Broad/workspace'
batchName='2016_04_01_a549_48hr_batch1_Mito_Project'



############################################# Mito Project rootDir
rootDirMito='/home/ubuntu/bucket/projects/2016_08_01_RadialMitochondriaDistribution_donna/workspace/'
# plateName='RC4_IF_15'
# SQL_file_path=rootDir+'/backend/PILOT_1_maxproj/'+plateName+'/'+plateName+'.sqlite'
# batchName=SQL_file_path.split('backend/')[1].split('/')[0]

############################# form ticker 1 
# dfRank4bokeh1=meta_lincs.groupby(['Metadata_broad_sample','Metadata_mmoles_per_liter']).size().reset_index()
dfRank4bokeh1=drug_list_rank.copy()
dfRank4bokeh1['ticker1']=dfRank4bokeh1['Metadata_broad_sample'].astype(str)+'__'+dfRank4bokeh1['Metadata_mmoles_per_liter'].astype(str)+"   "+"p-val: "+dfRank4bokeh1['phenotype_abundance_pval'].astype(str)#+'__'+dfRank4bokeh1['Metadata_Location']


DEFAULT_TICKERS=dfRank4bokeh1['ticker1'].tolist()

ticker1 = Select(title="Sample-Dose Ranked List", value=DEFAULT_TICKERS[0], options=DEFAULT_TICKERS)
stats2 = PreText(text='', width=600)
###################


# dataPath='/home/ubuntu/bucket/projects/2018_04_20_Rosetta/workspace/results/sample_cells/results/'


# # DATA_DIR=dataPath
# imageslist0=os.listdir(dataPath)
# imageslist=[s for s in imageslist0 if 'synth' in s or 'real' in s]

# n_rows=len(imageslist)

# final_df=pd.DataFrame(index=range(n_rows), columns=['im_list','user_inputs','true_class'])
# final_df['im_list']=imageslist


n_rows=1;

@lru_cache()
def read_rgb(imagePath):
    img = skimage.io.imread(imagePath)
    img=np.flip(img,0)
    rgb_img = img.astype(np.uint8)
    return rgb_img

@lru_cache()
def read_rgb2(imagePath):
#     img = skimage.io.imread(imagePath)
    img = Image.open(imagePath)
    # print(img.shape,img.max(),img.min())
    img = img.convert('RGBA')
    img = np.array(img)
    img=np.flip(img,0)
    rgb_img = img.astype(np.uint8)
    return rgb_img


# @lru_cache()
# def get_images():
#     print('get_images')

#     imss=[]
#     for i in range(n_rows):
# #                print(i,len(clusNum),clusNum[i])
#         imss.append(read_rgb2(dataPath+'/'+imageslist[i]))   
# #         print(read_rgb(dataPath+'/'+imageslist[i]).shape)
# #         imss.append(read_rgb(folder_path+'/cluster'+str(clusNum[i])+'_barImpFeatures.png'))
                
#     return imss


# @lru_cache()
def get_images(df_to_disp):
    print('get_images')

    wss=['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P']
    imss=[]
    for r in range(df_to_disp.shape[0]):
        p,w=df_to_disp.loc[r,["Metadata_Plate","Metadata_Well"]].values
    #     print(p,w)
        # select a random site
        site=10

        row=format(wss.index(w[0])+1, '02d')
        wellCode='r'+row+'c'+w[1:]    

        ch=1
        img_all_chans=[]
#         for ch in range(1,6):
#             png_addr=ImFilesRoot+"/"+p+"/"+wellCode+"f08p01-ch"+str(ch)+"sk1fk1fl1.png"
#             img_all_chans.append(read_rgb2(png_addr)[:,:,np.newaxis]) 
# #             img_all_chans.append(skimage.io.imread(png_addr)[:,:,np.newaxis])

#         all_chans_arr=np.concatenate(img_all_chans,axis=2)
#         rgb_im=visualize_data.CP_to_RGB_single(all_chans_arr)

        png_addr=ImFilesRoot+"/"+p+"/"+wellCode+"f08p01-ch"+str(5)+"sk1fk1fl1.png"
        rgb_im=read_rgb2(png_addr)
        imss.append(rgb_im)
                
    return imss

# # set up widgets
# nameBox=TextInput(title="Name", value="",width=200,height=40)
# bt = Button(label='Save as CSV',width=200,height=40)

# set up plots
def create_figure(df_to_disp):
    images = get_images(df_to_disp)
    kw = dict()
#     #kw['title'] = "%s vs %s" % (x_title, y_title)
#     N, M=int(480/2), int(640/2)
    N, M=int(1000), int(1000)
#     N2, M2=int(1200/2), int(800/2)
    N2, M2=int(1000), int(1000)
    kw['x_range'] = (10,M-10)
    kw['y_range'] = (10,N-10)
    kw['plot_width'] = M
    kw['plot_height'] = N
#     kw['title'] = 'Clustering'
#        kws=[kw]*len(images)
    kws=[];
    for i in range(0,len(images)):
#            print(i,clusterss[i])
#            print(kws[i])
        kw1=kw.copy()
#         kw1['title']='Cluster - '+str(clusterss[i])

        kw2=kw.copy()
#         kw2['title']='Cluster - '+str(clusterss[i])
        kw2['plot_width'] = M2
        kw2['plot_height'] = N2
        kws.append(kw2);
#         kws.append(kw1);
            
#            print(kws)
#        print(len(kws),len(images),len(clusterss))
    Ms=[M]*len(images)
    Ns=[N]*len(images)
#        kws=[kw,kw,kw,kw,kw]        
#        Ms=[M,M,M,M,M];Ns=[N,N,N,N,N]
        
    ts=[];
#     notes=[]
    for l in range(len(images)):   
        new_data1 = dict();        
        ts1 = figure(tools='pan,box_zoom,reset', **kws[l])
        ts1.axis.visible = False
        ts1.xgrid.visible = False
        ts1.ygrid.visible = False
#        print(images[l].shape,Ms[l],Ns[l])
        r1=ts1.image_rgba(image=[], x=0, y=0, dw=Ms[l], dh=Ns[l])
        ds1 = r1.data_source;
        new_data1['image'] = [images[l]]
        ds1.data = new_data1
        
#         if 'real' in imageslist[l]:
#             title0="Write down suggested class of each cell in the below figure:"
#         else:
#             title0="Write down suggested class of each cell (except first cell) in the below figure:"
            
#         nts=TextInput(title=title0, value="",width=500,height=40)
#         notes.append(nts);
#         ts.append(column(children=[nts,ts1]));
        ts.append(column(children=[ts1]));
        
#         ts.append(ts1);
        
#     if t3=='Subpopulation Analysis':
#     colOfRows=[row(children=ts[0:1])];
# #        print(len(ts))
#     for c in range(1,len(ts),2):
# #            print(c)
#         colOfRows.append(row(children=ts[c:c+2]))        
# #        colOfRows
#     t3SpecificLayout=column(children=colOfRows)
#     else:    
    t3SpecificLayout=column(children=ts)
        
    return t3SpecificLayout




# # # set up callbacks
# def change_click():
#     print('Save was clicked')
#     for i in range(n_rows):
#         print(notes[i].value);
#         final_df.loc[i,'user_inputs']=str(notes[i].value)
#         final_df.loc[i,'true_class']=''.join(imageslist[i].split('_')[1].split('.')[0].split('-'))
    
#     user_name=nameBox.value
#     final_df.to_csv(dataPath+user_name+'_eval.csv',index=False)
    
# bt.on_click(change_click)


def ticker1_change(attrname, old, new):
    """ Loads WT-MT data again and starts ... """
    
    tickVal = ticker1.value
    X_Metadata_broad_sample=tickVal.split('__')[0]
    X_Metadata_mmoles_per_liter=tickVal.split('__')[1].split(' ')[0]  
    
    
    df_to_disp=meta_lincs2[(meta_lincs2["Metadata_broad_sample"]==X_Metadata_broad_sample) &\
               (meta_lincs2["Metadata_mmoles_per_liter"]==float(X_Metadata_mmoles_per_liter))].reset_index(drop=True)

    fS=create_figure(df_to_disp)
    
    
    stats2_df=drug_list_rank[(drug_list_rank["Metadata_broad_sample"]==X_Metadata_broad_sample) &\
           (drug_list_rank["Metadata_mmoles_per_liter"]==float(X_Metadata_mmoles_per_liter))].T
    
    stats2.text = str(stats2_df)
    series.children=[column(children=[ticker1,stats2])]+[fS]

    
ticker1.on_change('value', ticker1_change)


X_Metadata_broad_sample=ticker1.value.split('__')[0]
X_Metadata_mmoles_per_liter=ticker1.value.split('__')[1].split(' ')[0]     


df_to_disp=meta_lincs2[(meta_lincs2["Metadata_broad_sample"]==X_Metadata_broad_sample) &\
           (meta_lincs2["Metadata_mmoles_per_liter"]==float(X_Metadata_mmoles_per_liter))].reset_index(drop=True)

fS=create_figure(df_to_disp)

stats2_df=drug_list_rank[(drug_list_rank["Metadata_broad_sample"]==X_Metadata_broad_sample) &\
       (drug_list_rank["Metadata_mmoles_per_liter"]==float(X_Metadata_mmoles_per_liter))].T
stats2.text = str(stats2_df)
# set up layout
# fS,notes=create_figure()
# fS=[]

# series = row(children=[ticker1]+[fS])
series = row(children=[column(children=[ticker1,stats2])]+[fS])
layoutt = column(series)


# initialize
# update()

curdoc().add_root(layoutt)
# curdoc().title = "Rosetta Eval"

