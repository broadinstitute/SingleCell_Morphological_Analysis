"""
@author: mhaghigh

* edge cells:
  - Are the cells which their nuclei center is less than a specified margin,
    - The margin could be cells an optional input or if not specified, 
    - By default, it is the 90th percentile of the "Cells_AreaShape_MajorAxisLength" feature devided by 2

"""

import pandas as pd
import numpy as np

########################function to remove cells on the border
def edgeCellFilter(df_input):   
    # remove cells on the border
#     imgSize=1080
#     imgSize=2048
    imgSize=df_input.Metadata_ImageSizeX.values[0]
    edge_margin=int(np.percentile(df_input.Cells_AreaShape_MajorAxisLength.values, 90)/2);
    df_1_centCells=df_input.loc[~((df_input['Nuclei_Location_Center_X']>(imgSize-edge_margin)) | \
                              (df_input['Nuclei_Location_Center_X']<(edge_margin))\
                            | (df_input['Nuclei_Location_Center_Y']>(imgSize-edge_margin)) | \
                              (df_input['Nuclei_Location_Center_Y']<(edge_margin))),:].reset_index(drop=True)
    return df_1_centCells,edge_margin



########################function to remove cells on the border
def edgeCellFilter2(df_1,imgSize,edge_margin):   

    df_1_centCells=df_1.loc[~((df_1['Nuclei_Location_Center_X']>(imgSize-edge_margin)) | \
                              (df_1['Nuclei_Location_Center_X']<(edge_margin))\
                            | (df_1['Nuclei_Location_Center_Y']>(imgSize-edge_margin)) | \
                              (df_1['Nuclei_Location_Center_Y']<(edge_margin))),:].reset_index(drop=True)
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