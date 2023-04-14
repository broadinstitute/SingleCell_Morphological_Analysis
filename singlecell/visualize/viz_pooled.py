"""
@author: mhaghigh
"""
import pandas as pd
import numpy as np
import seaborn as sns

# from sklearn import preprocessing
# import pickle
import matplotlib.pyplot as plt
import os
from functools import reduce
from sklearn.cluster import KMeans

# import skimage
from skimage import exposure
from skimage.transform import resize
import skimage.io
from skimage.segmentation import flood_fill
from . import crop_cell
from . import cp_to_rgb

# from . import read_from_gallery
from ..read import read_from_gallery
from skimage.util import img_as_ubyte, img_as_float


# def visualize_n_SingleCell(channels , sc_df , boxSize , outline = False, color=False,\
#                            info_columns=[], title="",compressed=False,compressed_im_size=None):


def visualize_n_SingleCell_pooled2(
    channels,
    sc_df,
    boxSize,
    im_size,
    outline=False,
    color=False,
    info_columns=[],
    title="",
    neigh_center_ls=[],
):
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

    halfBoxSize = int(boxSize / 2)
    #     print(channels)

    #     plt.ioff()
    f, axarr = plt.subplots(
        sc_df.shape[0], len(channels), figsize=(len(channels) * 2, sc_df.shape[0] * 2)
    )
    if len(title) > 0:
        print(title)
        f.suptitle(title)

    f.subplots_adjust(hspace=0, wspace=0)

    align_column_ch_name_map = {
        "DNA": "DAPI_Painting",
        "ER": "ConA",
        "Mito": "Mito",
        "Phalloidin": "Phalloidin",
        "WGA": "WGA",
    }

    for index in range(sc_df.shape[0]):

        xCenter = int(sc_df.loc[index, "Nuclei_Location_Center_X"])
        yCenter = int(sc_df.loc[index, "Nuclei_Location_Center_Y"])

        cpi = 0
        for c in channels:
            if c == "Outline":
                imPath = sc_df.loc[index, "Path_Outlines"]
                #                 im_size=sc_df["Width_CorrDNA"].values[0]   #cp220
                #                 im_size= 5500
                # print(imPath)

                #                 if os.path.exists(imPath):
                #                 im_out=read_from_gallery.read_image(imPath)
                ov_im = read_resize_overlay_pooled(imPath, im_size)

                cell_outline = extract_cell_outline(ov_im)
                colored_masks = flood_fill_cells(
                    cell_outline, [xCenter], [yCenter], color_new=0.7
                )

                if neigh_center_ls:
                    #                     print(neigh_center_ls[index][0],neigh_center_ls[index][1])
                    colored_masks = flood_fill_cells(
                        colored_masks,
                        neigh_center_ls[index][0],
                        neigh_center_ls[index][1],
                        color_new=0.6,
                    )
                imD2 = colored_masks[
                    yCenter - halfBoxSize : yCenter + halfBoxSize,
                    xCenter - halfBoxSize : xCenter + halfBoxSize,
                ]

                #                 else:
                #                     imD2 = flood_fill_cells(cell_outline,neigh_center_ls[index][0],neigh_center_ls[index][1])\
                #                     [yCenter-halfBoxSize:yCenter+halfBoxSize,\
                #                                                xCenter-halfBoxSize:xCenter+halfBoxSize]

                #                 print(im_out.max)
                #                     imD2 = resize(skimage.io.imread(imPath),(im_size,im_size),\
                #                 imD2 = resize(im_out,(im_size,im_size),\
                #                                  mode='constant',preserve_range=True,order=0)\

                #                 print(imD2.max)
                #                 else:
                #                     imD2=np.zeros((boxSize,boxSize))
                clim_max = imD2.max()
            else:
                #                 ch_D=sc_df.loc[index,'Image_FileName_Orig'+c];
                ch_D = sc_df.loc[index, "FileName_Corr" + c]
                #                 print(ch_D)
                #         imageDir=imDir+subjectID+' Mito_Morphology/'
                #                 imageDir=sc_df.loc[index,'Image_PathName_Orig'+c]+'/'
                imageDir = sc_df.loc[index, "PathName_Corr" + c] + "/"
                imPath = imageDir + ch_D
                #                 print(imPath)
                #                 print(yCenter,xCenter)

                xCenterC = int(sc_df.loc[index, "Nuclei_Location_Center_X"]) + int(
                    sc_df.loc[index, "Align_Xshift_" + align_column_ch_name_map[c]]
                )
                yCenterC = int(sc_df.loc[index, "Nuclei_Location_Center_Y"]) + int(
                    sc_df.loc[index, "Align_Yshift_" + align_column_ch_name_map[c]]
                )

                #                 print(yCenterC,xCenterC)
                #                 print(imPath)
                #                 print('im_size',np.squeeze(skimage.io.imread(imPath)).shape)

                #                 imD=np.squeeze(skimage.io.imread(imPath))
                print(imPath)
                imD = np.squeeze(read_from_gallery.read_image(imPath))

                imD1 = exposure.rescale_intensity(
                    imD, in_range=(imD.min(), np.percentile(imD, 99.95))
                )  # output is  unit16
                #                 and scaled to range 0,65535
                clim_max = imD1.max()
                imD2 = crop_cell.crop_single_cell_image(
                    imD1, xCenterC, yCenterC, halfBoxSize
                )

            #                 imD2=imD1[yCenterC-halfBoxSize:yCenterC+halfBoxSize,xCenterC-halfBoxSize:xCenterC+halfBoxSize]
            #                 print(imD1.min(),imD1.max())
            #                 print(np.squeeze(skimage.io.imread(imPath)).shape)
            #             axarr[index,cpi].imshow(imD,cmap='gray',clim=(0, maxRanges[c]));axarr[0,cpi].set_title(c);
            #             print(imD.shape,'h')
            axarr[index, cpi].imshow(imD2, cmap="gray", clim=(0, clim_max))
            axarr[0, cpi].set_title(c)
            axarr[index, cpi].axes.xaxis.set_ticks([])
            axarr[index, cpi].axes.yaxis.set_ticks([])
            cpi += 1

        #         Well=sc_df.loc[index,'Metadata_Well']
        #         Site=str(sc_df.loc[index,'Metadata_Site'])
        #         imylabel=Well+'\n'+Site

        if "label" in sc_df.columns:
            imylabel = (
                sc_df.loc[index, "label"]
                + "\n"
                + sc_df.loc[index, "Metadata_Foci_Barcode_MatchedTo_Barcode"][0:9]
            )
        else:
            imylabel = sc_df.loc[index, "Metadata_Foci_Barcode_MatchedTo_Barcode"][0:12]
        #         print(imylabel)
        axarr[index, 0].set_ylabel(imylabel)

    #     plt.ion()
    return f


def visualize_n_SingleCell_pooled(
    channels,
    sc_df,
    boxSize,
    im_size,
    outline=False,
    color=False,
    info_columns=[],
    title="",
    neigh_center_ls=[],
):
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

    halfBoxSize = int(boxSize / 2)
    #     print(channels)

    #     plt.ioff()

    columns_count = len(channels) + int(outline) + int(color)
    rows_count = sc_df.shape[0]

    f, axarr = plt.subplots(
        rows_count, columns_count, figsize=(columns_count * 2, rows_count * 2)
    )
    if len(title) > 0:
        print(title)
        f.suptitle(title)

    f.subplots_adjust(hspace=0, wspace=0)

    align_column_ch_name_map = {
        "DNA": "DAPI_Painting",
        "ER": "ConA",
        "Mito": "Mito",
        "Phalloidin": "Phalloidin",
        "WGA": "WGA",
    }

    for index in range(rows_count):

        xCenter = int(sc_df.loc[index, "Nuclei_Location_Center_X"])
        yCenter = int(sc_df.loc[index, "Nuclei_Location_Center_Y"])

        sc_collage_row = np.zeros((boxSize, boxSize, columns_count))

        clim_max = []
        for ci in range(len(channels)):
            if index == 0:
                axarr[index, ci].set_title(channels[ci])

            ch_fName = sc_df.loc[index, "FileName_Corr" + channels[ci]]
            ch_pName = sc_df.loc[index, "PathName_Corr" + channels[ci]]

            xCenterC = xCenter + int(
                sc_df.loc[
                    index, "Align_Xshift_" + align_column_ch_name_map[channels[ci]]
                ]
            )
            yCenterC = yCenter + int(
                sc_df.loc[
                    index, "Align_Yshift_" + align_column_ch_name_map[channels[ci]]
                ]
            )

            imPath = ch_pName + "/" + ch_fName
            # print(imPath)
            image = np.squeeze(read_from_gallery.read_image(imPath))
            if 1:
                image = exposure.rescale_intensity(
                    image, in_range=(image.min(), np.percentile(image, 99.95))
                )

            clim_max.append(image.max())
            image_cropped = crop_cell.crop_single_cell_image(
                image, xCenterC, yCenterC, halfBoxSize
            )
            # print(image_cropped.shape, xCenter, yCenter, image_cropped.max())

            sc_collage_row[:, :, ci] = image_cropped

        for c in range(len(channels)):
            axarr[index, c].imshow(
                sc_collage_row[:, :, c],
                interpolation=None,
                cmap="gray",
                clim=(0, clim_max[c]),
            )
        #             axarr[index,c].axes.xaxis.set_ticks([]);axarr[index,c].axes.yaxis.set_ticks([])

        if color:
            #             print('not available here yet')
            color_im, colorImagesList = cp_to_rgb.CP_to_RGB_single(
                sc_collage_row[:, :, : len(channels)], channels
            )
            #             print(color_im.max(),color_im.min())
            # for c in range(len(channels)):
            #     axarr[index, c].imshow(colorImagesList[c], interpolation=None)

            axarr[index, c + 1].imshow(color_im)
            #             axarr[index,c+1].axes.xaxis.set_ticks([]);axarr[index,c+1].axes.yaxis.set_ticks([]);

            if index == 0:
                axarr[index, c + 1].set_title("composite")

        if outline:

            imPath = sc_df.loc[index, "Path_Outlines"]
            # print(imPath)
            #           im_out=read_from_gallery.read_image(imPath)
            ov_im = read_resize_overlay_pooled(imPath, im_size)

            cell_outline = extract_cell_outline(ov_im)

            if neigh_center_ls:
                #                     print(neigh_center_ls[index][0],neigh_center_ls[index][1])
                cell_outline = flood_fill_cells(
                    cell_outline,
                    neigh_center_ls[index][0],
                    neigh_center_ls[index][1],
                    color_new=0.4,
                )

            colored_masks = flood_fill_cells(
                cell_outline, [xCenter], [yCenter], color_new=0.9
            )

            imD2 = colored_masks[
                yCenter - halfBoxSize : yCenter + halfBoxSize,
                xCenter - halfBoxSize : xCenter + halfBoxSize,
            ]

            clim_max = imD2.max()
            #             axarr[index,cpi].imshow(imD2,cmap='gray',clim=(0, clim_max));axarr[0,cpi].set_title(c);
            axarr[index, columns_count - 1].imshow(
                imD2, cmap="gray", clim=(0, clim_max)
            )
            axarr[0, columns_count - 1].set_title("cell outline")
        #             axarr[index,columns_count-1].axes.xaxis.set_ticks([]);axarr[index,columns_count-1].axes.yaxis.set_ticks([]);

        for c in range(columns_count):
            axarr[index, c].axes.xaxis.set_ticks([])
            axarr[index, c].axes.yaxis.set_ticks([])
        ## add info in the columns specified in the input next to to each single cell
        #         if info_columns:
        #             image_y_label = form_image_y_label_string(sc_df.loc[index,info_columns],info_columns);
        #             axarr[index,0].set_ylabel(image_y_label);

        if "label" in sc_df.columns:
            imylabel = (
                sc_df.loc[index, "label"]
                + "\n"
                + sc_df.loc[index, "Metadata_Foci_Barcode_MatchedTo_Barcode"][0:9]
            )
        else:
            imylabel = sc_df.loc[index, "Metadata_Foci_Barcode_MatchedTo_Barcode"][0:12]
        #         print(imylabel)
        axarr[index, 0].set_ylabel(imylabel)

    #     plt.ion()
    return f


def read_resize_overlay_pooled(overlay_dir, orig_im_w):
    """
    This function read and resizes overlay to the size of the original image
    """

    #     import skimage.io
    #     from skimage.transform import rescale, resize, downscale_local_mean

    ov_im = resize(
        read_from_gallery.read_image(overlay_dir),
        (orig_im_w, orig_im_w),
        mode="constant",
        preserve_range=True,
        order=0,
    )

    #     ov_im=resize(skimage.io.imread(overlay_dir), (orig_im_w,orig_im_w),\
    #                  mode='constant',preserve_range=True,order=0).astype('uint8')
    return img_as_ubyte(ov_im)


def extract_cell_outline(ov_im):
    """
    This function read and resizes overlay to the size of the original image
    """

    #     import skimage.io
    #     from skimage.transform import rescale, resize, downscale_local_mean

    #     ov_im = resize(read_from_gallery.read_image(overlay_dir),(orig_im_w,orig_im_w),\
    #                                  mode='constant',preserve_range=True,order=0)
    if ov_im.dtype == np.uint8:
        cell_color = (255, 255, 255)
    elif ov_im.dtype == np.uint16:
        cell_color = (65535, 65535, 65535)
    elif ov_im.dtype == np.float32:
        cell_color = (1, 1, 1)

    cell_bound = np.copy(ov_im)
    indices_not_w = np.where(~np.all(cell_bound == cell_color, axis=-1))
    cell_bound[indices_not_w] = 0

    colored_cells = cell_bound[:, :, 0]

    colored_cells = img_as_float(colored_cells)
    return colored_cells


def flood_fill_cells(colored_cells, cent_x_ls, cent_y_ls, color_new=0.7):

    for p in range(len(cent_x_ls)):
        cent_x, cent_y = cent_x_ls[p], cent_y_ls[p]
        colored_cells = flood_fill(
            colored_cells, (cent_y, cent_x), color_new, connectivity=1
        )

    return colored_cells


def colormask_cells_to_parent_cell_number(ov_im, cell_centers_df):
    """
    given an input overlay image, and a mapping of cell centers to object numbers
        -> color each cell to its object number

    Inputs:
       ov_im: overlay image, is an RGB image with outline of nuclei and cells in the first channel
       cell_centers_df: a dataframe containing cell_centers and map to object numbers
               required columns in the df:
                  - Parent_Cells
                  - Nuclei_Location_Center_X , Nuclei_Location_Center_Y

    """

    from skimage.segmentation import flood_fill

    cell_color = (65535, 65535, 65535)
    cell_bound = np.copy(ov_im)
    indices_not_w = np.where(~np.all(cell_bound == cell_color, axis=-1))
    cell_bound[indices_not_w] = 0

    colored_cells = cell_bound[:, :, 0].astype(int)

    parent_cells = cell_centers_df.Parent_Cells.unique().tolist()
    for p in parent_cells:
        cent_x, cent_y = cell_centers_df.loc[
            cell_centers_df["Parent_Cells"] == p,
            ["Nuclei_Location_Center_X", "Nuclei_Location_Center_Y"],
        ].values[0]
        colored_cells = flood_fill(colored_cells, (cent_y, cent_x), p, connectivity=1)

    return colored_cells


def extract_neighbor_cells_center(df_p_s, df_samples, max_dist_of_neigh):
    neigh_center_ls = []
    for i in range(df_samples.shape[0]):
        df_samples_0 = df_samples.loc[i]

        df_samples_0_neighs = df_p_s[
            (df_p_s["PathName_CorrER"] == df_samples_0["PathName_CorrER"])
            & (
                df_p_s["Metadata_Foci_Barcode_MatchedTo_Barcode"]
                == df_samples_0["Metadata_Foci_Barcode_MatchedTo_Barcode"]
            )
        ]

        df_s0_neighs = df_samples_0_neighs[
            (
                abs(
                    df_samples_0_neighs["Nuclei_Location_Center_X"]
                    - df_samples_0["Nuclei_Location_Center_X"]
                )
                < max_dist_of_neigh
            )
            & (
                abs(
                    df_samples_0_neighs["Nuclei_Location_Center_Y"]
                    - df_samples_0["Nuclei_Location_Center_Y"]
                )
                < max_dist_of_neigh
            )
        ]

        x_cent_ls = df_s0_neighs.Nuclei_Location_Center_X.astype(
            int
        ).tolist()  # +[df_samples_0['Nuclei_Location_Center_X']]
        y_cent_ls = df_s0_neighs.Nuclei_Location_Center_Y.astype(
            int
        ).tolist()  # +[df_samples_0['Nuclei_Location_Center_Y']]

        print(df_samples_0["Nuclei_Location_Center_Y"])
        x_cent_ls.remove(df_samples_0["Nuclei_Location_Center_X"].astype(int))
        y_cent_ls.remove(df_samples_0["Nuclei_Location_Center_Y"].astype(int))

        neigh_center_ls.append([x_cent_ls, y_cent_ls])
    return neigh_center_ls
