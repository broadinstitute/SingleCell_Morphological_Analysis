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
#     img=abs(img.min())+img
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


# def normalize(img):
#     # normalize to [0,1]
#     percentile = 99.95
#     high = np.percentile(img, percentile)
#     low = np.percentile(img, 100-percentile)

# #     img = np.minimum(high, img)
# #     img = np.maximum(low, img)
    
# #     high=np.max(img)
# #     low=np.min(img)

# #     img = (img - low) / (high - low) # gives float64, thus cast to 8 bit later
# #     vmin, vmax = scipy.stats.scoreatpercentile(image, (0.05, 99.95))
# #     vmax = min(vmax, pmax)
# #     image = skimage.exposure.rescale_intensity(img, in_range=(low, high))
#     image = skimage.exposure.rescale_intensity(img, in_range=(low, high),out_range=(0,255)).astype(np.uint8)
# #     image = skimage.exposure.rescale_intensity(img, in_range=(-1, 1))
# #     image[image>1]=1
# #     image[image<-1]=-1
# #     print(image.max(),image.min())
# #     img = skimage.img_as_ubyte(image)
# #     print(image.max(),image.min())
#     return imageJab