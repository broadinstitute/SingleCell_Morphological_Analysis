from io import BytesIO
from matplotlib import image as mpimg
import boto3
import pyarrow.parquet as pq
import pandas as pd

s3_client = boto3.client("s3")



def read_parquet(file_name):

    # specify the S3 bucket and file you want to read
    bucket_name = "cellpainting-gallery"

    # use the client to download the file to a temporary file object
    tmp = io.BytesIO()
    s3.download_fileobj(bucket_name, file_name, tmp)

    # read the file from the temporary file object into a pandas dataframe
    tmp.seek(0)
    df = pd.read_parquet(tmp)    
    
    return df


def read_csv(path):

    data = s3_client.get_object(
        Bucket="cellpainting-gallery",\
        Key=path)# Key="/".join(image_url.split("/")[3:])
    df_p_s = pd.read_csv(data['Body'])    
    
    return df_p_s


def read_csv_gzip(path):

    data = s3_client.get_object(
        Bucket="cellpainting-gallery",\
        Key=path)# Key="/".join(image_url.split("/")[3:])
    df_p_s = pd.read_csv(data['Body'], compression='gzip')    
    
    return df_p_s



def read_image(path):
#     import os
#     import requests

    response = s3_client.get_object(
        Bucket="cellpainting-gallery", Key=path)
    
#     example Key="cpg0021-periscope/broad/images/20200805_A549_WG_Screen/images_corrected_cropped/CP186N_Well3/CorrDNA/CorrDNA_Site_95.tiff"
    format_ext=path.split('.')[-1]
    image = mpimg.imread(BytesIO(response["Body"].read()), format=format_ext)
    
    return image