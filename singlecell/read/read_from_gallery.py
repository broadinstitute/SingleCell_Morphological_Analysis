from io import BytesIO
from matplotlib import image as mpimg
import boto3
import pyarrow.parquet as pq
import pandas as pd
import sqlite3
from functools import reduce

s3_client = boto3.client("s3")
bucket_name = "cellpainting-gallery"
compartments = ["cells", "cytoplasm", "nuclei"]


def read_sqlite(file_name, wells_subset_param="all"):
    #     chunk_size = 1000000
    #     offset = 0
    #     while True:

    #     local_file = '/tmp/tmp_loc.sqlite'
    #     s3.download_file(bucket_name, file_name, local_file)

    data = s3_client.get_object(Bucket=bucket_name, Key=file_name)

    buffer = BytesIO(data["Body"].read())
    conn = sqlite3.connect(buffer)

    img_query = "select * from {}".format("Image")

    if wells != "all":
        wells, meta_well_col_str = wells_subset_param
        list_str = "('"
        for i in wells:
            list_str = list_str + str(i) + "','"
        list_str = list_str[:-2] + ")"

        im_query_well_subset = " WHERE {} IN {};".format(meta_well_col_str, list_str)

        img_query = img_query + im_query_well_subset

    plateImageDf = pd.read_sql_query(img_query, conn)

    if wells != "all":
        img_nums = plateImageDf.ImageNumber.unique().tolist()

        img_nums_str = "("
        for i in img_nums:
            img_nums_str = img_nums_str + str(i) + ","
        img_nums_str = img_nums_str[:-1] + ")"

    plateDf_list = []
    for compartment in compartments:
        compartment_query = "select * from {}".format(compartment)
        if wells != "all":
            compartment_query = compartment_query + " WHERE {} IN {};".format(
                "ImageNumber", img_nums_str
            )

        #         plateDf_list.append(pd.read_sql_query(compartment_query, conn))
        plateDf_list.append(read_chunk(compartment_query, conn))

    plateDf = reduce(
        lambda left, right: pd.merge(
            left, right, on=["TableNumber", "ImageNumber", "ObjectNumber"]
        ),
        plateDf_list,
    )

    plateDfwMeta = pd.merge(plateDf, plateImageDf, on=["TableNumber", "ImageNumber"])
    plateDfwMeta = plateDfwMeta.loc[:, ~plateDfwMeta.columns.duplicated()]

    return plateDfwMeta


# def read_chunk(file_name, wells_subset_param = 'all'):
#     # create a sqlite connection
#     conn = sqlite3.connect(':memory:')

#     # create an iterator to read the file in chunks
#     chunk_size = 1000000
#     offset = 0
#     while True:
#         # use the client to download the chunk of the file to a BytesIO object
#         obj = s3.get_object(Bucket=bucket_name, Key=file_name,
#                             Range=f'bytes={offset}-{offset+chunk_size-1}')
#         buffer = BytesIO(obj['Body'].read())

#         # read the chunk into a temporary SQLite database
#         with sqlite3.connect(buffer) as chunk_conn:
#             for chunk in pd.read_sql_query(query, chunk_conn, chunksize=10000):
#                 # process the chunk as needed
#                 # for example, you can append it to a list of dataframes
#                 df_ls.append(chunk)
#                 # or write it to a new SQLite database
#                 pass

#         # if the chunk was smaller than the requested chunk size, we've reached the end of the file
#         if len(buffer.getbuffer()) < chunk_size:
#             break

#         # update the offset to read the next chunk
#         offset += chunk_size

#     # concatenate the results of the processing into a final pandas dataframe
#     df = pd.concat(df_ls)
#     return df


def read_parquet(file_name):
    # use the client to download the file to a temporary file object
    tmp = BytesIO()
    s3_client.download_fileobj(bucket_name, file_name, tmp)

    # read the file from the temporary file object into a pandas dataframe
    tmp.seek(0)
    df = pd.read_parquet(tmp)

    return df


def read_csv(path):
    data = s3_client.get_object(
        Bucket=bucket_name, Key=path
    )  # Key="/".join(image_url.split("/")[3:])
    df_p_s = pd.read_csv(data["Body"])

    return df_p_s


def read_csv_gzip(path):
    data = s3_client.get_object(
        Bucket=bucket_name, Key=path
    )  # Key="/".join(image_url.split("/")[3:])
    df_p_s = pd.read_csv(data["Body"], compression="gzip")

    return df_p_s


def read_image(path):
    #     import os
    #     import requests

    response = s3_client.get_object(Bucket=bucket_name, Key=path)

    #     example Key="cpg0021-periscope/broad/images/20200805_A549_WG_Screen/images_corrected_cropped/CP186N_Well3/CorrDNA/CorrDNA_Site_95.tiff"
    format_ext = path.split(".")[-1]
    image = mpimg.imread(BytesIO(response["Body"].read()), format=format_ext)

    return image
