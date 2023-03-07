import numpy as np
import scipy.spatial
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from random import sample, choices
from scipy.stats import pearsonr


def abundance_rawFeature(inDf, groupingColumns, featureName, thrsh):

    """
    Calculates abundance of an input feature in a well/site. Defiend as:
    Ratio of Numebr of single cells in that well/site that their feature value is more than a certain
    threshold over total number of single cells in that well/site.

    How the threshold is usually defined:
    - For data in different plates, for each plate threshold can be 90 percentile of that feature across
      all single cells in that plate
    - For patient data: 90 percentile of the feature in all patient single cells


    This function takes the input dataframe and output/plot Abundance of featureName.

    Parameters:
    inDf   (pandas df): input dataframe
    featureName  (str): The column name corresponding to the feature we want to calculate its abundance
    groupingColumns (list): list of columns used for grouping single cells
    thrsh (float):

    Returns:
    abundanceRes (pandas df): as defined above

    """
    abundanceRes = (
        inDf.groupby(groupingColumns)
        .apply(
            lambda x: np.sum(x[featureName].values >= thrsh) / (x[featureName].shape[0])
        )
        .reset_index()
    )

    abundanceRes = abundanceRes.rename(columns={0: "Abundance"})
    return abundanceRes


def abundance_classificationProbabiltyOutputs(
    inDf, groupingColumns, featureName, thrsh
):

    """
    Calculates abundance of an input feature in a well/site. Defiend as:
    Ratio of Numebr of single cells in that well/site that their class prediction probability is more than a certain
    threshold over number of certain+uncertain number of samples.

    How the threshold is usually defined:
    - Based on how certain we need to be on the phenotype


    This function takes the input dataframe and output/plot Abundance of featureName.

    Parameters:
    inDf   (pandas df): input dataframe
    featureName  (str): The column name corresponding to the feature we want to calculate its abundance
    groupingColumns (list): list of columns used for grouping single cells
    thrsh (float):  range:[0,1] probaility that a of certian decision

    Returns:
    abundanceRes (pandas df): as defined above

    """
    #     abundanceResCon=df_T.groupby([groupingColumns]).apply(lambda x: \
    #      np.sum(x[featureName]values>=thrsh)/np.sum((x[featureName].values>=thrsh)|(x[featureName].values<=1-thrsh))).reset_index();

    abundanceResCon = (
        df_T.groupby(["subject", "label"])
        .apply(
            lambda x: np.sum(x.P1.values >= thrsh)
            / np.sum((x.P1.values >= thrsh) | (x.P1.values <= 1 - thrsh))
        )
        .reset_index()
    )

    abundanceRes = abundanceRes.rename(columns={0: "Abundance"})
    return abundanceRes
