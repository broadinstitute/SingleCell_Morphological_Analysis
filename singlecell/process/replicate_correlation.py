"""
@author: mhaghigh
"""

import numpy as np
import scipy.spatial
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from random import sample, choices
from scipy.stats import pearsonr
import WeightedCorr
import itertools as it

# sns.set_style("whitegrid")
sns.set(rc={"lines.linewidth": 2})

# from singlecell.process.replicate_correlation import replicate_null_corr_coefs


def replicate_null_corr_coefs(
    inDf, pertColName, featColNames, plotEnabled, title="", hist_bins=100
):

    """
    Calculates replicate correlation versus across purtburtion correlations

    This function takes the input dataframe and output/plot replicate correlations.

    Parameters:
    inDf   (pandas df): input dataframe contains metadata and features
    pertColName  (str): The column based on which we define replicates of a purturbation
    featColNames(list): The list of all columns corresponding to features
    plotEnabled (bool): If True or 1, plots the curves

    Returns:
    repCorrDf   (list):

    """

    #     inDf=df_meanPrepCorr.copy()
    #     pertColName='Metadata_Sample_uniqe'
    #     featColNames=cpFeats_P_reduced

    df = inDf.copy()

    # df[featColNames]=inDf[featColNames].interpolate();
    uniqPert = df[pertColName].unique().tolist()
    repC = []

    repCorrDf = pd.DataFrame(index=uniqPert, columns=["RepCor"])

    repSizeDF = df.groupby([pertColName]).size().reset_index()
    highRepComp = repSizeDF[repSizeDF[0] > 1][pertColName].tolist()

    #     corr_mat_df=df.set_index([pertColName,'Metadata_Plate']).loc[:,featColNames].T.corr()
    corr_mat_df = df.set_index(pertColName).loc[:, featColNames].T.corr()

    rep_corr_ls = []
    null_corr_ls = []

    for u in highRepComp:

        not_u = highRepComp[:]
        not_u.remove(u)

        repCorrPurtbs = corr_mat_df.loc[u, u]
        repCorr = np.nanmean(
            repCorrPurtbs.values[np.triu_indices(repCorrPurtbs.shape[0], k=1)]
        )
        repCorrDf.loc[u, "RepCor"] = repCorr

        #         corr_rest=corr_mat_df.loc[u,not_u].mean().mean()
        corr_rest = (
            corr_mat_df.loc[not_u, u].sample(repCorrPurtbs.shape[0]).mean().mean()
        )

        rep_corr_ls.append(repCorr)
        null_corr_ls.append(corr_rest)

    rep_corr_ls = [
        rep_corr_ls for rep_corr_ls in rep_corr_ls if str(rep_corr_ls) != "nan"
    ]
    null_corr_ls = [
        null_corr_ls for null_corr_ls in null_corr_ls if str(null_corr_ls) != "nan"
    ]

    perc95 = np.percentile(null_corr_ls, 90)
    rep10 = np.percentile(rep_corr_ls, 10)

    if plotEnabled:
        sns.set(font_scale=0.7)
        fig, axes = plt.subplots(figsize=(5, 4))
        #         hist_bins=int(len(null_corr_ls)/10)
        #         print(hist_bins)
        sns.distplot(
            null_corr_ls,
            kde=True,
            hist=True,
            bins=hist_bins,
            label="random pairs",
            ax=axes,
            norm_hist=True,
        )
        sns.distplot(
            rep_corr_ls,
            kde=True,
            hist=True,
            bins=hist_bins,
            label="replicate pairs",
            ax=axes,
            norm_hist=True,
            color="r",
        )
        #         perc5=np.percentile(repCC, 50);axes.axvline(x=perc5,linestyle=':',color='darkorange');
        axes.axvline(x=perc95, linestyle=":", label="Null 90th percentile")
        axes.axvline(
            x=rep10,
            linestyle=":",
            color="r",
            label="rep corr 10th percentile",
            markersize=2,
        )
        axes.legend(loc=2, fontsize=8)
        axes.set_title(title)
        axes.set_xlim(-1, 1)
        plt.tight_layout()

    repCorrDf["Rand90Perc"] = perc95
    repCorrDf["Rep10Perc"] = rep10
    #     highRepPertbs=repCorrDf[repCorrDf['RepCor']>perc95].index.tolist()
    #     return repCorrDf

    return fig, repCorrDf


#     return null_corr_ls,rep_corr_ls,repCorrDf


def weighted_replicate_enrichement_profiles(
    inDf, df_counts, pertColName, featColNames, plotEnabled, title=""
):

    """
    Calculates replicate correlation versus across purtburtion correlations

    This function takes the input dataframe and output/plot replicate correlations.

    Parameters:
    inDf   (pandas df): input dataframe contains metadata and features
    pertColName  (str): The column based on which we define replicates of a purturbation
    featColNames(list): The list of all columns corresponding to features
    plotEnabled (bool): If True or 1, plots the curves

    Returns:
    repCorrDf   (list):

    """

    #     inDf=df_meanPrepCorr.copy()
    #     pertColName='Metadata_Sample_uniqe'
    #     featColNames=cpFeats_P_reduced

    df = inDf.copy()

    # df[featColNames]=inDf[featColNames].interpolate();
    uniqPert = df[pertColName].unique().tolist()
    repC = []

    repCorrDf = pd.DataFrame(index=uniqPert, columns=["RepCor"])

    repSizeDF = df.groupby([pertColName]).size().reset_index()
    highRepComp = repSizeDF[repSizeDF[0] > 1][pertColName].tolist()

    rep_corr_ls = []
    null_corr_ls = []

    for u in highRepComp:
        # for u in ['GNRHR E90K-DMSO']:

        not_u = highRepComp[:]
        not_u.remove(u)

        scaledvals_df = df.loc[df[pertColName] == u, featColNames].reset_index(
            drop=True
        )
        weights_df = df_counts.loc[
            df_counts[pertColName] == u, featColNames
        ].reset_index(drop=True)

        u_not_u_df_idx = (
            df.loc[
                df[pertColName].isin(
                    [u] + list(np.random.choice(not_u, scaledvals_df.shape[0]))
                )
            ]
            .groupby(pertColName)
            .sample(1)
            .index
        )
        scaledvals_df_rest = df.loc[u_not_u_df_idx, featColNames].reset_index(drop=True)
        weights_df_rest = df_counts.loc[u_not_u_df_idx, featColNames].reset_index(
            drop=True
        )

        repCorr = calculate_wcorr_row_pairs(scaledvals_df, weights_df)
        corr_rest = calculate_wcorr_row_pairs(scaledvals_df_rest, weights_df_rest)

        repCorrDf.loc[u, "RepCor"] = repCorr

        rep_corr_ls.append(repCorr)
        null_corr_ls.append(corr_rest)

    perc95 = np.percentile(null_corr_ls, 90)
    rep10 = np.percentile(rep_corr_ls, 10)

    if plotEnabled:
        sns.set(font_scale=0.7)
        fig, axes = plt.subplots(figsize=(5, 4))
        sns.distplot(
            null_corr_ls,
            kde=True,
            hist=True,
            bins=100,
            label="random pairs",
            ax=axes,
            norm_hist=True,
        )
        sns.distplot(
            rep_corr_ls,
            kde=True,
            hist=True,
            bins=100,
            label="replicate pairs",
            ax=axes,
            norm_hist=True,
            color="r",
        )
        #         perc5=np.percentile(repCC, 50);axes.axvline(x=perc5,linestyle=':',color='darkorange');
        axes.axvline(x=perc95, linestyle=":", label="Null 90th percentile")
        axes.axvline(
            x=rep10,
            linestyle=":",
            color="r",
            label="rep corr 10th percentile",
            markersize=2,
        )
        axes.legend(loc=2, fontsize=8)
        axes.set_title(title)
        axes.set_xlim(-1, 1)
        plt.tight_layout()

    repCorrDf["Rand90Perc"] = perc95
    repCorrDf["Rep10Perc"] = rep10
    #     highRepPertbs=repCorrDf[repCorrDf['RepCor']>perc95].index.tolist()
    #     return repCorrDf
    return [null_corr_ls, rep_corr_ls, repCorrDf]


def calculate_wcorr_row_pairs(scaledvals_df, weights_df):
    n_reps = scaledvals_df.shape[0]

    u_wcorr_ls = []
    combinations = it.combinations(range(n_reps), 2)
    for c in combinations:
        vec0 = pd.Series(data=scaledvals_df.loc[c[0], :].values)
        vec1 = pd.Series(data=scaledvals_df.loc[c[1], :].values)

        weights_vec = pd.Series(data=weights_df.loc[c, :].sum().values)
        #         weights_vec=pd.Series(data=np.ones((vec0.shape[0])))

        #         wcorr=matthews_corrcoef(vec0, vec1,sample_weight=weights_vec)
        wcorr = WeightedCorr.WeightedCorr(x=vec0, y=vec1, w=weights_vec)(
            method="pearson"
        )
        u_wcorr_ls.append(wcorr)

    repCorr = np.nanmean(u_wcorr_ls)
    return repCorr


def replicateCorrs(inDf, pertColName, featColNames, plotEnabled):

    """
    Calculates replicate correlation versus across purtburtion correlations

    This function takes the input dataframe and output/plot replicate correlations.

    Parameters:
    inDf   (pandas df): input dataframe contains metadata and features
    pertColName  (str): The column based on which we define replicates of a purturbation
    featColNames(list): The list of all columns corresponding to features
    plotEnabled (bool): If True or 1, plots the curves

    Returns:
    repCorrDf   (list):

    """

    df = inDf.copy()
    df[featColNames] = inDf[featColNames].interpolate()
    uniqPert = df[pertColName].unique().tolist()
    repC = []

    repCorrDf = pd.DataFrame(index=uniqPert, columns=["RepCor"])

    repSizeDF = df.groupby([pertColName]).size().reset_index()
    highRepComp = repSizeDF[repSizeDF[0] > 1][pertColName].tolist()

    for u in highRepComp:
        df1 = df[df[pertColName] == u].drop_duplicates().reset_index(drop=True)

        repCorrPurtbs = df1.loc[:, featColNames].T.corr()
        repCorr = list(
            repCorrPurtbs.values[np.triu_indices(repCorrPurtbs.shape[0], k=1)]
        )
        #         print(repCorr)
        repCorrDf.loc[u, "RepCor"] = np.nanmean(repCorr)

        repC = repC + [np.nanmedian(repCorr)]

    randC_v2 = []

    for i in range(1):
        uniqeSamplesFromEachPurt = inDf.groupby(pertColName)[featColNames].apply(
            lambda s: s.sample(1)
        )
        corrMatAcrossPurtbs = uniqeSamplesFromEachPurt.loc[:, featColNames].T.corr()
        randCorrVals = list(
            corrMatAcrossPurtbs.values[
                np.triu_indices(corrMatAcrossPurtbs.shape[0], k=1)
            ]
        )
    randC_v2 = randC_v2 + randCorrVals

    if 0:
        fig, axes = plt.subplots(figsize=(5, 3))
        sns.kdeplot(randC, bw=0.1, label="random pairs", ax=axes)
        sns.kdeplot(repC, bw=0.1, label="replicate pairs", ax=axes)
        axes.set_xlabel("CC")
        sns.kdeplot(randC_v2, bw=0.1, label="random v2 pairs", ax=axes)
        axes.set_xlabel("CC")
        #         perc5=np.percentile(repCC, 50);axes.axvline(x=perc5,linestyle=':',color='darkorange');
        #         perc95=np.percentile(randCC, 90);axes.axvline(x=perc95,linestyle=':');
        axes.legend()
        # axes.set_title('');
        axes.set_xlim(-1.1, 1.1)

    repC = [repC for repC in repC if str(repC) != "nan"]
    randC_v2 = [randC_v2 for randC_v2 in randC_v2 if str(randC_v2) != "nan"]

    perc95 = np.percentile(randC_v2, 90)
    rep10 = np.percentile(repC, 10)

    if plotEnabled:
        fig, axes = plt.subplots(figsize=(5, 4))
        #         sns.kdeplot(randC_v2, bw=.1, label="random pairs",ax=axes);axes.set_xlabel('CC');
        #         sns.kdeplot(repC, bw=.1, label="replicate pairs",ax=axes,color='r');axes.set_xlabel('CC');
        sns.distplot(
            randC_v2,
            kde=True,
            hist=True,
            bins=100,
            label="random pairs",
            ax=axes,
            norm_hist=True,
        )
        sns.distplot(
            repC,
            kde=True,
            hist=True,
            bins=100,
            label="replicate pairs",
            ax=axes,
            norm_hist=True,
            color="r",
        )

        #         perc5=np.percentile(repCC, 50);axes.axvline(x=perc5,linestyle=':',color='darkorange');
        axes.axvline(x=perc95, linestyle=":")
        axes.axvline(x=0, linestyle=":")
        axes.legend(loc=2)
        # axes.set_title('');
        axes.set_xlim(-1, 1)
        plt.tight_layout()

    repCorrDf["Rand90Perc"] = perc95
    repCorrDf["Rep10Perc"] = rep10
    #     highRepPertbs=repCorrDf[repCorrDf['RepCor']>perc95].index.tolist()
    #     return repCorrDf
    return [randC_v2, repC, repCorrDf]


# input is a list of dfs--> [cp,l1k,cp_cca,l1k_cca]
#######
def plotRepCorrs(allData, pertName):
    corrAll = []
    for d in range(len(allData)):
        df = allData[d][0]
        features = allData[d][1]
        uniqPert = df[pertName].unique().tolist()
        repC = []
        randC = []
        for u in uniqPert:
            df1 = df[df[pertName] == u].drop_duplicates().reset_index(drop=True)
            df2 = df[df[pertName] != u].drop_duplicates().reset_index(drop=True)
            repCorr = np.sort(np.unique(df1.loc[:, features].T.corr().values))[
                :-1
            ].tolist()
            #             print(repCorr)
            repC = repC + repCorr
            randAllels = (
                df2[pertName]
                .drop_duplicates()
                .sample(df1.shape[0], replace=True)
                .tolist()
            )
            df3 = pd.concat(
                [
                    df2[df2[pertName] == i].reset_index(drop=True).iloc[0:1, :]
                    for i in randAllels
                ],
                ignore_index=True,
            )
            randCorr = df1.corrwith(df3, axis=1, method="pearson").values.tolist()
            randC = randC + randCorr

        corrAll.append([randC, repC])
    return corrAll


#####################################################################
def replicateCorrs2(inDf, pertColName, featColNames, plotEnabled):

    """
    Calculates replicate correlation versus across purtburtion correlations

    This function takes the input dataframe and output/plot replicate correlations.

    Parameters:
    inDf   (pandas df): input dataframe contains metadata and features
    pertColName  (str): The column based on which we define replicates of a purturbation
    featColNames(list): The column based on which we define replicates of a purturbation
    plotEnabled (bool): If True or 1, plots the curves

    Returns:
    int: Description of return value

    """
    df = inDf.copy()
    uniqPert = df[pertColName].unique().tolist()
    repC = []
    randC = []
    df_repCor = pd.DataFrame(
        columns=[pertColName, "R1", "R2", "Plate1", "Plate2", "cc"]
    )
    for u in uniqPert:
        df1 = df[df[pertColName] == u].drop_duplicates().reset_index(drop=True)
        df2 = df[df[pertColName] != u].drop_duplicates().reset_index(drop=True)

        repCorrPurtbs = df1.loc[:, featColNames].T.corr()
        triuIndices = np.triu_indices(repCorrPurtbs.shape[0], k=1)
        repCorr = list(repCorrPurtbs.values[triuIndices])

        #         repCorr=np.sort(np.unique(df1.loc[:,featColNames].T.corr().values))[:-1].tolist()
        #         repC=repC+repCorr
        #         print(np.median(repCorr))
        repC = repC + [np.median(repCorr)]
        #         randPertbs=df2[pertColName].drop_duplicates().sample(df1.shape[0],replace=True).tolist()
        nS = np.min([len(df2[pertColName].unique().tolist()), df1.shape[0]])
        #         nS=df1.shape[0]

        #         print(nS,[len(df2[pertColName].unique().tolist()),df1.shape[0]])
        for pa in range(len(repCorr)):
            i1 = triuIndices[0][pa]
            i2 = triuIndices[1][pa]
            Plate1 = df1.loc[i1 : i1 + 1, "Metadata_Plate"].tolist()[0]
            Plate2 = df1.loc[i2 : i2 + 1, "Metadata_Plate"].tolist()[0]
            Rep1 = df1.loc[i1 : i1 + 1, "rep"].tolist()[0]
            Rep2 = df1.loc[i2 : i2 + 1, "rep"].tolist()[0]

            dfTemp = pd.DataFrame(
                {
                    pertColName: [u],
                    "R1": [Rep1],
                    "R2": [Rep2],
                    "Plate1": [Plate1],
                    "Plate2": [Plate2],
                    "cc": [repCorr[pa]],
                }
            )
            df_repCor = df_repCor.append(dfTemp, ignore_index=True)

        ##################### correlation with random

        randPertbs = sample(df2[pertColName].unique().tolist(), k=nS)
        #         print(randPertbs)
        df3 = pd.concat(
            [df2[df2[pertColName] == i].sample(1, replace=True) for i in randPertbs],
            ignore_index=True,
        )
        #         print(df1.sample(df3.shape[0],replace=False).shape,df3.shape)
        randCorr = (
            df1[featColNames]
            .sample(df3.shape[0], replace=False)
            .reset_index(drop=True)
            .corrwith(df3[featColNames], axis=1, method="pearson", drop=True)
            .values.tolist()
        )

        #         x1=df1.sample(df3.shape[0],replace=False).values

        #         randCorr=pearsonr()
        #         randCorr = [x for x in randCorr if str(x) != 'nan']
        randC = randC + randCorr
    #     print(randC)
    randC_v2 = []
    for i in range(30):
        uniqeSamplesFromEachPurt = inDf.groupby(pertColName)[featColNames].apply(
            lambda s: s.sample(1)
        )
        corrMatAcrossPurtbs = uniqeSamplesFromEachPurt.loc[:, featColNames].T.corr()
        randCorrVals = list(
            corrMatAcrossPurtbs.values[
                np.triu_indices(corrMatAcrossPurtbs.shape[0], k=1)
            ]
        )
    randC_v2 = randC_v2 + randCorrVals

    if plotEnabled:
        fig, axes = plt.subplots(figsize=(6, 3))
        sns.kdeplot(randC, bw=0.1, label="random pairs", ax=axes)
        sns.kdeplot(repC, bw=0.1, label="replicate pairs", ax=axes)
        axes.set_xlabel("CC")
        sns.kdeplot(randC_v2, bw=0.1, label="random v2 pairs", ax=axes)
        axes.set_xlabel("CC")
        #         perc5=np.percentile(repCC, 50);axes.axvline(x=perc5,linestyle=':',color='darkorange');
        #         perc95=np.percentile(randCC, 90);axes.axvline(x=perc95,linestyle=':');
        axes.legend()
        # axes.set_title('');
    return randC, repC, randC_v2, df_repCor


def categoricalRepCor(inDf, pertColName, featColNames, plotEnabled):

    """
    Calculates replicate correlation versus across purtburtion correlations

    This function takes the input dataframe and output/plot replicate correlations.

    Parameters:
    inDf   (pandas df): input dataframe contains metadata and features
    pertColName  (str): The column based on which we define replicates of a purturbation
    featColNames(list): The column based on which we define replicates of a purturbation
    plotEnabled (bool): If True or 1, plots the curves

    Returns:
    list: [repC,randC]

    """
    df = inDf.copy()
    uniqPert = df[pertColName].unique().tolist()
    repC = []
    #     df_repCor= pd.DataFrame(columns=[pertColName,'R1','R2','Plate1','Plate2','cc'])
    for u in uniqPert:
        df1 = df[df[pertColName] == u].drop_duplicates().reset_index(drop=True)
        #         df2=df[df[pertColName]!=u].drop_duplicates().reset_index(drop=True)

        repCorrPurtbs = df1.loc[:, featColNames].T.corr()
        triuIndices = np.triu_indices(repCorrPurtbs.shape[0], k=1)
        repCorr = list(repCorrPurtbs.values[triuIndices])

        #         repCorr=np.sort(np.unique(df1.loc[:,featColNames].T.corr().values))[:-1].tolist()
        #         repC=repC+repCorr
        #         print(np.median(repCorr))
        repC = repC + [np.median(repCorr)]

    ##################### correlation with random
    #     print(randC)
    randC_v2 = []
    for i in range(10):
        uniqeSamplesFromEachPurt = inDf.groupby(pertColName)[featColNames].apply(
            lambda s: s.sample(1)
        )
        corrMatAcrossPurtbs = uniqeSamplesFromEachPurt.loc[:, featColNames].T.corr()
        randCorrVals = list(
            corrMatAcrossPurtbs.values[
                np.triu_indices(corrMatAcrossPurtbs.shape[0], k=1)
            ]
        )
    randC_v2 = randC_v2 + randCorrVals

    return repC, randC_v2
