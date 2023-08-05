"""
This script contains usefull functions used in the notebooks

@author: mhaghigh
"""
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import os

# from functools import reduce
from sklearn.cluster import KMeans, SpectralClustering
from . import visualize_n_SingleCell


def add_clustering_index_column(
    sc_df, cp_features_analysis, clustering_method, n_clusters
):
    if clustering_method == "kmeans":
        cls_model = KMeans(n_clusters=n_clusters)  # , random_state=0)

    elif clustering_method == "spectral":
        cls_model = SpectralClustering(
            n_clusters=n_clusters, affinity="nearest_neighbors"
        )

    else:
        print("Add clustering method to this function here!")

    cluster_indices = cls_model.fit(sc_df[cp_features_analysis].values).labels_
    distance_to_cluster_cent = cls_model.fit_transform(
        sc_df[cp_features_analysis].values
    )

    sc_df["cluster_label"] = cluster_indices
    sc_df["distance_2cluster_center"] = np.min(distance_to_cluster_cent, 1)

    return sc_df


def subpopulation_clustering_visualization(
    sc_df, cp_features_analysis, params, pooled=False
):
    """
    This function select cells based on input cell_selection_method

    Inputs:
    ++ DirsDict (dict) a dictionary containing all the paths for reading the images and saving the results
       keys:
           - resDir: results directory
    ++ df_sc   (pandas df) input dataframe which contains single cells profiles for the conditions under
                           comparision
    ++ contLabel (str): value for reference perturbation for comparision
    ++ d (str): value for target perturbation for comparision
    ++ nClus (int): number of clusters
    ++ cp_features_analysis (list): list of features to use for clustering
    ++ channels (list): list of channels for single cell visualization
    ++ boxSize (int): single cell box size  for visualizations

    """

    resultsDir = params["results_dir"]
    saveFormat = params["img_save_format"]
    contLabel = params["control_label"]
    d = params["pert_label"]
    nClus = params["n_clusters"]
    nSampleSCs = params["n_cells_4single_cell_plots"]
    channels = params["channels_4single_cell_plots"]
    boxSize = params["boxSize_4single_cell_plots"]
    compressed_ims_flag = params["compressed_ims_flag"]

    if compressed_ims_flag:
        compressed_im_size = params["compressed_im_size"]
    else:
        compressed_im_size = None

    fig, axes = plt.subplots(1, 2)

    sc_pert = sc_df[(sc_df["label"] == d)]["cluster_label"].values
    sc_control = sc_df[sc_df["label"] == contLabel]["cluster_label"].values
    # print(sc_pert)

    print("sc shapes pert/control: ", sc_pert.shape, sc_control.shape)

    hist_pert, bin_edges = np.histogram(sc_pert, range(nClus + 1), density=True)
    hist_control, bin_edges = np.histogram(sc_control, range(nClus + 1), density=True)

    histDiff = hist_pert - hist_control

    map_4ordered_clusters = dict(zip(np.argsort(histDiff), range(nClus)))
    # print(map_4ordered_clusters)
    sc_df["clusterLabels2"] = sc_df["cluster_label"].map(map_4ordered_clusters)

    sc_pert = sc_df[~(sc_df["label"] == contLabel)]["clusterLabels2"].values
    #     print(sc_pert)
    sc_control = sc_df[sc_df["label"] == contLabel]["clusterLabels2"].values

    hist_pert, bin_edges = np.histogram(sc_pert, range(nClus + 1), density=True)
    hist_control, bin_edges = np.histogram(sc_control, range(nClus + 1), density=True)
    histDiff = hist_pert - hist_control

    #     print(hist_pert, bin_edges)
    sns.distplot(
        sc_pert,
        kde=False,
        norm_hist=True,
        bins=bin_edges,
        label=d,
        ax=axes[0],
        color="r",
        hist_kws=dict(edgecolor="k"),
    )
    sns.distplot(
        sc_control,
        kde=False,
        norm_hist=True,
        bins=bin_edges,
        label=contLabel,
        ax=axes[0],
        hist_kws=dict(edgecolor="k"),
    )
    sns.distplot(
        sc_pert,
        kde=False,
        hist=True,
        norm_hist=False,
        bins=bin_edges,
        label=d,
        ax=axes[1],
        color="r",
        hist_kws=dict(edgecolor="k"),
    )
    sns.distplot(
        sc_control,
        kde=False,
        hist=True,
        norm_hist=False,
        bins=bin_edges,
        label=contLabel,
        ax=axes[1],
        hist_kws=dict(edgecolor="k"),
    )

    #   axes.xaxis.set_ticklabels(range(0,20,2));
    axes[0].set_ylabel("Density")
    axes[0].set_xlabel("cell category index")
    axes[1].set_ylabel("Histogram")
    axes[1].set_xlabel("cell category index")
    axes[0].legend()
    axes[1].legend()
    plt.tight_layout()

    os.system("mkdir -p " + resultsDir)
    fig.savefig(resultsDir + "/clusterDensity" + saveFormat)

    mean_control = sc_df.loc[sc_df["label"] == contLabel, cp_features_analysis].mean()

    diff_threshold = 1 / (nClus * 100)

    for c in range(len(hist_pert)):
        if abs(histDiff[c]) > diff_threshold and (
            hist_pert[c] > 0.001 or hist_control[c] > 0.001
        ):
            #         c=3;
            #             print(histDiff[c])
            clusterDF = sc_df[sc_df["clusterLabels2"] == c].reset_index(drop=True)
            mean_pert_cluster = clusterDF.loc[
                ~(clusterDF["label"] == contLabel), cp_features_analysis
            ].mean()
            diffOfMutMeanAndWTMean = pd.DataFrame(
                data=mean_pert_cluster.values - mean_control.values,
                columns=["Diff"],
                index=mean_pert_cluster.index,
            )
            diffOfMutMeanAndWTMean.loc[:, "Diff2"] = diffOfMutMeanAndWTMean.loc[
                :, "Diff"
            ].abs()
            absFeatureImportanceSS = diffOfMutMeanAndWTMean.sort_values(
                "Diff2", ascending=False
            )[:10]
            fig, axes = plt.subplots()
            sns.barplot(
                x="Diff",
                y=absFeatureImportanceSS.index,
                data=absFeatureImportanceSS,
                ax=axes,
            )
            sns.despine()
            plt.tight_layout()
            fig.savefig(
                resultsDir + "/cluster" + str(c) + "_barImpFeatures" + saveFormat
            )
            plt.close("all")
            #             nSampleSCs=6
            if clusterDF.shape[0] > nSampleSCs:
                samples2plot = (
                    clusterDF.sort_values("distance_2cluster_center", ascending=True)
                    .sample(nSampleSCs)
                    .reset_index(drop=True)
                )
                title_str = "Cluster " + str(c)
                if pooled:
                    im_size = 5500
                    f = visualize_n_SingleCell.visualize_n_SingleCell_pooled(
                        channels, samples2plot, boxSize, im_size, title=title_str
                    )
                else:
                    f = visualize_n_SingleCell.visualize_n_SingleCell(
                        channels,
                        samples2plot,
                        boxSize,
                        info_columns=["label"],
                        title=title_str,
                        compressed=compressed_ims_flag,
                        compressed_im_size=compressed_im_size,
                    )
                f.savefig(resultsDir + "/cluster" + str(c) + "_examplar" + saveFormat)
                plt.close("all")

    return


#################################
def subpopulation_triplet_clustering_visualization(
    sc_df, cp_features_analysis, params, pooled=False
):
    resultsDir = params["results_dir"]
    saveFormat = params["img_save_format"]
    contLabel = params["control_label"]
    d = params["pert_1_label"]
    tr = params["pert_2_label"]
    nClus = params["n_clusters"]
    nSampleSCs = params["n_cells_4single_cell_plots"]
    channels = params["channels_4single_cell_plots"]
    boxSize = params["boxSize_4single_cell_plots"]
    compressed_ims_flag = params["compressed_ims_flag"]

    if compressed_ims_flag:
        compressed_im_size = params["compressed_im_size"]
    else:
        compressed_im_size = None

    fig, axes = plt.subplots(1, 2)
    # wtANDmtDf['clusterLabels']
    data2plotMut = sc_df[(sc_df["label"] == d)]["cluster_label"].values
    data2plotWT = sc_df[sc_df["label"] == contLabel]["cluster_label"].values
    data2plotTR = sc_df[sc_df["label"] == tr]["cluster_label"].values

    histMut, bin_edges = np.histogram(data2plotMut, range(nClus + 1), density=True)
    histWT, bin_edges = np.histogram(data2plotWT, range(nClus + 1), density=True)
    histTR, bin_edges = np.histogram(data2plotTR, range(nClus + 1), density=True)
    histMaxWTTR = np.max((histWT, histTR), axis=0)

    #     histDiff=histMut-histWT;
    histDiff = histMut - histMaxWTTR
    #     print(histDiff)
    sortedDiff = np.sort(histDiff)
    # ind=[np.where(histDiff[i]==sortedDiff)[0][0] for i in range(len(histDiff))]
    ind = []
    for i in range(len(histDiff)):
        iinndd = np.where(histDiff[i] == sortedDiff)[0].tolist()
        for j in range(len(iinndd)):
            if iinndd[j] not in ind:
                ind = ind + [iinndd[j]]
                break

    sc_df["clusterLabels2"] = sc_df["cluster_label"].replace(range(nClus), ind)
    data2plotMut = sc_df[sc_df["label"] == d]["clusterLabels2"].values
    data2plotWT = sc_df[sc_df["label"] == contLabel]["clusterLabels2"].values
    data2plotTR = sc_df[sc_df["label"] == tr]["clusterLabels2"].values

    histMut, bin_edges = np.histogram(data2plotMut, range(nClus + 1), density=True)
    histWT, bin_edges = np.histogram(data2plotWT, range(nClus + 1), density=True)
    histTR, bin_edges = np.histogram(data2plotTR, range(nClus + 1), density=True)
    histMaxWTTR = np.max((histWT, histTR), axis=0)
    # def mapToNewLabelCats(histMut,histWT):

    sns.distplot(
        data2plotMut,
        kde=False,
        norm_hist=True,
        bins=bin_edges,
        label=d,
        ax=axes[0],
        color="r",
        hist_kws=dict(edgecolor="k"),
    )
    sns.distplot(
        data2plotWT,
        kde=False,
        norm_hist=True,
        bins=bin_edges,
        label=contLabel,
        ax=axes[0],
        hist_kws=dict(edgecolor="k"),
    )
    sns.distplot(
        data2plotTR,
        kde=False,
        norm_hist=True,
        bins=bin_edges,
        label=tr,
        ax=axes[0],
        color="g",
        hist_kws=dict(edgecolor="k"),
    )

    sns.distplot(
        data2plotMut,
        kde=False,
        hist=True,
        norm_hist=False,
        bins=bin_edges,
        label=d,
        ax=axes[1],
        color="r",
        hist_kws=dict(edgecolor="k"),
    )
    sns.distplot(
        data2plotWT,
        kde=False,
        hist=True,
        norm_hist=False,
        bins=bin_edges,
        label=contLabel,
        ax=axes[1],
        hist_kws=dict(edgecolor="k"),
    )
    sns.distplot(
        data2plotTR,
        kde=False,
        hist=True,
        norm_hist=False,
        bins=bin_edges,
        label=tr,
        ax=axes[1],
        color="g",
        hist_kws=dict(edgecolor="k"),
    )

    #   axes.xaxis.set_ticklabels(range(0,20,2));
    axes[0].set_ylabel("Density")
    axes[0].set_xlabel("cell category index")
    axes[1].set_ylabel("Histogram")
    axes[1].set_xlabel("cell category index")
    axes[0].legend()
    axes[1].legend()
    plt.tight_layout()
    os.system("mkdir -p " + resultsDir)
    fig.savefig(resultsDir + "/clusterDensity" + saveFormat)

    #     meanWT=wtANDmtDf_scaled.loc[wtANDmtDf_scaled['label'] == contLabel,feats2use].mean()
    meanMT = sc_df.loc[sc_df["label"] == d, cp_features_analysis].mean()

    for c in range(len(histMut)):
        if histMut[c] > 0.001 or histMaxWTTR[c] > 0.001:
            if abs(histMut[c] - histMaxWTTR[c]) > 0.01:
                #             print(c,histMut[c]-histMaxWTTR[c])
                #         c=3;
                clusterDF = sc_df[sc_df["clusterLabels2"] == c].reset_index(drop=True)
                #             meanMutCluster=clusterDF.loc[~(clusterDF['label'] == contLabel),feats2use].mean();
                meanNonMutCluster = clusterDF.loc[
                    ~(clusterDF["label"] == d), cp_features_analysis
                ].mean()
                #             diffOfMutMeanAndWTMean=pd.DataFrame(data=meanMutCluster.values-meanWT.values,columns=['Diff'],index=meanMutCluster.index);
                diffOfMutMeanAndWTMean = pd.DataFrame(
                    data=meanNonMutCluster.values - meanMT.values,
                    columns=["Diff"],
                    index=meanNonMutCluster.index,
                )
                diffOfMutMeanAndWTMean.loc[:, "Diff2"] = diffOfMutMeanAndWTMean.loc[
                    :, "Diff"
                ].abs()
                absFeatureImportanceSS = diffOfMutMeanAndWTMean.sort_values(
                    "Diff2", ascending=False
                )[:10]
                fig, axes = plt.subplots()
                sns.barplot(
                    x="Diff",
                    y=absFeatureImportanceSS.index,
                    data=absFeatureImportanceSS,
                    ax=axes,
                )
                sns.despine()
                plt.tight_layout()
                fig.savefig(
                    resultsDir + "/cluster" + str(c) + "_barImpFeatures" + saveFormat
                )
                plt.close("all")
                #                 nSampleSCs=6

                if clusterDF.shape[0] > nSampleSCs:
                    samples2plot = (
                        clusterDF.sort_values(
                            "distance_2cluster_center", ascending=True
                        )
                        .sample(nSampleSCs)
                        .reset_index(drop=True)
                    )
                    title_str = "Cluster " + str(c)
                    if pooled:
                        im_size = 5500
                        f = visualize_n_SingleCell.visualize_n_SingleCell_pooled(
                            channels, samples2plot, boxSize, im_size, title=title_str
                        )
                    else:
                        f = visualize_n_SingleCell.visualize_n_SingleCell(
                            channels,
                            samples2plot,
                            boxSize,
                            info_columns=["label"],
                            title=title_str,
                            compressed=compressed_ims_flag,
                            compressed_im_size=compressed_im_size,
                        )
                    f.savefig(
                        resultsDir + "/cluster" + str(c) + "_examplar" + saveFormat
                    )
                    plt.close("all")

    return


# def visualize_n_SingleCell(channels , sc_df , boxSize , outline = False, color=False,\
#                            info_columns=[], title="",compressed=False,compressed_im_size=None):


def clusteringHists(
    DirsDict,
    wtANDmtDf_scaled,
    contLabel,
    d,
    nClus,
    feats2use,
    compartments,
    boxSize,
    pooled=False,
):
    """
    This function select cells based on input cell_selection_method

    Inputs:
    ++ DirsDict (dict) a dictionary containing all the paths for reading the images and saving the results
       keys:
           - resDir: results directory
    ++ df_sc   (pandas df) input dataframe which contains single cells profiles for the conditions under
                           comparision
    ++ contLabel (str): value for reference perturbation for comparision
    ++ d (str): value for target perturbation for comparision
    ++ nClus (int): number of clusters
    ++ feats2use (list): list of features to use for clustering
    ++ compartments (list): list of channels for single cell visualization
    ++ boxSize (int): single cell box size  for visualizations

    """

    #     rootDir=DirsDict['root']
    resultsDir = DirsDict["resDir"]
    #     DirsDict['imDir']=rootDir+"Mito_Morphology_input/images/"

    #     d1=d.split(" ")[0]
    saveFormat = ".png"
    #'.png'
    #     plt.ioff()
    fig, axes = plt.subplots(1, 2)
    # wtANDmtDf['clusterLabels']
    data2plotMut = wtANDmtDf_scaled[(wtANDmtDf_scaled["label"] == d)][
        "clusterLabels"
    ].values
    data2plotWT = wtANDmtDf_scaled[wtANDmtDf_scaled["label"] == contLabel][
        "clusterLabels"
    ].values
    print(data2plotMut.shape, data2plotWT.shape)
    histMut, bin_edges = np.histogram(data2plotMut, range(nClus + 1), density=True)
    histWT, bin_edges = np.histogram(data2plotWT, range(nClus + 1), density=True)

    histDiff = histMut - histWT
    sortedDiff = np.sort(histDiff)
    # ind=[np.where(histDiff[i]==sortedDiff)[0][0] for i in range(len(histDiff))]
    ind = []
    for i in range(len(histDiff)):
        iinndd = np.where(histDiff[i] == sortedDiff)[0].tolist()
        for j in range(len(iinndd)):
            if iinndd[j] not in ind:
                ind = ind + [iinndd[j]]
                break

    wtANDmtDf_scaled["clusterLabels2"] = wtANDmtDf_scaled["clusterLabels"].replace(
        range(nClus), ind
    )
    data2plotMut = wtANDmtDf_scaled[~(wtANDmtDf_scaled["label"] == contLabel)][
        "clusterLabels2"
    ].values
    data2plotWT = wtANDmtDf_scaled[wtANDmtDf_scaled["label"] == contLabel][
        "clusterLabels2"
    ].values

    histMut, bin_edges = np.histogram(data2plotMut, range(nClus + 1), density=True)
    histWT, bin_edges = np.histogram(data2plotWT, range(nClus + 1), density=True)
    # def mapToNewLabelCats(histMut,histWT):

    sns.distplot(
        data2plotMut,
        kde=False,
        norm_hist=True,
        bins=bin_edges,
        label=d,
        ax=axes[0],
        color="r",
        hist_kws=dict(edgecolor="k"),
    )
    sns.distplot(
        data2plotWT,
        kde=False,
        norm_hist=True,
        bins=bin_edges,
        label=contLabel,
        ax=axes[0],
        hist_kws=dict(edgecolor="k"),
    )
    sns.distplot(
        data2plotMut,
        kde=False,
        hist=True,
        norm_hist=False,
        bins=bin_edges,
        label=d,
        ax=axes[1],
        color="r",
        hist_kws=dict(edgecolor="k"),
    )
    sns.distplot(
        data2plotWT,
        kde=False,
        hist=True,
        norm_hist=False,
        bins=bin_edges,
        label=contLabel,
        ax=axes[1],
        hist_kws=dict(edgecolor="k"),
    )

    #   axes.xaxis.set_ticklabels(range(0,20,2));
    axes[0].set_ylabel("Density")
    axes[0].set_xlabel("cell category index")
    axes[1].set_ylabel("Histogram")
    axes[1].set_xlabel("cell category index")
    axes[0].legend()
    axes[1].legend()
    plt.tight_layout()
    os.system("mkdir -p " + resultsDir)
    fig.savefig(resultsDir + "/clusterDensity" + saveFormat)

    meanWT = wtANDmtDf_scaled.loc[
        wtANDmtDf_scaled["label"] == contLabel, feats2use
    ].mean()

    for c in range(len(histMut)):
        if histMut[c] > 0.001 or histWT[c] > 0.001:
            #         c=3;
            clusterDF = wtANDmtDf_scaled[
                wtANDmtDf_scaled["clusterLabels2"] == c
            ].reset_index(drop=True)
            meanMutCluster = clusterDF.loc[
                ~(clusterDF["label"] == contLabel), feats2use
            ].mean()
            diffOfMutMeanAndWTMean = pd.DataFrame(
                data=meanMutCluster.values - meanWT.values,
                columns=["Diff"],
                index=meanMutCluster.index,
            )
            diffOfMutMeanAndWTMean.loc[:, "Diff2"] = diffOfMutMeanAndWTMean.loc[
                :, "Diff"
            ].abs()
            absFeatureImportanceSS = diffOfMutMeanAndWTMean.sort_values(
                "Diff2", ascending=False
            )[:10]
            fig, axes = plt.subplots()
            sns.barplot(
                x="Diff",
                y=absFeatureImportanceSS.index,
                data=absFeatureImportanceSS,
                ax=axes,
            )
            sns.despine()
            plt.tight_layout()
            fig.savefig(
                resultsDir + "/cluster" + str(c) + "_barImpFeatures" + saveFormat
            )
            plt.close("all")
            nSampleSCs = 6
            if clusterDF.shape[0] > nSampleSCs:
                samples2plot = (
                    clusterDF.sort_values("dist2Mean", ascending=True)
                    .sample(nSampleSCs)
                    .reset_index(drop=True)
                )
                title_str = "Cluster " + str(c)
                if pooled:
                    im_size = 5500
                    f = visualize_n_SingleCell.visualize_n_SingleCell_pooled(
                        compartments, samples2plot, boxSize, im_size, title=title_str
                    )
                else:
                    f = visualize_n_SingleCell.visualize_n_SingleCell(
                        compartments, samples2plot, boxSize, title=title_str
                    )
                f.savefig(resultsDir + "/cluster" + str(c) + "_examplar" + saveFormat)
                plt.close("all")

    return


def check_feature_similarity_dendrogram(data, feature_names, figsize):
    from scipy.cluster import hierarchy
    from scipy.spatial.distance import squareform

    #     import hdmedians as hd

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)
    #     corr = spearmanr(X).correlation
    corr = data[feature_names].corr().values

    # Ensure the correlation matrix is symmetric
    corr = (corr + corr.T) / 2
    np.fill_diagonal(corr, 1)

    # We convert the correlation matrix to a distance matrix before performing
    # hierarchical clustering using Ward's linkage.
    distance_matrix = 1 - np.abs(corr)
    dist_linkage = hierarchy.ward(squareform(distance_matrix))
    dendro = hierarchy.dendrogram(
        dist_linkage, labels=feature_names, ax=ax1, leaf_rotation=90
    )
    dendro_idx = np.arange(0, len(dendro["ivl"]))

    pos = ax2.imshow(corr[dendro["leaves"], :][:, dendro["leaves"]], vmin=-1, vmax=1)
    ax2.grid(False)
    fig.colorbar(pos, ax=ax2)
    ax2.set_xticks(dendro_idx)
    ax2.set_yticks(dendro_idx)
    ax2.set_xticklabels(dendro["ivl"], rotation="vertical")
    ax2.set_yticklabels(dendro["ivl"])
    fig.tight_layout()

    return fig
