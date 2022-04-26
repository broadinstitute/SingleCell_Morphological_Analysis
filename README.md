# SingleCell_Morphological_Analysis
Basic analysis and visualization of CellProfiler single cell profiles.


Reading Single Cell profiles into the memory¶
All the information about single cells are stored in a sqlite file for each plate
sqlite files are huge (up to 50 GB) and loading them to memory may cause memory errors
Here are alternative ways of handling this issue:


## Reading the whole or a subset of an input .sqlite file
- All the information about single cells are stored in a sqlite file for each plate
sqlite files are huge (up to 50 GB) and loading them to memory may cause memory errors. 
We can load subset of the huge per plate tables to handle the memory issue.


The following table shows different available functions at this repo for reading a subset of .sql files:

| Function               | Subset | 
| :-------------------- | :------------------------------- |
| readSingleCellData_sqlalch     | Reading All the Single Cells of a plate             |
| readSingleCellData_sqlalch_random_image_subset  | Reading random images or defind subset of the plate images             |
| readSingleCellData_sqlalch_well_subset | Reading a subset of wells from the plate            |
| readSingleCellData_sqlalch_features_subset | Reading a subset of features from the plate                | 
| readSingleCellData_sqlalch_FeatureAndWell_subset         | Reading a subset of features and a subset of wells of a plate            |
| readSingleCellData_sqlalch_wellAndObject_subset         | Reading a subset of objects from a subset of wells plate           |


## Visualize single cell examples for both arrayed and pooled data:

  - **cell_selection_methods**
     - **random** - generate n randomly selected cells
     - **representative** - clusters the data and sample from the "closest to mean cluster"
     - **geometric_median** - plots single sample than is the geometric median of samples
        

## Subpopulation analysis


   * (Based on Rohban's elife, and modified to adapt this dataset limitations)
#### A. 1. Morphological Features Group level Feature importance Based on Mean profiles
What groups of morphological features are distinguishing in the Mutant well relative to the WT well?

* [generateFeatureGroupsMap_medianZscore]
   - Z-scored MT feature vectors with respect to the WT samples are calculated first, then averaged to have one vector for MT and then np.median(abs(zScored_averageAcrossSamples)) will be the feature importance score
* [generateFeatureGroupsMap_l2ofDiff]
   - L2 norm of the absolute difference between treatment level of WT and MT 

#### A. 2. Individual Morphological Feature importance Map Based on Mean profiles
Which individual morphological features are distinguishing in the Mut well profile relative to the WT well profile? Red/Blue means the feature has a positive/negative relative difference. Size is proportional to the difference absolute value.
Red/Blue--->negative-positive


#### Subpopulation analysis:
* preprocessing:
   - For each WT- MT pair we have:
      - X transfected cell for WT 
      - X transfected cell for Mutant
      - (If X is more than 2000, it will be randomly subsampled to 2000 single cells)
   - We concatenate single cells for both wells and apply standard scalerization

* How we choose number of clusters:
   - If total WT+MT pairs are more than 60 ---> nClus=20
   - If total number is less than 60 → nClus=⅓ (WT+MT single cells)

#### B. 1. WT+ MUT Single Cells subpopulation enrichment/suppression in the cluster
- Bar plot shows which of 20 subpopulations of cells are enriched and suppressed relative to WT single cells. 
- Plots in the subsequent pages are representing subpopulations whose occurrence differs from WT significantly, whether enriched or suppressed (subpopulations which are very small in both WT/Mut clusters are omitted). 

#### B. 2. Each Cluster Behaviour 
- For each subpopulation, a bar plot shows the top 10 most-distinguishing feature names for single cell mutant wrt all WT single cells 
- Sample images for the cluster are shown
