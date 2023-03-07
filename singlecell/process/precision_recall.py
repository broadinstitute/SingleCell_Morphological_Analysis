"""
@author: mhaghigh
"""
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import average_precision_score
from sklearn.metrics import PrecisionRecallDisplay

################################################################################


def precision_recall_curve_multivar(df, score_col, label_cols):
    """
    This function takes a pandas dataframe and plots precision recall curve multiple columns as reference

    Inputs:
    ++ df_p_s   (pandas df) size --> (number of single cells)x(columns):
    input dataframe contains single cells profiles as rows

    ++ n_cells (dtype: int): number of cells to extract from the input dataframe

    ++ cell_selection_method (str):
        - 'random' - generate n randomly selected cells
        - 'representative' - clusters the data and sample from the "closest to mean cluster"
        - 'geometric_median' - plots single sample than is the geometric median of samples
           -> method 1 (hdmedians package): faster but can handle up to 800 rows without memory error
           -> method 2 (skfda package): slower but can handle up to 1500 rows without memory error

    Returns:
    dff (pandas df): sampled input dataframe
    """

    import matplotlib.pyplot as plt
    from sklearn.metrics import PrecisionRecallDisplay

    _, ax = plt.subplots(figsize=(8, 6))
    plt.axhline(1, 0, 0.955, linestyle="--", label="Ideal")
    plt.axvline(1, 0, 0.955, linestyle="--")
    plt.title("Precision-Recall Curve")
    plt.plot()
    for li in label_cols:
        scores = df[score_col].values
        binary_labels = df[li].values

        ap = average_precision_score(binary_labels, scores, average=None)
        precision, recall, thresholds = precision_recall_curve(binary_labels, scores)
        plt = PrecisionRecallDisplay(precision, recall, average_precision=ap)
        plt.plot(ax=ax, name=li)
    #         plt.legend()

    return
