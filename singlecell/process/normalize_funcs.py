"""
@author: mhaghigh
"""
import pandas as pd
import sklearn.preprocessing as spr


# def standardize_per_catX(df, column_name, cp_features):

#     # column_name='Metadata_Plate'
#     #     cp_features=df.columns[df.columns.str.contains("Cells_|Cytoplasm_|Nuclei_")]
#     df_scaled_perPlate = df.copy()
#     df_scaled_perPlate[cp_features] = (
#         df[cp_features + [column_name]]
#         .groupby(column_name)
#         .transform(lambda x: (x - x.mean()) / x.std())
#         .values
#     )
#     return df_scaled_perPlate


def standardize_per_catX(df, column_name, cp_features):
    # """

    df_scaled_perPlate = df.copy()
    group_means = df.groupby(column_name)[cp_features].mean()
    group_stds = df.groupby(column_name)[cp_features].std()
    df_scaled_perPlate[cp_features] = (
        df[cp_features] - group_means.loc[df[column_name]].values
    ) / group_stds.loc[df[column_name]].values
    return df_scaled_perPlate


def standardize_df_columns(sc_df, cp_features_analysis, scaling_method):
    """
    scaling method can be any method supported by sklearn package, ex. Robust, Standard, ..

    """

    #     scaling_string = "scaler = sp."+scaling_method+"Scaler()"
    #     exec(scaling_string)
    #     scaler = sp.RobustScaler()
    #     scaler = sp.StandardScaler()

    if scaling_method == "MinMax":
        scaler = spr.MinMaxScaler(feature_range=(-1, 1))

    elif scaling_method == "Robust":
        scaler = spr.RobustScaler()

    elif scaling_method == "Standard":
        scaler = spr.StandardScaler()

    sc_df_output = sc_df.copy()
    sc_df_output[cp_features_analysis] = scaler.fit_transform(
        sc_df.loc[:, cp_features_analysis].values
    )

    return sc_df_output


def zscore_df_columns_by_control(
    sc_control_df, sc_df, cp_features_analysis, scaling_method
):
    """
    scaling method can be any method supported by sklearn package, ex. Robust, Standard, ..

    """

    #     scaling_string = "scaler = sp."+scaling_method+"Scaler()"
    #     exec(scaling_string)

    if scaling_method == "MinMax":
        scaler = spr.MinMaxScaler(feature_range=(0, 1))

    elif scaling_method == "Robust":
        scaler = spr.RobustScaler()

    elif scaling_method == "Standard":
        scaler = spr.StandardScaler()

    sc_df_output = sc_df.copy()
    scaler.fit(sc_control_df.loc[:, cp_features_analysis])
    sc_df_output.loc[:, cp_features_analysis] = scaler.transform(
        sc_df.loc[:, cp_features_analysis]
    )

    return sc_df_output


# zscore_df_columns_by_control_perPlate(sc_df,cp_features_analysis,scaling_method, ['X','negcon']):


def zscore_df_columns_by_control_perPlate(
    sc_df, cp_features_analysis, scaling_method, scale_per_col_indic, negcon_col_val
):
    """
    scaling method can be any method supported by sklearn package, ex. Robust, Standard, ..
    negcon_col: [column_name, key_value for negcon]
    """

    #     scaling_string = "scaler = sp."+scaling_method+"Scaler()"
    #     exec(scaling_string)

    if scaling_method == "MinMax":
        scaler = spr.MinMaxScaler(feature_range=(0, 1))

    elif scaling_method == "Robust":
        scaler = spr.RobustScaler()

    elif scaling_method == "Standard":
        scaler = spr.StandardScaler()

    sc_df_output = sc_df.copy()
    plates = sc_df[scale_per_col_indic].unique().tolist()
    for p in plates:
        sc_df_perts = sc_df_output.loc[
            (sc_df_output[scale_per_col_indic] == p)
            & (sc_df_output[negcon_col_val[0]] != negcon_col_val[1]),
            :,
        ]
        sc_df_control = sc_df_output.loc[
            (sc_df_output[scale_per_col_indic] == p)
            & (sc_df_output[negcon_col_val[0]] == negcon_col_val[1]),
            :,
        ]

        scaler.fit(sc_df_control.loc[:, cp_features_analysis])
        sc_df_output.loc[
            (sc_df_output[scale_per_col_indic] == p)
            & (sc_df_output[negcon_col_val[0]] != negcon_col_val[1]),
            cp_features_analysis,
        ] = scaler.transform(sc_df_perts.loc[:, cp_features_analysis])

    return sc_df_output
