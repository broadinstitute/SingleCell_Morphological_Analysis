import numpy as np
from sklearn import datasets
from scipy.stats import f


def TwoSampleT2Test(X, Y):
    """
    reference: https://www.r-bloggers.com/2020/10/hotellings-t2-in-julia-python-and-r/
    """
    nx, p = X.shape
    ny, _ = Y.shape
    delta = np.mean(X, axis=0) - np.mean(Y, axis=0)
    Sx = np.cov(X, rowvar=False)
    Sy = np.cov(Y, rowvar=False)
    S_pooled = ((nx - 1) * Sx + (ny - 1) * Sy) / (nx + ny - 2)
    S_pooled = S_pooled + np.eye(S_pooled.shape[0]) * 1e-6
    t_squared = (
        (nx * ny)
        / (nx + ny)
        * np.matmul(np.matmul(delta.transpose(), np.linalg.inv(S_pooled)), delta)
    )
    statistic = t_squared * (nx + ny - p - 1) / (p * (nx + ny - 2))
    F = f(p, nx + ny - p - 1)
    p_value = 1 - F.cdf(statistic)
    #     print(f"Test statistic: {statistic}\nDegrees of freedom: {p} and {nx+ny-p-1}\np-value: {p_value}")
    return statistic, p_value


# from hotelling.stats import hotelling_t2


# def TwoSampleT2Test_hotelling(X, Y):
#     dif_ctl_t, _, dif_ctl_p, _ = hotelling.stats.hotelling_t2(
#         per_site_df_pert_plate[uncorr_feats_condese], control_df[uncorr_feats_condese]
#     )

#     return dif_ctl_t, dif_ctl_p
