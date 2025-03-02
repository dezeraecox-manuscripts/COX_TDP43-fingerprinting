from scipy import stats
import statsmodels as sm
import pandas as pd
import numpy as np

def distribution_ks(dataframe, group_cols, cat_col, cat1, cat2, parameter_col, num_samples=1000, seed=None):
    # Set seed for reproducible random sampling by pandas
    ks_stats = []
    for group, df in dataframe.groupby(group_cols):
        ks, pval = stats.ks_2samp(
        df[df[cat_col] == cat1][parameter_col].sample(
            num_samples, random_state=seed).values,
        df[df[cat_col] == cat2][parameter_col].sample(num_samples, random_state=seed).values
    )
        ks_stats.append([*group, parameter_col, f'{cat1}_{cat2}', ks, pval])

    ks_summary = pd.DataFrame(ks_stats, columns=group_cols+['parameter', 'comparison', 'ks', 'pval'])
    ks_summary['corr_pval'] = sm.stats.multitest.multipletests(
    ks_summary['pval'].tolist(), alpha=0.05, method='bonferroni', is_sorted=False, returnsorted=False)[1]

    ks_summary['sign'] = [
    '***' if pval < 0.001 else (
        '**' if pval < 0.01 else (
            '*' if pval < 0.05 else (
                'ns'
            )
        )
    )
    for pval in ks_summary['corr_pval']]
        
    return ks_summary

def fit_ecdf(x):
    """Function to fit ecdfs for cumulative distributions

    Args:
        x (array): Array to be sorted for fitting the ecdf

    Returns:
        Sorted array
    """
    x = np.sort(x)

    def result(v):
        return np.searchsorted(x, v, side='right') / x.size
    return result


def sample_ecdf(df, value_cols, num_points=100, method='nearest', order=False):
    """Function to interpolate cumulative distributions

    Args:
        df (dataframe): dataframe to compute cumulative distribution
        value_cols (list): Column names to fit ecdf on
        num_points (int, optional): Percentage. Defaults to 100.
        method (str, optional): Fitting method. Defaults to 'nearest'.
        order (bool, optional): Defaults to False.

    Returns:
        dataframe
    """

    test_vals = pd.DataFrame(
        np.arange(0, 1.01, (1/num_points)), columns=['ecdf'])
    test_vals['type'] = 'interpolated'
    interpolated = test_vals.copy()
    for col in value_cols:
        test_df = df.dropna(subset=[col])
        ecdf = fit_ecdf(test_df[col])
        test_df['ecdf'] = ecdf(
            test_df.dropna(subset=[col])[col])
        combined = pd.concat([test_df.sort_values(
            'ecdf').dropna(subset=[col]), test_vals])
        combined = combined.set_index('ecdf').interpolate(
            method=method, order=order).reset_index()
        interpolated[col] = combined[combined['type'] == 'interpolated'].copy()[
            col].tolist()

    return interpolated


def fitting_ecdfs(dataframe, group_cols, val_col):
    """Fitting ecdf for cumulative distrubutions with multiple replicates
    """
    fitted_ecdfs = []
    for group, df in dataframe.groupby(group_cols):
        fitted_ecdf = sample_ecdf(df, value_cols=[
            val_col], method='nearest', order=False)
        fitted_ecdf[group_cols] = group
        fitted_ecdfs.append(fitted_ecdf)

    fitted_ecdfs = pd.concat(fitted_ecdfs)
    return fitted_ecdfs
