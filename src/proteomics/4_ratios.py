import os, re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import ttest_1samp
from collections import defaultdict

from loguru import logger

logger.info('Import OK')

input_path = 'results/proteomics/initial_cleanup/normalised_proteins.csv'
output_folder = 'results/proteomics/ratios/'

if not os.path.exists(output_folder):
    os.makedirs(output_folder)
    
    # ---------Visualisation---------
palette = {
    'Control': '#04A777',
    'TDP-43': '#FB8B24',
    'SOD1': '#820263',

}
# Visualise distriubution of average ratios per protein per cohort


def plot_distribution(ratios, col, palette):
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    for i, (cat, df) in enumerate(ratios.groupby(['Category'])):
        sns.distplot(np.log2(df.groupby(['Majority protein IDs']).mean()[col]), color=palette[cat], ax=axes[i])
        axes[i].set_title(f'{cat}')
    plt.show()


# Visualise volcano plot
def plot_volcano(ttest_results, xcol, ycol, palette, xlim=(-2.1, 2.1), ylim=(-0.01, 4.1)):
    fig, axes = plt.subplots(1, len(ttest_results['Category'].unique()), figsize=(6*len(ttest_results['Category'].unique()), 5))
    for i, (cat, df) in enumerate(ttest_results.groupby(['Category'])):

        sns.scatterplot(data=df, x=xcol,
                        y=ycol, color=palette[cat], ax=axes[i])
        axes[i].set_ylim(*ylim)
        axes[i].set_xlim(*xlim)
        axes[i].set_title(cat)
        axes[i].axvline(1, linestyle='--', color='lightgrey')
        axes[i].axvline(-1, linestyle='--', color='lightgrey')
        axes[i].axhline(1.3, linestyle='--', color='lightgrey')
    plt.show()



# PERCEPT as defined in Eq. 1
def percept(m0, m1, F, p):
    return m0 + ((m0 - m1) / -(F**p))
    

# define a simple function to apply the PERCEPT scaling
def apply_percept(data, hypothethical_mean, penalty):
    # 1. Calculate p-value    
    tval, pval = ttest_1samp(data, popmean=hypothethical_mean, nan_policy='omit')
    # 2. Calculate sample mean
    sample_mean = np.mean(data)    
    # 3. Apply percept, returning scaled mean value
    return percept(m0=hypothethical_mean, m1=sample_mean, F=penalty, p=pval)


# Read in dataset
raw_data = pd.read_csv(input_path)
raw_data.drop([col for col in raw_data.columns.tolist() if 'Unnamed: ' in col], axis=1, inplace=True)

# Calculate ratios for each normalisation type for each sample type
ratios_combined = []
for sample_type, data in raw_data.groupby('sample_type'):

    for norm in ['tnorm', 'anorm', 'cnorm']:
        # Calculate ratio to control samples for TDP-43 and SOD1 independently - this is done using the average of the control replicates

        controls = data[data['Category'] == 'Control'].copy()
        controls = controls.groupby(['Majority protein IDs', 'replicate']).mean()[[f'{norm}_corr-intensity']].copy()
        controls.rename(columns={f'{norm}_corr-intensity':f'{norm}_control'}, inplace=True)

        ratios = data[data['tissue_id'] != 'reference'].copy()
        ratios = pd.merge(
            ratios[['Category', 'sample_type', 'sample', 'replicate', 'tissue_id', 'Majority protein IDs', 'Gene names', 'Protein names', 'Number of proteins', 'Peptides', 'Razor + unique peptides', 'Unique peptides', 'Intensity', f'{norm}_corr-intensity']],
            controls.reset_index(),
            on=['Majority protein IDs','replicate'],
            how='left'
        )
        # Remove proteins for which no control mean can be assigned
        ratios.dropna(subset=[f'{norm}_control'], inplace=True, how='any') # loose ~1000 entries

        ratios[f'ratio'] = ratios[f'{norm}_corr-intensity'] / ratios[f'{norm}_control']

        # Calculate mean of technical replicates for each tissue sample
        ratios = ratios.groupby(['Category', 'sample_type', 'tissue_id', 'Majority protein IDs', 'Gene names', 'Protein names']).mean()[['ratio']].copy().reset_index()
        ratios['norm_type'] = norm
        ratios_combined.append(ratios)
ratios_combined = pd.concat(ratios_combined)

# Perform t-test
ratio_summary = []
# Perform ttest comparing ratios in each Category to hypothetical value of 1
for (norm, category, protein, sample_type), df in ratios_combined.groupby(['norm_type', 'Category', 'Majority protein IDs', 'sample_type']):
    vals = df['ratio'].tolist()
    tval, pval = ttest_1samp(df['ratio'].tolist(), popmean=1, nan_policy='omit')
    # Calculate mean per cohort (per protein), then log2 for plotting
    df = df.groupby(['norm_type', 'Category', 'sample_type', 'Majority protein IDs', 'Gene names', 'Protein names']).mean().reset_index().rename(columns={'ratio': 'mean_ratio'})
    df['ratio_vals'] = [vals]
    df['tval'] = tval
    df['pval'] = pval
    
    # Apply PERCEPT
    df['percept'] = apply_percept(vals, hypothethical_mean=1, penalty=50)
    
    ratio_summary.append(df)
ratio_summary = pd.concat(ratio_summary)
ratio_summary[f'log2_ratio'] = np.log2(ratio_summary['mean_ratio'])
ratio_summary[f'-log10_pval'] = - np.log10(ratio_summary['pval'])

# ---------Visualise ratios of interest---------
# For total ratios, want to consider the tnorm ratio as we assume that all samples started with an equivalent total protein concentration

plot_volcano(ratio_summary[(ratio_summary['sample_type'] == 'Total') & (ratio_summary['norm_type'] == 'tnorm')], xcol='log2_ratio', ycol='-log10_pval', palette=palette, xlim=(-2.1, 2.1), ylim=(-0.01, 4.1))

plot_volcano(ratio_summary[(ratio_summary['sample_type'] == 'Beads') & (ratio_summary['norm_type'] == 'tnorm')], xcol='log2_ratio', ycol='-log10_pval', palette=palette, xlim=(-2.1, 2.1), ylim=(-0.01, 4.1))


poi = ratio_summary[
    (ratio_summary['sample_type'] == 'Beads') &
    (ratio_summary['Category'] == 'TDP-43') &
    (ratio_summary['norm_type'] == 'tnorm') &
    (ratio_summary['-log10_pval'] >= 1.3) &
    (ratio_summary['log2_ratio'].abs() >= 1)
    ].copy().sort_values(['mean_ratio'])
poi.head(50)
poi.tail(48)


ratio_summary[(ratio_summary['sample_type'] == 'Total') & (ratio_summary['norm_type'] == 'tnorm') &(ratio_summary['Gene names'].str.contains('TARDBP'))]

# ------------Ratio of ratios to normalise Bead/Depleted ratios------------

# Calculate ratio-of-ratios to normalise changes in beads/depleted samples according to existing differences in the total proteome
normalised_ratios = []
for (norm, sample_type, cat), data in ratios_combined.groupby(['norm_type', 'sample_type', 'Category']):
    if sample_type == 'Total':
        continue
    if cat == 'Control':
        continue
    totals = ratios_combined[(ratios_combined['norm_type'] == norm) & (ratios_combined['Category'] == cat) & (ratios_combined['sample_type'] == 'Total') ][['tissue_id', 'Majority protein IDs', 'ratio']].copy()
    totals['key'] = [f'{tid}_{pro}' for tid, pro in totals[['tissue_id', 'Majority protein IDs']].values]
    total_ratios = defaultdict(lambda: np.nan)
    total_ratios.update(dict(totals[['key', 'ratio']].values))
    
    data['total_ratio'] = [total_ratios[f'{tid}_{pro}'] for tid, pro in data[['tissue_id', 'Majority protein IDs']].values]
    
    data['ror_ratio'] = data['ratio'] / data['total_ratio'] #a fair few proteins are potentialy lost at this point due to not being quantified in the total samples (~below abundance)

    # Perform ttest comparing ratios in each Category to hypothetical value of 1
    ttest_results = []
    for (protein), df in data.groupby(['Majority protein IDs']):
        vals = df['ror_ratio'].tolist()
        tval, pval = ttest_1samp(df['ror_ratio'].tolist(), popmean=1, nan_policy='omit')
        # Calculate mean per cohort (per protein), then log2 for plotting
        df = df.groupby(['Category', 'sample_type', 'Majority protein IDs', 'Gene names', 'Protein names']).mean().reset_index().rename(columns={'ror_ratio': 'mean_ror_ratio'})
        df['ratio_vals'] = [vals]
        df['tval'] = tval
        df['pval'] = pval
        ttest_results.append(df)

        # Apply PERCEPT
        df['percept'] = apply_percept(vals, hypothethical_mean=1, penalty=50)

    ttest_results = pd.concat(ttest_results)
    ttest_results['norm_type'] = norm
    ttest_results[f'log2_ratio'] = np.log2(ttest_results['mean_ror_ratio'])
    ttest_results[f'-log10_pval'] = - np.log10(ttest_results['pval'])    
    
    normalised_ratios.append(ttest_results)
normalised_ratios = pd.concat(normalised_ratios)


plot_volcano(normalised_ratios[(normalised_ratios['sample_type'] == 'Depleted') & (normalised_ratios['norm_type'] == 'cnorm')], xcol='log2_ratio', ycol='-log10_pval', palette=palette, xlim=(-3.1, 3.1), ylim=(-0.01, 6.1))

norm_poi = normalised_ratios[
    (normalised_ratios['sample_type'] == 'Depleted') &
    (normalised_ratios['Category'] == 'TDP-43') &
    (normalised_ratios['norm_type'] == 'cnorm') &
    (normalised_ratios['-log10_pval'] >= 1.3) &
    (normalised_ratios['log2_ratio'].abs() >= 1)
    ]['Majority protein IDs'].tolist()
poi.head(50)
poi.tail(15)

# -----------Save dataframes-----------
ratios_combined.to_csv(f'{output_folder}ratios_combined.csv')
ratio_summary.to_csv(f'{output_folder}ratio_summary.csv')
normalised_ratios.to_csv(f'{output_folder}normalised_ratios.csv')

# Save dataframe for STRING
normalised_ratios  = pd.read_csv(f'{output_folder}normalised_ratios.csv')
normalised_ratios.drop([col for col in normalised_ratios.columns.tolist() if 'Unnamed: ' in col], axis=1, inplace=True)
normalised_ratios['-log10_pval'] = normalised_ratios['-log10_pval'].replace(np.inf, np.nan).replace(-np.inf, np.nan)

ms_ratios = normalised_ratios[
    (normalised_ratios['Category'] == 'TDP-43') & 
    (normalised_ratios['sample_type'] == 'Beads') & 
    (normalised_ratios['norm_type'] == 'cnorm')
    ].copy().dropna(subset=['log2_ratio', '-log10_pval'])

ms_ratios['thresh_cat'] = ['no change' if ((abs(ratio) < 1) | (abs(pval) < 1.3)) else ('up' if ratio > 0.5 else 'down') for ratio, pval in ms_ratios[['log2_ratio', '-log10_pval']].values]

ms_ratios[ms_ratios['thresh_cat'] == 'up'].to_csv(f'{output_folder}Beads_TDP_UP-for_STRING.csv')