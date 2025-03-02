import os, re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from loguru import logger

logger.info('Import OK')

input_folder = 'results/proteomics/initial_cleanup/'
sample_map = 'data/proteomics/sample_map.xlsx'
output_folder = 'results/proteomics/initial_cleanup/'

if not os.path.exists(output_folder):
    os.makedirs(output_folder)
 
# Read in sample map
sample_map = pd.read_excel(sample_map)

# Read in dataset
peptides = []
proteins = []
file_list = [filename for filename in os.listdir(input_folder) if filename.split('_')[1] in ['Peptides.csv', 'Proteins.csv']]
for filename in file_list:
    data_source = filename.split('_')[0]
    data = pd.read_csv(f'{input_folder}{filename}')
    data.drop([col for col in data.columns.tolist() if 'Unnamed: ' in col], axis=1, inplace=True)
    data['sample'] = data_source
    data.columns = [col.replace(f' {data_source}', '') for col in data.columns.tolist()]
    
    if 'Peptides' in filename:
        peptides.append(data)
    elif 'Proteins' in filename:
        proteins.append(data)    
peptides = pd.concat(peptides)
proteins = pd.concat(proteins)
for data in [peptides, proteins]:
    data[['sample_type', 'replicate']] = data['sample'].str.split('', expand=True)[[1, 2]]
    data['sample_type'] = data['sample_type'].map({'B': 'Beads', 'D': 'Depleted', 'T':'Total',})

peptides.to_csv(f'{output_folder}compiled_peptides.csv')
proteins.to_csv(f'{output_folder}compiled_proteins.csv')

# Melt df and add sample information
summary_proteins = pd.melt(
    proteins,
    id_vars=['Majority protein IDs', 'Gene names', 'Protein names', 'Number of proteins', 'Peptides', 'Razor + unique peptides', 'Unique peptides', 'Intensity', 'sample', 'sample_type', 'replicate'],
    value_vars=[col for col in proteins.columns if 'Reporter intensity corrected' in col],
    value_name='Reporter intensity corrected',
    var_name='channel'
    )
summary_proteins['channel'] = summary_proteins['channel'].str.replace('Reporter intensity corrected ', '')

summary_proteins['tissue_id'] = [f'TIS0{10+int(tis_id)}' if tis_id != '11' else 'reference' for tis_id in summary_proteins['channel']]

summary_proteins = pd.merge(summary_proteins, sample_map.drop('channel', axis=1), how='left', left_on=['tissue_id'], right_on=['Sample_ID'])

# ---------------------Reference normalisation---------------------
# Calculate normalisation factors as the mean sum of the reference channel of each plex, either applied to sample types (tnorm), all samples (anorm), or across all individual channels (cnorm)
norm_factors = summary_proteins[summary_proteins['channel'] == '11'].groupby(['sample_type', 'sample']).sum()['Reporter intensity corrected'].reset_index()
norm_factors['type_mean'] = norm_factors['sample_type'].map(dict(norm_factors.groupby('sample_type').mean().reset_index().values))
norm_factors['all_mean'] = norm_factors['Reporter intensity corrected'].mean()
norm_factors['tnorm_factor'] = norm_factors['type_mean'] / norm_factors['Reporter intensity corrected']
norm_factors['anorm_factor'] = norm_factors['all_mean'] / norm_factors['Reporter intensity corrected']
tnorm_factors = dict(norm_factors[['sample', 'tnorm_factor']].values)
anorm_factors = dict(norm_factors[['sample', 'anorm_factor']].values)

norm_factors = summary_proteins.groupby(['sample_type', 'sample', 'channel']).sum()['Reporter intensity corrected'].reset_index()
norm_factors['channel_mean'] = norm_factors['Reporter intensity corrected'].mean()
norm_factors['cnorm_factor'] = norm_factors['channel_mean'] / \
    norm_factors['Reporter intensity corrected']
norm_factors['key'] = [f'{s}_{c}' for s, c in norm_factors[['sample', 'channel']].values]
cnorm_factors = dict(norm_factors[['key', 'cnorm_factor']].values)

# Add normalisation factors to dataset
summary_proteins['tnorm_factor'] = summary_proteins['sample'].map(tnorm_factors)
summary_proteins['tnorm_corr-intensity'] = summary_proteins['Reporter intensity corrected'] * summary_proteins['tnorm_factor']

summary_proteins['anorm_factor'] = summary_proteins['sample'].map(anorm_factors)
summary_proteins['anorm_corr-intensity'] = summary_proteins['Reporter intensity corrected'] * summary_proteins['anorm_factor']

summary_proteins['cnorm_factor'] = [cnorm_factors[f'{s}_{c}'] for s, c in summary_proteins[['sample', 'channel']].values]
summary_proteins['cnorm_corr-intensity'] = summary_proteins['Reporter intensity corrected'] * summary_proteins['cnorm_factor']

# check normalisation prforms as expected - channel 11 should be equal in either each sample_type (tnorm) or across all sampes (anorm)
test = summary_proteins.groupby(['sample', 'channel']).sum().reset_index()
test[test['channel'] == '11']

# ----------------Save summary to excel----------------
summary_proteins.to_csv(f'{output_folder}normalised_proteins.csv')

# ------------Visualise normalisation process------------
for sample, df in summary_proteins.groupby('sample_type'):
    fig, axes = plt.subplots(3, 1, figsize=(22, 16))
    for i, norm in enumerate(['tnorm', 'anorm', 'cnorm']):
        sns.boxplot(
            data=df, 
            x='channel',
            y=f'{norm}_corr-intensity',
            hue='replicate',
            ax=axes[i],
            fliersize=0,
            palette='magma'
        )
        
        axes[i].set_ylim(0, 500000)
        if i == 0:
            axes[i].legend(bbox_to_anchor=(1.0, 1.0), frameon=False)
        else:
            axes[i].legend('', frameon=False)
    axes[0].set_title(sample)
    plt.show()
