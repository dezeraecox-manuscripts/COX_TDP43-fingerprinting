import os, re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from loguru import logger

logger.info('Import OK')

input_path = f'results/proteomics/txt/proteinGroups.txt'
sample_map = 'data/proteomics/labels_peptides.xlsx'
output_folder = 'results/proteomics/summary/'

if not os.path.exists(output_folder):
    os.makedirs(output_folder)

# Read in sample map
sample_map = pd.read_excel(sample_map)
sample_map = sample_map[['ID', 'tissue']].copy()

# Read in datasets
raw_data = pd.read_table(input_path)
raw_data.drop([col for col in raw_data.columns.tolist() if 'Unnamed: ' in col], axis=1, inplace=True)

# Collect useful columns
columns = [
    'Protein IDs',
    'Protein names',
    'Gene names',
    'Number of proteins',
    'Sequence coverage [%]',
    'Mol. weight [kDa]',
    'Sequence length',
    'Intensity B1',
    'Intensity B2',
    'Intensity B3',
    'Intensity B4',
    'Intensity T1',
    'Intensity T2',
    'Intensity T3',
    'Intensity T4',
    'Intensity D1',
    'Intensity D2',
    'Intensity D3',
    'Intensity D4',
    'Reverse',
    'Potential contaminant',
    ]


raw_data = raw_data[columns].copy()

# Remove reverse, contaminant proteins
cleaned = raw_data[(raw_data['Reverse'] != '+') & (raw_data['Potential contaminant'] != '+')].copy()
cleaned.drop(['Reverse', 'Potential contaminant'], axis=1, inplace=True)
[columns.remove(val) for val in ['Reverse', 'Potential contaminant']]

# Melt intensity columns 
cleaned = pd.melt(
    cleaned,
    id_vars=[col for col in columns if 'Intensity' not in col],
    value_vars=[col for col in columns if 'Intensity' in col],
    value_name='intensity',
    var_name='sample'
)
cleaned['intensity'] = cleaned['intensity'].replace(0.0, np.nan)
cleaned['sample'] = cleaned['sample'].str.replace('Intensity ', '')

# Add sample info

# Save summary
cleaned.to_csv(f'{output_folder}id_summary.csv')

