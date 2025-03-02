import os
import pandas as pd
import numpy as np

from loguru import logger
logger.info('Import OK')

input_beads = 'results/proteomics/ratios/normalised_ratios.csv'
input_total = 'results/proteomics/ratios/ratio_summary.csv'
background_path = 'results/proteomics/initial_cleanup/compiled_proteins.csv'
output_folder = 'results/proteomics/go_enrichment/'

if not os.path.exists(output_folder):
    os.mkdir(output_folder)


def apply_enrichment(df, searches=None, obo_path='PANTHERGOslim.obo', organism='10090', refOrganism='10090', enrichmentTestType='FISHER', correction='BONFERONNI', min_proteins=5, reference=None):
    """
    Worker function to apply GO enrichment for df of proteins against the reference list (background) on a per-column basis. Column labels are returned in the 'column' key
    """

    if not searches:
        searches = {'Bio Process': 'ANNOT_TYPE_ID_PANTHER_GO_SLIM_BP', 'Mol Function': 'ANNOT_TYPE_ID_PANTHER_GO_SLIM_MF', 'Cell Component': 'ANNOT_TYPE_ID_PANTHER_GO_SLIM_CC'}

    enrichment_dict = {}
    for col in df.columns.tolist():
        if len(df[col].dropna()) >= min_proteins:
            for abbr, search_type in searches.items():
                enrichment_dict[(col, abbr)] = enrichment(organism=organism, refOrganism=refOrganism, annotDataSet=search_type, enrichmentTestType=enrichmentTestType, correction=correction, genes=df[col].dropna(), obo_path=obo_path, reference=reference)

    summary = []
    for (col, search_type), result in enrichment_dict.items():
        result['column'] = col
        result['search_type'] = search_type
        summary.append(result)
    summary = pd.concat(summary)

    return summary


# -----------------Read in standard components-----------------
protein_lists = {}
for label, path in zip(['total', 'beads', 'background'], [input_total, input_beads, background_path]):
    data = pd.read_csv(f'{path}')
    data.drop([col for col in data.columns.tolist()
                if 'Unnamed: ' in col], axis=1, inplace=True)
    data['Proteins'] = data['Majority protein IDs'].str.split(';')
    if label == 'total':
        data = data[data['norm_type'] == 'tnorm'].copy()
    elif label == 'beads':
        data = data[data['norm_type'] == 'cnorm'].copy()

    if label != 'background':
        data = data[(data['log2_ratio'].abs() >=1 ) & (data['-log10_pval'] >= 1.3 )].copy()        
        
        for (cat, sample_type), df in data.groupby(['Category', 'sample_type']):
            proteins = df.explode('Proteins')['Proteins'].unique().tolist()
            protein_lists[f'{label}_{cat}_{sample_type}'] = proteins
            
            proteins_up = df[(df['log2_ratio'] >=1 )].explode('Proteins')['Proteins'].unique().tolist()  
            protein_lists[f'{label}_{cat}_{sample_type}_up'] = proteins_up
            
            proteins_down = df[(df['log2_ratio'] <= -1 )].explode('Proteins')['Proteins'].unique().tolist()  
            protein_lists[f'{label}_{cat}_{sample_type}_down'] = proteins_down
            
        else:    
            background = data.explode('Proteins')['Proteins'].unique().tolist()

# ----------------Perform Panther enrichment test----------------
combinations = pd.DataFrame(
    dict([(k, pd.Series(v)) for k, v in protein_lists.items()]))
combinations.to_csv(f'{output_folder}regulated_combinations.csv')

# perform enrichment test
enrichment = apply_enrichment(combinations, searches=None, obo_path='resources/PANTHERGOslim.obo', organism='9606', refOrganism='9606', enrichmentTestType='FISHER', correction='FDR', min_proteins=2, reference=background)

# Save all to excel
enrichment.to_csv(f'{output_folder}go_enrichment.csv')
