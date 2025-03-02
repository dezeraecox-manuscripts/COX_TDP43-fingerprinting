import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from src.utilities.visualise import bardotplot

from loguru import logger
logger.info('Import OK')

output_folder = 'figures/'
if not os.path.exists(output_folder):
    os.makedirs(output_folder)
       
# =================Panel B: flTDP spot count=================
input_path = 'data/microscopy-tissue/flTDP_colabelling-counts.csv'

spots = pd.read_csv(input_path)
spots.drop([col for col in spots.columns.tolist() if 'Unnamed: ' in col], axis=1, inplace=True)

# samples
spots = spots[
    (spots['protein_detected'] != 'None') &
    (spots['sample_type'] == 'h')
    ].copy().dropna(subset='Category')

spots = spots.groupby(['Sample_ID',	'sample_type',	'Region',	'Category',	'Subcategory', 'protein_detected']).mean().reset_index()

spots = spots[(spots['Region'] == 'FCx') & (spots['protein_detected'] != 'flTDP')].copy()

# ================Compile figure================

import matplotlib
import matplotlib.transforms as mtransforms
import subprocess
font = {'family' : 'arial',
'weight' : 'normal',
'size'   : 8 }
matplotlib.rc('font', **font)
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['figure.dpi'] = 500
cm = 1/2.54  # centimeters in inches

# Visualise average fluorescence
palette = {
    'Control': '#04A777',
    'SOD1': '#820263',
    'TDP-43': '#FB8B24',
    'Control_hatch': '#025038',
    'SOD1_hatch': '#50013D',
    'TDP-43_hatch': '#B45904',
}
region_map = {
    'Cerebellum': 'Cerebellum',
    'FCx': 'Frontal cortex',
    'Fcx': 'Frontal cortex',
}


fig = plt.figure(figsize=(12.53*cm, 6.1*2*cm))
gs = fig.add_gridspec(nrows=2, ncols=8, wspace=4.0, hspace=0.35)

axA = fig.add_subplot(gs[0:1, 0:2])
axB = fig.add_subplot(gs[0:1, 2:4])
axC = fig.add_subplot(gs[0:1, 4:6])
axD = fig.add_subplot(gs[0:1, 6:8])
axE = fig.add_subplot(gs[1:2, 0:2])
axF = fig.add_subplot(gs[1:2, 2:4])
axG= fig.add_subplot(gs[1:2, 4:6])
axH = fig.add_subplot(gs[1:2, 6:8])

for ax, label, xpos in zip([axA, axB, axC, axD, axE, axF, axG, axH], ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H',], [-30, -30, -30, -30, -30, -30, -30, -30, ]):
# label physical distance to the left and up:
    trans = mtransforms.ScaledTranslation(xpos/72, -8/72, fig.dpi_scale_trans)
    ax.text(0.0, 1.1, label, transform=ax.transAxes + trans,
            fontsize=12, va='bottom', fontweight='bold')

# -------counts-------
# The Holm–Bonferroni method also controls the FWER at \alpha , but with a lower increase of type II error risk than the classical Bonferroni method.
for (ax, (protein, df)) in zip([axA, axB, axC, axD, axE, axF, axG, axH], spots.groupby('protein_detected')):
    bardotplot(
        df, xcol='Category', ycol='norm_spots_count', order=['Control', 'SOD1', 'TDP-43'], 
        hue=None, hue_order=None, scat_hue=None, scat_hue_order=None, 
        palette=palette, 
        xlabel='', ylabel='Number of puncta', 
        pairs=[('Control', 'TDP-43')],
        correction='holm-bonferroni', 
        xticks=None, 
        groups=None, group_label_y=-0.18, group_line_y=-0.05, 
        ax=ax, legend='', 
        dot_size=5, cap_size=0.2, cap_width=1
    )

    ax.set_title(protein.upper(), fontsize=8)
    ax.set_xticks(np.arange(0, 3), labels=['CRL   ', 'SOD', '   TDP'], fontsize=8)
    ax.set_ylabel('Number of puncta', labelpad=-1, fontsize=8)
    
    for category, patch in zip(['Control', 'SOD1', 'TDP-43'], ax.patches):
        patch.set_edgecolor(palette[f'{category}_hatch'])
        patch.set_hatch('///')
        patch.set_linewidth(0)


# ------------Figure admin------------
for ax in [axA, axB, axC, axD, axE, axF, axG, axH]:
    ax.spines[['right', 'top']].set_visible(False)
plt.savefig(f'{output_folder}figure_S7.svg', transparent=True, dpi=1000)
plt.savefig(f'{output_folder}figure_S7.png', transparent=True, dpi=1000)
plt.show()
    