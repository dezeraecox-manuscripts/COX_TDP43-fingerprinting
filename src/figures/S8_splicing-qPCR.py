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

# =================Panel A-B: Expression=================
input_path = 'data/orthogonals/PCR_delta-summary.csv'

# Read in summary data
delta_data = pd.read_csv(input_path)
delta_data.drop([col for col in delta_data.columns.tolist() if 'Unnamed: ' in col], axis=1, inplace=True)

# =================Panel C-F: Splicing=================
input_path = 'data/orthogonals/PCR_uni-normed.csv'

# Read in summary data
splice_data = pd.read_csv(input_path)
splice_data.drop([col for col in splice_data.columns.tolist() if 'Unnamed: ' in col], axis=1, inplace=True)

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
}

fig = plt.figure(figsize=(12.53*cm, 6.1*2*cm))
gs = fig.add_gridspec(nrows=2, ncols=3, wspace=0.5, hspace=0.35)

axA = fig.add_subplot(gs[0:1, 0:1])
axC = fig.add_subplot(gs[0:1, 1:2])
axE = fig.add_subplot(gs[0:1, 2:3])
axB = fig.add_subplot(gs[1:2, 0:1])
axD = fig.add_subplot(gs[1:2, 1:2])
axF = fig.add_subplot(gs[1:2, 2:3])

for ax, label, xpos in zip([axA, axB, axC, axD, axE, axF], ['A', 'B', 'C', 'D', 'E', 'F'], [-30, -30, -30, -30, -30, -30]):
# label physical distance to the left and up:
    trans = mtransforms.ScaledTranslation(xpos/72, -8/72, fig.dpi_scale_trans)
    ax.text(0.0, 1.1, label, transform=ax.transAxes + trans,
            fontsize=12, va='bottom', fontweight='bold')

palette = {
    'Control': '#04A777',
    'SOD1': '#820263',
    'TDP-43': '#FB8B24',
}

# =================Panel A-B: Expression=================

genes = ['TARDBP', 'SOD1']

axes = [axA, axB]
for i, gene in enumerate(genes):
    data = delta_data[delta_data['target_gene'] == gene].copy()
    bardotplot(
        data, xcol='Category', ycol='norm_mean', order=['Control', 'SOD1', 'TDP-43'], 
        hue=None, hue_order=None, scat_hue=None, scat_hue_order=None, 
        palette=palette, 
        xlabel='', ylabel='Relative expression', 
        pairs=[('Control', 'TDP-43')],
        correction='holm-bonferroni', 
        xticks=None, 
        groups=None, group_label_y=-0.18, group_line_y=-0.05, 
        ax=axes[i], legend='', 
        dot_size=5, cap_size=0.2, cap_width=1
    )

    axes[i].set_title(gene.upper(), fontsize=8)
    axes[i].set_xticks(np.arange(0, 3), labels=['CRL   ', 'SOD', '   TDP'], fontsize=8)
    axes[i].set_ylabel('Relative expression', labelpad=0.5, fontsize=8)
    

# =================Panel C-F: Splicing=================
genes = ['MADD E31+', 'MADD E31-', 'POLDIP E3+', 'POLDIP E3-']

axes = [axC, axD, axE, axF]
for i, gene in enumerate(genes):
    data = splice_data[splice_data['Target'] == gene].copy()
    bardotplot(
        data, xcol='Sample', ycol='/Uni', order=['Control', 'SOD1', 'TDP-43'], 
        hue=None, hue_order=None, scat_hue=None, scat_hue_order=None, 
        palette=palette, 
        xlabel='', ylabel='Relative expression', 
        pairs=[('Control', 'TDP-43')],
        correction='holm-bonferroni', 
        xticks=None, 
        groups=None, group_label_y=-0.18, group_line_y=-0.05, 
        ax=axes[i], legend='', 
        dot_size=5, cap_size=0.2, cap_width=1
    )

    axes[i].set_title(gene.upper(), fontsize=8)
    axes[i].set_xticks(np.arange(0, 3), labels=['CRL   ', 'SOD', '   TDP'], fontsize=8)
    axes[i].set_ylabel('Relative expression', labelpad=0.5, fontsize=8)
    
    
# ------------Figure admin------------
for ax in [axA, axB, axC, axD, axE, axF]:
    ax.spines[['right', 'top']].set_visible(False)
plt.savefig(f'{output_folder}figure_S8.svg', transparent=True, dpi=1000)
plt.savefig(f'{output_folder}figure_S8.png', transparent=True, dpi=1000)
plt.show()
    