import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from src.utilities.visualise import plot_volcano
from skimage import io

from loguru import logger
logger.info('Import OK')

output_folder = 'figures/'
if not os.path.exists(output_folder):
    os.makedirs(output_folder)
       
# ================Panel A: Bead Volcano ================
input_path = 'data/proteomics/ratio-summary.csv'

ms_ratios = pd.read_csv(input_path)
ms_ratios.drop([col for col in ms_ratios.columns.tolist() if 'Unnamed: ' in col], axis=1, inplace=True)
ms_ratios['-log10_pval'] = ms_ratios['-log10_pval'].replace(np.inf, np.nan).replace(-np.inf, np.nan)

bead_ratios = ms_ratios[
    (ms_ratios['Category'] == 'TDP-43') & 
    (ms_ratios['sample_type'] == 'Beads') & 
    (ms_ratios['norm_type'] == 'cnorm')
    ].copy().dropna(subset=['log2_ratio', '-log10_pval'])

bead_ratios['thresh_cat'] = ['no change' if ((abs(ratio) < 1) | (abs(pval) < 1.3)) else ('up' if ratio > 0.5 else 'down') for ratio, pval in bead_ratios[['log2_ratio', '-log10_pval']].values]

highlights = {
    # 'Q13148': 'TDP-43',
    'P35637': 'FUS',
    'P21926': 'CD9',
    'P0DMV9;P0DMV8' : 'HSP70s', # (HSPA1A, HSPA1B)
    'O60763': 'USO1', # General vesicular transport factor p115
    'P27824': 'CANX', # Calnexin
    'P06748': 'NPM1', # Nucleophosmin
}

bead_ratios['thresh_cat'] = [col if prot not in highlights else 'highlight' for col, prot in bead_ratios[['thresh_cat', 'Majority protein IDs']].values]


# ================Panel B: Total Volcano ================
ms_ratios = ms_ratios[
    (ms_ratios['Category'] == 'TDP-43') & 
    (ms_ratios['sample_type'] == 'Total') & 
    (ms_ratios['norm_type'] == 'cnorm')
    ].copy().dropna(subset=['log2_ratio', '-log10_pval'])

ms_ratios['thresh_cat'] = ['no change' if ((abs(ratio) < 1) | (abs(pval) < 1.3)) else ('up' if ratio > 0.5 else 'down') for ratio, pval in ms_ratios[['log2_ratio', '-log10_pval']].values]

highlights = {
    # 'Q13148': 'TDP-43',
    'P35637': 'FUS',
    'P21926': 'CD9',
    'P0DMV9;P0DMV8' : 'HSP70s', # (HSPA1A, HSPA1B)
    'O60763': 'USO1', # General vesicular transport factor p115
    'P27824': 'CANX', # Calnexin
    'P06748': 'NPM1', # Nucleophosmin
}
ms_ratios['thresh_cat'] = [col if prot not in highlights else 'highlight' for col, prot in ms_ratios[['thresh_cat', 'Majority protein IDs']].values]


# ================Panel C: GO Enrichment ================

input_path = 'data/proteomics/go-enrichment.csv'
go_enrichment = pd.read_csv(input_path)
go_enrichment.drop([col for col in go_enrichment.columns.tolist() if 'Unnamed: ' in col], axis=1, inplace=True)

go_enrichment['log2_enrichment'] = np.log2(go_enrichment['fold_enrichment'] ).replace(-np.inf, 0).replace(np.inf, 0)
go_enrichment['-log10_pval'] = -np.log10(go_enrichment['pValue'])
go_enrichment[go_enrichment['column'] == 'total_TDP-43_Total_up']


# ================Compile figure================

import matplotlib
import matplotlib.transforms as mtransforms
import subprocess
font = {'family' : 'arial',
'weight' : 'normal',
'size'   : 8 }
matplotlib.rc('font', **font)
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['figure.dpi'] = 1000
cm = 1/2.54  # centimeters in inches

# Visualise average fluorescence
palette = {
    'Control': '#04A777',
    'SOD1': '#820263',
    'TDP-43': '#FB8B24',
    
    # Volcano
    'up': '#E36414',
    'down': '#E36414',
    'no change': '#FEC086',
    'highlight': '#291720',
}
    
region_map = {
    'Cerebellum': 'Cerebellum',
    'FCx': 'Frontal cortex',
    'Fcx': 'Frontal cortex',
}


fig = plt.figure(figsize=(18.1*cm, 6.1*3*cm))

gs = fig.add_gridspec(nrows=3, ncols=3, wspace=0.5, hspace=0.5)
axA = fig.add_subplot(gs[0:1, 0:1])
axB = fig.add_subplot(gs[1:2, 0:1])
axC = fig.add_subplot(gs[0:2, 2:3])
axD = fig.add_subplot(gs[2:3, 0:3])
axDi = fig.add_subplot(gs[2:3, 0:3])

region_map = {
    'Cerebellum': 'Cerebellum',
    'FCx': 'Frontal cortex',
    'Fcx': 'Frontal cortex',
}

# -----------Volcanoes-----------

plot_volcano(bead_ratios[bead_ratios['thresh_cat'] == 'no change'], xcol='log2_ratio', ycol='-log10_pval', palette=palette, hue_col='thresh_cat', ax=axA, ythresh=1.3, xthresh=[-1, 1], alpha=0.7)
plot_volcano(bead_ratios[bead_ratios['thresh_cat'].isin(['up', 'down'])], xcol='log2_ratio', ycol='-log10_pval', palette=palette, hue_col='thresh_cat', ax=axA, ythresh=1.3, xthresh=[-1, 1], alpha=0.7)
plot_volcano(bead_ratios[bead_ratios['thresh_cat'] == 'highlight'], xcol='log2_ratio', ycol='-log10_pval', palette=palette, hue_col='thresh_cat', ax=axA, ythresh=1.3, xthresh=[-1, 1])
axA.legend('', frameon=False)
axA.set_xlabel('Log$_2$ Abundance ratio', labelpad=-1)
axA.set(ylabel='-Log$_{10}$ p-value', xlim=(-3.25, 3.25))

# Add annotations for highlighted spots
positions = bead_ratios[bead_ratios['thresh_cat'] == 'highlight'][['Majority protein IDs', 'log2_ratio', '-log10_pval']].copy()
positions = {key: (round(pos1, 3), round(pos2, 3)) for key, pos1, pos2 in positions.values}

label_positions = {
    # 'Q13148': (0, 0), #TDP
    'P35637': (2.0, 4.5), #FUS
    'P21926': (-2, 0.5), #CD9
    'P0DMV9;P0DMV8' : (-1.7, 4.1),  #HSPs
    'O60763' : (2.1, 2.5),  #USO1
    'P27824': (2, 3.2), #CANX
    'P06748': (-1.9, 1.7), # Nucleophosmin
    
}

for label in label_positions:
    label_pos = label_positions[label]
    arrow_pos = positions[label]
    ha = 'left' if label_pos[0] > 0 else 'right'
    offset = 0.1 if label_pos[0] < 0 else -0.1
    axA.annotate(highlights[label], xy=label_positions[label], fontsize=5, ha=ha, va='center', zorder=100000)
    axA.plot([label_pos[0]+offset, arrow_pos[0]], [label_pos[1], arrow_pos[1]], color=palette['highlight'], linewidth=0.5, zorder=100000)
    
    
# -----------Volcanoes-----------

plot_volcano(ms_ratios[ms_ratios['thresh_cat'] == 'no change'], xcol='log2_ratio', ycol='-log10_pval', palette=palette, hue_col='thresh_cat', ax=axB, ythresh=1.3, xthresh=[-1, 1], alpha=0.7)
plot_volcano(ms_ratios[ms_ratios['thresh_cat'].isin(['up', 'down'])], xcol='log2_ratio', ycol='-log10_pval', palette=palette, hue_col='thresh_cat', ax=axB, ythresh=1.3, xthresh=[-1, 1], alpha=0.7)
plot_volcano(ms_ratios[ms_ratios['thresh_cat'] == 'highlight'], xcol='log2_ratio', ycol='-log10_pval', palette=palette, hue_col='thresh_cat', ax=axB, ythresh=1.3, xthresh=[-1, 1])
axB.legend('', frameon=False)
axB.set_xlabel('Log$_2$ Abundance ratio', labelpad=-1)
axB.set(ylabel='-Log$_{10}$ p-value', xlim=(-3.25, 3.25))

# Add annotations for highlighted spots
positions = ms_ratios[ms_ratios['thresh_cat'] == 'highlight'][['Majority protein IDs', 'log2_ratio', '-log10_pval']].copy()
positions = {key: (round(pos1, 3), round(pos2, 3)) for key, pos1, pos2 in positions.values}

label_positions = {
    # 'Q13148': (0, 0), #TDP
    'P35637': (-2.0, 0.95), #FUS
    'P21926': (-2, 0), #CD9
    'P0DMV9;P0DMV8' : (-1.95, 4),  #HSPs
    'O60763' : (-2, 2.75),  #USO1
    'P27824': (1.9, 3.5), #CANX
    'P06748': (1.2, 5), # Nucleophosmin
    
}

for label in label_positions:
    label_pos = label_positions[label]
    arrow_pos = positions[label]
    ha = 'left' if label_pos[0] > 0 else 'right'
    offset = 0.1 if label_pos[0] < 0 else -0.1
    axB.annotate(highlights[label], xy=label_positions[label], fontsize=5, ha=ha, va='center', zorder=100000)
    axB.plot([label_pos[0]+offset, arrow_pos[0]], [label_pos[1], arrow_pos[1]], color=palette['highlight'], linewidth=0.5, zorder=100000)
    
# -----------enrichment-----------

df = go_enrichment[(go_enrichment['column'] == 'total_TDP-43_Total')].copy()
df = df[df['term_label'] != 'UNCLASSIFIED'].copy()
sns.scatterplot(
    data=df.sort_values(['search_type', 'log2_enrichment']),
    y='term_label',
    x='log2_enrichment',
    size='-log10_pval',
    color='#E36414',
    ax=axC
)
axC.set_ylabel('GO term name')
axC.set_xlabel('Fold enrichment')

handles, labels = axC.get_legend_handles_labels()
by_label = dict(zip(reversed(labels), reversed(handles)))
axC.legend([handles[0], handles[-1]], [' ', ' '], frameon=False, loc='upper right', handlelength=1., handletextpad=0.25, ncol=2, title='$-Log_{10}$ p-value',)
axC.plot([0.4, 2.5], [1.6, 1.35], color='lightgrey', linestyle='--', linewidth=0.5, zorder=0)
axC.plot([0.4, 2.5], [2.42, 2.72], color='lightgrey', linestyle='--', linewidth=0.5, zorder=0)

        
# ------------Figure admin------------
for ax in [axA, axB, axC]:
    ax.spines[['right', 'top']].set_visible(False)

for ax, label, (xpos, ypos) in zip([axA, axB, axC, axD, ], ['A', 'B', 'C', 'D', ], [(-30, -3), (-30, -3), (-170, -18), (-30, -3)]):
# label physical distance to the left and up:
    trans = mtransforms.ScaledTranslation(xpos/72, ypos/72, fig.dpi_scale_trans)
    ax.text(0.0, 1.1, label, transform=ax.transAxes + trans,
            fontsize=12, va='bottom', fontweight='bold')    

plt.savefig(f'{output_folder}figure_S6.svg', transparent=True)
plt.savefig(f'{output_folder}figure_S6.png', transparent=True, dpi=1000)
plt.show()
