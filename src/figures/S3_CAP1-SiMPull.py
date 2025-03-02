
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from statannotations.Annotator import Annotator
from src.utilities.visualise import bardotplot
from skimage import io
import subprocess
import matplotlib.patches as mpatches
from statannotations.Annotator import Annotator
from pingouin  import anova, pairwise_tests

from loguru import logger
logger.info('Import OK')

output_folder = 'figures/'

if not os.path.exists(output_folder):
    os.makedirs(output_folder)
    

# ================Panel B: CAP1 pulldown================

input_path = 'data/microscopy-tissue/CAP1_pulldown-spots-summary.csv'

cap_pd = pd.read_csv(input_path)
cap_pd.drop([col for col in cap_pd.columns.tolist() if 'Unnamed: ' in col], axis=1, inplace=True)

# Controls
cap_pd_controls = cap_pd[(cap_pd['detect'] == 'IgG') & (cap_pd['sample'].str.contains('MIX')) &(cap_pd['channel'] == '638')].copy()

# samples
cap_pd = cap_pd[
    (cap_pd['channel'] == 638)
    ].copy().dropna(subset='Category')

cap_pd = cap_pd[['Sample_ID', 'sample_type', 'Region',	'Category',	'Subcategory', 'spots_count']].groupby(['Sample_ID', 'sample_type', 'Region',	'Category',	'Subcategory']).mean().reset_index()

# Perform ANOVA
anova_res = anova(data=cap_pd, dv='spots_count', between=['Category', 'sample_type'])
posthoc_res = pairwise_tests(data=cap_pd, dv='spots_count', between=['Category', 'sample_type'], padjust='holm')
posthoc_res = posthoc_res[(posthoc_res['Category'] != '-') & (posthoc_res['Category'] != 'SOD1')][['Category', 'p-corr']].copy()
posthoc_res['symbol'] = ['***' if pval < 0.001 else ('**' if pval < 0.01 else ('*' if pval < 0.05 else ('ns'))) for pval in posthoc_res['p-corr']]

# =================Panel D-E: CAP1 spot count=================
input_path = 'data/microscopy-tissue/CAP1_normalised-spot-counts.csv'

cap_count = pd.read_csv(input_path)
cap_count.drop([col for col in cap_count.columns.tolist() if 'Unnamed: ' in col], axis=1, inplace=True)

# Controls
cap_controls = cap_count[(cap_count['detect'] == 'IgG') & (cap_count['sample'].str.contains('MIX')) & (cap_count['sample_type'] == 'homogenate')].copy()

# samples
cap_count = cap_count[
    (cap_count['channel'] == 638)
    ].copy().dropna(subset='Category')

cap_count = cap_count[['Sample_ID', 'sample_type', 'Region',	'Category',	'Subcategory', 'protein_detected', 'spots_count_cmbnorm']].groupby(['Sample_ID', 'sample_type', 'Region',	'Category',	'Subcategory', 'protein_detected']).mean().reset_index()

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
    
    'total': 'black',
    'depleted': 'lightgrey',
}

region_map = {
    'Cerebellum': 'Cerebellum',
    'FCx': 'Frontal cortex',
    'Fcx': 'Frontal cortex',
}


fig = plt.figure(figsize=(18.8*cm, 6.1*2*cm))

gs = fig.add_gridspec(nrows=2, ncols=5, wspace=0.55, hspace=0.35)

axC = fig.add_subplot(gs[1:2, 0:1])
axCi = fig.add_subplot(gs[1:2, 0:1])
axDi = fig.add_subplot(gs[1:2, 1:2])
axDii = fig.add_subplot(gs[1:2, 2:3])
axEi = fig.add_subplot(gs[1:2, 3:4])
axEii = fig.add_subplot(gs[1:2, 4:5])

axA = fig.add_subplot(gs[0:1, 0:3])
axAi = fig.add_subplot(gs[0:1, 0:3])
axB = fig.add_subplot(gs[0:1, 3:5])

for ax, label, pos in zip([axA, axB, axC, axDi, axEi, ], ['A', 'B', 'C', 'D', 'E'], [-25, -30, -25, -30, -30]):
# label physical distance to the left and up:
    trans = mtransforms.ScaledTranslation(pos/72, -8/72, fig.dpi_scale_trans)
    ax.text(0.0, 1.1, label, transform=ax.transAxes + trans,
            fontsize=12, va='bottom', fontweight='bold')

# -------CAP1 pulldown schematic-------
axA.axis('off')
axAi.axis('off')
axAi.set_position([0.075, 0.5, 0.48, 0.5])

# -------CAP1 pulldown counts-------
pairs = [
    (('Control', 'total'), ('Control', 'depleted')),
    (('TDP-43', 'total'), ('TDP-43', 'depleted')),
]
data=cap_pd
xcol='Category'
ycol='spots_count'
order=['Control', 'SOD1', 'TDP-43']
hue='sample_type'
hue_order=['total', 'depleted']
scat_hue='sample_type'
scat_hue_order=['total', 'depleted']
correction='holm-bonferroni'
dot_size=5
cap_size=0.2
cap_width=1

sns.barplot(
    data=data,
    x=xcol,
    y=ycol,
    hue=hue,
    palette=palette,
    capsize=cap_size,
    errwidth=cap_width,
    ax=axB,
    dodge=True,
    order=order,
    hue_order=hue_order,
    edgecolor='white'
)
for cat in order:
    sns.stripplot(
        data=data[data['Category'] == cat],
        x=xcol,
        y=ycol,
        hue=scat_hue,
        palette={'total': palette[cat], 'depleted': palette[cat]},
        ax=axB,
        edgecolor='#fff',
        linewidth=1,
        s=dot_size,
        order=order,
        hue_order=scat_hue_order,
        dodge=True,
    )

annotator = Annotator(
    ax=axB, pairs=pairs, data=data, x=xcol, y=ycol, order=order, hue=hue, hue_order=hue_order)
annotator.configure(test=None, loc='inside', line_width=0.5)
annotator.set_custom_annotations(posthoc_res['symbol'].tolist())
annotator.annotate()

for cat, stype, patch in zip(['Control', 'SOD1', 'TDP-43']*2, ['total', 'total', 'total', 'depleted', 'depleted', 'depleted'], axB.patches):
    if stype == 'depleted':
        patch.set_facecolor('#fff')
        patch.set_edgecolor(palette[cat])
        patch.set_linewidth(1)
    else:
        patch.set_facecolor(palette[f'{cat}'])
        patch.set_edgecolor(palette[f'{cat}_hatch'])
        patch.set_hatch('///')
        patch.set_linewidth(0)
d_patch = mpatches.Patch(facecolor='white', label='Depleted', edgecolor='black', linewidth=1)
t_patch = mpatches.Patch(facecolor='black', label='Total', edgecolor='black', linewidth=1)
axB.legend(handles=[t_patch, d_patch], frameon=False, bbox_to_anchor=(0.45, 1.1))
        
axB.axhline(cap_pd_controls['spots_count'].mean(), linestyle='--', color='lightgrey')
axB.set_ylabel('Number of puncta', labelpad=-0.4, fontsize=8)
axB.set_xlabel('')
axB.set_xticks(np.arange(0, 3), labels=['CRL', 'SOD', 'TDP'], fontsize=8)

# -------CAP1 SiMPull schematic-------
axC.axis('off')
axCi.axis('off')
axCi.set_position([0.075, 0.075, 0.15, 0.4])

# -------CAP1 counts-------
for sample_type, axes in zip(['soak', 'homogenate'], [(axDi, axDii), (axEi, axEii)]):
    dataframe = cap_count[cap_count['sample_type'] == sample_type].copy()
    for (ax, (region, df)) in zip(axes, dataframe.groupby('Region')):
        bardotplot(
            df, xcol='Category', ycol='spots_count_cmbnorm', order=['Control', 'SOD1', 'TDP-43'], 
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

        if sample_type == 'homogenate':
            for category, patch in zip(['Control', 'SOD1', 'TDP-43'], ax.patches):
                patch.set_edgecolor(palette[f'{category}_hatch'])
                patch.set_hatch('///')
                patch.set_linewidth(0)
                
    # ------------Figure admin------------
        ax.axhline(cap_controls['spots_count'].mean(), linestyle='--', color='lightgrey')
        ax.set_title(region_map[region.capitalize()], fontsize=8)
        ax.set_xticks(np.arange(0, 3), labels=['CRL   ', 'SOD', '   TDP'], fontsize=8)
    axes[0].set_ylabel('Number of puncta', labelpad=-0.2, fontsize=8)
    axes[1].set_ylabel('', labelpad=-0.2)
    
for ax in [axB, axDi, axDii, axEi, axEii]:
    ax.spines[['right', 'top']].set_visible(False)
plt.tight_layout()
plt.savefig(f'{output_folder}figure_S3.svg', transparent=True)
plt.savefig(f'{output_folder}figure_S3.png', transparent=True)
plt.show()
    
    