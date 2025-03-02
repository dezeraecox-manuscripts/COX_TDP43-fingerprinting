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
input_path = 'data/orthogonals/ASC_summary.csv'

# Read in summary data
asc_speck = pd.read_csv(input_path)
asc_speck.drop([col for col in asc_speck.columns.tolist() if 'Unnamed: ' in col], axis=1, inplace=True)

asc_speck = asc_speck[['Region', 'Category', 'Sample_ID', 'sample_type','AgregateCounts', 'AgregateArea', 'AgregateBrightness']].groupby(['Region', 'Category', 'Sample_ID', 'sample_type']).mean().reset_index()

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

fig = plt.figure(figsize=(12.53*cm, 6.1*cm))
gs = fig.add_gridspec(nrows=1, ncols=4, wspace=0.6, hspace=0.35)

axAi = fig.add_subplot(gs[0:1, 0:1])
axAii = fig.add_subplot(gs[0:1, 1:2])
axBi = fig.add_subplot(gs[0:1, 2:3])
axBii = fig.add_subplot(gs[0:1, 3:4])

for ax, label, xpos in zip([axAi, axBi], ['A', 'B'], [-30, -25]):
# label physical distance to the left and up:
    trans = mtransforms.ScaledTranslation(xpos/72, -8/72, fig.dpi_scale_trans)
    ax.text(0.0, 1.1, label, transform=ax.transAxes + trans,
            fontsize=12, va='bottom', fontweight='bold')

palette = {
    'Control': '#04A777',
    'SOD1': '#820263',
    'TDP-43': '#FB8B24',
    'Control_hatch': '#025038',
    'SOD1_hatch': '#50013D',
    'TDP-43_hatch': '#B45904',
}
reg_map = {
    'Cerebellum': 'Cerebellum',
    'FCx': 'Frontal cortex'
}

axes = [axAi, axAii, axBi, axBii]

for i, ((sample_type, region), df) in enumerate(asc_speck.groupby(['sample_type', 'Region'])):
    print(f'{sample_type} {region}')
    
    bardotplot(
        data=df, 
        xcol='Category', 
        ycol='AgregateCounts', 
        order=['Control', 'SOD1', 'TDP-43'], 
        hue=None, 
        hue_order=None, 
        scat_hue=None, 
        scat_hue_order=None, 
        palette=palette, 
        xlabel='', 
        ylabel='Number of puncta', 
        pairs=[('Control', 'TDP-43')], 
        correction='holm-bonferroni', 
        xticks=None, 
        groups=None, 
        group_label_y=-0.18, 
        group_line_y=-0.05, 
        ax=axes[i], 
        legend='', 
        dot_size=5, 
        cap_size=0.2, 
        cap_width=1)
    axes[i].set_title(reg_map[region], fontsize=8)
    axes[i].set_xticks(np.arange(0, 3), labels=['CRL   ', 'SOD', '   TDP'], fontsize=8)
    axes[i].set_ylabel('Number of puncta', labelpad=-1, fontsize=8)
    
    if sample_type == 'homog':
        for category, patch in zip(['Control', 'SOD1', 'TDP-43'], axes[i].patches):
            patch.set_edgecolor(palette[f'{category}_hatch'])
            patch.set_hatch('///')
            patch.set_linewidth(0)
    
# ------------Figure admin------------
for ax in [axAi, axAii, axBi, axBii]:
    ax.spines[['right', 'top']].set_visible(False)
plt.savefig(f'{output_folder}figure_S10.svg', transparent=True, dpi=1000)
plt.savefig(f'{output_folder}figure_S10.png', transparent=True, dpi=1000)
plt.show()
    