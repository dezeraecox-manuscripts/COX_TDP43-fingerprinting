import os
import pandas as pd
import matplotlib.pyplot as plt
from src.utilities.visualise import plot_lda_scatter

from loguru import logger
logger.info('Import OK')

output_folder = 'figures/'
if not os.path.exists(output_folder):
    os.makedirs(output_folder)
       
# ================Panel A: homogenate ================

sample_type = 'homogenate'
input_path = f'data/microscopy-tissue/LDA_homogenate-fitted.csv'
homog_lda = pd.read_csv(input_path)

input_path = f'data/microscopy-tissue/LDA_complete-fitted.csv'
all_lda = pd.read_csv(input_path)

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
}
    
region_map = {
    'Cerebellum': 'Cerebellum',
    'FCx': 'Frontal cortex',
    'Fcx': 'Frontal cortex',
}


fig = plt.figure(figsize=(12.5*cm, 6.1*1*cm))

gs = fig.add_gridspec(nrows=1, ncols=2, wspace=0.5, hspace=0.5)
axA = fig.add_subplot(gs[0:1])
axB = fig.add_subplot(gs[1:2])


for ax, label, xpos in zip([axA, axB], ['A', 'B'], [-30, -30, -30]):
# label physical distance to the left and up:
    trans = mtransforms.ScaledTranslation(xpos/72, -3/72, fig.dpi_scale_trans)
    ax.text(0.0, 1.1, label, transform=ax.transAxes + trans,
            fontsize=12, va='bottom', fontweight='bold')    

region_map = {
    'Cerebellum': 'Cerebellum',
    'FCx': 'Frontal cortex',
    'Fcx': 'Frontal cortex',
}

# ------------LDA------------
for ax, data in zip([axA, axB], [homog_lda, all_lda]):
    plot_lda_scatter(data, hue='Category', style='Region', palette=palette, ax=ax, s=50)

    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))

    legend1 = ax.legend([by_label[entry] for entry in ['Control', 'SOD1', 'TDP-43']], ['CRL', 'SOD', 'TDP'], ncol=3, bbox_to_anchor=(0.5, 1.12), frameon=False, columnspacing=0.2, loc='center', handletextpad=0.0025)
    legend2 = ax.legend([by_label[entry] for entry in ['Cerebellum', 'FCx']], ['Cerebellum', 'Frontal cortex'], ncol=3, bbox_to_anchor=(0.5, 1.05), frameon=False, columnspacing=0.2, loc='center', handletextpad=0.0025)
    ax.add_artist(legend1)
    ax.add_artist(legend2)

# ------------Figure admin------------
for ax in [axA, axB]:
    ax.spines[['right', 'top']].set_visible(False)
plt.savefig(f'{output_folder}figure_S5.svg', transparent=True, dpi=1000, bbox_inches="tight")
plt.savefig(f'{output_folder}figure_S5.png', transparent=True, dpi=1000, bbox_inches="tight")
plt.show()
