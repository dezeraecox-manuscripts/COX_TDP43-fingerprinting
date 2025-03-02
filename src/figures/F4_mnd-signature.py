import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from src.utilities.visualise import plot_lda_scatter, plot_radar, radar_factory

from loguru import logger

logger.info('Import OK')

output_folder = 'figures/'

if not os.path.exists(output_folder):
    os.makedirs(output_folder)
       
# ================Panel A-C: Spyder ================
input_path = 'data/microscopy-tissue/LDA_proportion-categories.csv'

# read in raw data
proportions = pd.read_csv(input_path)
proportions.drop([col for col in proportions.columns.tolist() if 'Unnamed: ' in col], axis=1, inplace=True)

categories = ['long', 'fibrillar',  'dense', 'bright', 'coloc', 'spot_count']
proportions = proportions[
    (proportions['threshold_category'].isin(categories))
    ].copy()
proportions['proportion'] = [val / 100 if cat == 'coloc' else val for val, cat in proportions[['proportion', 'threshold_category']].values]

proportions = proportions[
    (proportions['sample_type'] == 'soak')
].copy()

# ================Panel D: LDA ================

sample_type = 'soak'
input_path = f'data/microscopy-tissue/LDA_{sample_type}-fitted.csv'
for_LDA = pd.read_csv(input_path)

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

radar_factory(num_vars=6, frame='polygon')

fig = plt.figure(figsize=(18.8*cm, 6.1*1*cm))

gs = fig.add_gridspec(nrows=1, ncols=3, wspace=0.5, hspace=0)
axA = fig.add_subplot(gs[0:1], projection='radar')
axAi = fig.add_subplot(gs[0:1])
axB = fig.add_subplot(gs[1:2], projection='radar')
axBi = fig.add_subplot(gs[1:2])
axC = fig.add_subplot(gs[2:3])

for ax in [axAi, axBi]:
    ax.axis('off')

for ax, label, xpos in zip([axAi, axC], ['A', 'B'], [-30, -30, -30]):
# label physical distance to the left and up:
    trans = mtransforms.ScaledTranslation(xpos/72, -3/72, fig.dpi_scale_trans)
    ax.text(0.0, 1.1, label, transform=ax.transAxes + trans,
            fontsize=12, va='bottom', fontweight='bold')    

# ------------Spyders------------
labels = {
    'long': 'Long\n>100 nm',
    'fibrillar': 'Fibrillar\n>0.9 ecc.',
    'dense': 'Dense \n>0.01\nlocs/nm$^2$',
    'bright': 'Bright\n>0.9 A.U.',
    'coloc': 'pSER',
    'spot_count': 'Count'
}

region_map = {
    'Cerebellum': 'Cerebellum',
    'FCx': 'Frontal cortex',
    'Fcx': 'Frontal cortex',
}


for (region, df), ax, shadowax in zip(proportions.groupby('Region'), [axA, axB], [axAi, axBi]):
    plot_radar(dataframe=df, label_col='threshold_category', value_col='proportion', hue_col='Category', palette=palette, ax=ax, err_bars=True, ax_labels=[0, 0.2, 0.4, 0.6], label_dict=labels)
    
    shadowax.annotate(region_map[region], xy=(0.5, -0.158), xycoords='axes fraction', fontsize=8, ha='center', va='center')

handles, labels = axA.get_legend_handles_labels()
by_label = dict(zip(labels, handles))
legend = axAi.legend([by_label[entry] for entry in ['Control', 'SOD1', 'TDP-43']], ['CRL', 'SOD', 'TDP'], bbox_to_anchor=(0, -0.068), frameon=False, loc='center', handletextpad=0.25, handlelength=1.5)
axAi.add_artist(legend)

# ------------LDA------------
plot_lda_scatter(for_LDA, hue='Category', style='Region', palette=palette, ax=axC, s=50)
axC.set(xlim=(-26, 12), ylim=(-12, 6))

handles, labels = axC.get_legend_handles_labels()
by_label = dict(zip(labels, handles))

legend1 = axC.legend([by_label[entry] for entry in ['Control', 'SOD1', 'TDP-43']], ['CRL', 'SOD', 'TDP'], ncol=3, bbox_to_anchor=(0.5, 1.12), frameon=False, columnspacing=0.2, loc='center', handletextpad=0.0025)
legend2 = axC.legend([by_label[entry] for entry in ['Cerebellum', 'FCx']], ['Cerebellum', 'Frontal cortex'], ncol=3, bbox_to_anchor=(0.45, 1.05), frameon=False, columnspacing=0.2, loc='center', handletextpad=0.0025)
axC.add_artist(legend1)
axC.add_artist(legend2)

axC.set_yticks(ticks=np.arange(-10, 6, 5), labels=np.arange(-10, 6, 5))
axC.set_ylabel('Dimension 2', labelpad=-1)

# ------------Figure admin------------
for ax in [axC]:
    ax.spines[['right', 'top']].set_visible(False)
# plt.tight_layout()
# plt.subplots_adjust(wspace=2.5, hspace=5)
plt.savefig(f'{output_folder}figure_4.svg', transparent=True, dpi=1000, bbox_inches="tight")
plt.savefig(f'{output_folder}figure_4.png', transparent=True, dpi=1000, bbox_inches="tight")
plt.show()


"""
Signatures of MND found in nanoscopic TDP-43 particles. **A** Radar visualisation summarising the particle properties derived from soaked samples characterised via SiMPull. Shown is mean ± S.D. for each cohort (CRL *n*=3, SOD1 *n*=2, TDP-43 *n*=5) across six measures; number of puncta as a proportion of the maximum count in any individual donor (count), proportion colabelled with pSER (pSER), proportion with length > 120 nm (long), proportion with eccentricity > 0.9 (fibrillar), proportion with localisation density > 0.01 nm^2 (dense), and proportion whose intensity >0.9 A.U. (bright). **B** Two-dimensional projection of linear discriminant analysis performed on particle properties derived from soaked samples characterised via SiMPull. 
"""