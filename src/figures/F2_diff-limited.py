import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from src.utilities.visualise import bardotplot, plot_interpolated_ecdf
from src.utilities.analyze import distribution_ks, fitting_ecdfs
from skimage import io
import subprocess

from loguru import logger

logger.info('Import OK')

output_folder = 'figures/'

if not os.path.exists(output_folder):
    os.makedirs(output_folder)
       

# =================Panel B: flTDP spot count=================
input_path = 'data/microscopy-tissue/flTDP_normalised-counts.csv'

fl_spots = pd.read_csv(input_path)
fl_spots.drop([col for col in fl_spots.columns.tolist() if 'Unnamed: ' in col], axis=1, inplace=True)

# Controls
fl_controls = fl_spots[(fl_spots['detect'] == 'IgG') & (fl_spots['sample'].str.contains('h-MIX'))].copy()

# samples
fl_spots = fl_spots[
    (fl_spots['protein_detected'] != 'None') &
    (fl_spots['sample_type'] == 'soak') &
    (fl_spots['protein_detected'] == 'flTDP')
    ].copy().dropna(subset='Category')

fl_spots = fl_spots.groupby(['tissue_ID',	'sample_type',	'Region',	'Category',	'Subcategory', 'protein_detected']).mean().reset_index()

# =================Panel C-E: Brightness=================
input_path = 'data/microscopy-tissue/flTDP_brightness-summary.csv'

fl_brightness = pd.read_csv(input_path)
fl_brightness.drop([col for col in fl_brightness.columns.tolist() if 'Unnamed: ' in col], axis=1, inplace=True)

# samples
fl_brightness = fl_brightness[
    (fl_brightness['channel'] == 638) &
    (fl_brightness['sample_type'] == 'soak')
    ].copy().dropna(subset='Category')

fl_brightness['scaled_intensity'] = fl_brightness['scaled_intensity'] / 1000

fl_brightness_mean = fl_brightness.groupby(['Sample_ID','EXP-number', 'sample_type', 'Region',	'Category',	'Subcategory']).mean().reset_index().groupby(['Sample_ID', 'sample_type', 'Region',	'Category',	'Subcategory']).mean().reset_index()

fl_brightness_ecdf = fitting_ecdfs(fl_brightness, group_cols=['Sample_ID', 'sample_type', 'Region', 'Category',	'Subcategory'], val_col='scaled_intensity').dropna(subset=['scaled_intensity'])

fl_brightness_ecdf_ks = pd.concat([distribution_ks(dataframe=fl_brightness, group_cols=['sample_type', 'Region'], cat_col='Category', cat1='Control', cat2=comparison, parameter_col='scaled_intensity', num_samples=1000, seed=1) for comparison in ['SOD1', 'TDP-43']])



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
region_map = {
    'Cerebellum': 'Cerebellum',
    'FCx': 'Frontal cortex',
    'Fcx': 'Frontal cortex',
}


fig = plt.figure(figsize=(12.53*cm, 6.1*2*cm))

gs = fig.add_gridspec(nrows=2, ncols=8, wspace=3.0, hspace=0.35)


axBi = fig.add_subplot(gs[0:1, 0:2])
axBii = fig.add_subplot(gs[0:1, 2:4])
axCi = fig.add_subplot(gs[0:1, 4:6])
axCii = fig.add_subplot(gs[0:1, 6:8])

axD = fig.add_subplot(gs[1:2, 0:4])
axE = fig.add_subplot(gs[1:2, 4:8])

for ax, label, xpos in zip([axBi, axCi, axD, axE], ['A', 'B', 'C', 'D'], [-30, -28, -30, -28]):
# label physical distance to the left and up:
    trans = mtransforms.ScaledTranslation(xpos/72, -8/72, fig.dpi_scale_trans)
    ax.text(0.0, 1.1, label, transform=ax.transAxes + trans,
            fontsize=12, va='bottom', fontweight='bold')

# -------flTDP counts-------
# The Holm–Bonferroni method also controls the FWER at \alpha , but with a lower increase of type II error risk than the classical Bonferroni method.
for (ax, (region, df)) in zip([axBi, axBii], fl_spots.groupby('Region')):
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

    ax.axhline(fl_controls['spots_count'].mean(), linestyle='--', color='lightgrey')
    ax.set_title(region_map[region.capitalize()], fontsize=8)
    ax.set_xticks(np.arange(0, 3), labels=['CRL   ', 'SOD', '   TDP'], fontsize=8)
axBi.set_ylabel('Number of puncta', labelpad=0, fontsize=8)
axBii.set_ylabel('', labelpad=0)
    
# -------flTDP brightness mean-------
for (ax, (region, df)) in zip([axCi, axCii], fl_brightness_mean.groupby('Region')):
    bardotplot(
        df, xcol='Category', ycol='scaled_intensity', order=['Control', 'SOD1', 'TDP-43'], 
        hue=None, hue_order=None, scat_hue=None, scat_hue_order=None, 
        palette=palette, 
        xlabel='', ylabel='Intensity (AU)', 
        pairs=[('Control', 'TDP-43')], # ('Control', 'SOD1'), ], 
        correction='holm-bonferroni', 
        xticks=None, 
        groups=None, group_label_y=-0.18, group_line_y=-0.05, 
        ax=ax, legend='', 
        dot_size=5, cap_size=0.2, cap_width=1
    )
    ax.set_title(region_map[region.capitalize()], fontsize=8)
    ax.set_xticks(np.arange(0, 3), labels=['CRL   ', 'SOD', '   TDP'], fontsize=8)
axCi.set_ylabel('Intensity (AU)', labelpad=0, fontsize=8)
axCii.set_ylabel('', labelpad=0)

# -------flTDP brightness distribution-------
for (ax, (region, df)) in zip([axD, axE], fl_brightness_ecdf.groupby('Region')):

    plot_interpolated_ecdf(df, ycol='scaled_intensity', huecol='Category', palette=palette, ax=ax, orientation='h')
    ax.legend(fontsize=8, frameon=False, handlelength=0.75, labelspacing=0.1, handletextpad=0.2, borderaxespad=0.05)
    ax.set_ylim(0.84, 1.01)
    ax.set_yticks(np.linspace(0.85, 1.0, 4), labels=[0.85, 0.90, 0.95, 1.00])
    ax.set_xlim(0, 11)
    ax.set_title(region_map[region.capitalize()], fontsize=8)
    ax.set_ylabel('Proportion', fontsize=8, labelpad=0)
    ax.set_xlabel('Intensity (AU)', fontsize=8)
    
    for combination, lowy, xpos in zip(['Control_TDP-43'], [0.8525], [6.8]):
        topy = 0.878
        stat = fl_brightness_ecdf_ks[(fl_brightness_ecdf_ks['comparison'] == combination) & (fl_brightness_ecdf_ks['Region'] == region)]['sign'].tolist()[0]
        bar_tips = xpos + 0.2
        ax.plot(
            [bar_tips, xpos, xpos, bar_tips],
            [topy, topy, lowy, lowy], lw=0.5, c='black'
        )
        stat_offset = 0.5 if stat == 'ns' else 0.16
        ax.annotate(stat, xy=(xpos-stat_offset, (lowy+(topy-lowy)/2)), rotation=90, fontsize=8, ha= 'center', va='center')
    

# ------------Figure admin------------
for ax in [axBi, axBii, axCi, axCii, axD, axE]:
    ax.spines[['right', 'top']].set_visible(False)
plt.savefig(f'{output_folder}figure_2.svg', transparent=True, dpi=1000)
plt.savefig(f'{output_folder}figure_2.png', transparent=True, dpi=1000)
plt.show()
