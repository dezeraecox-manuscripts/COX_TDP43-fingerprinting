import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from src.utilities.visualise import bardotplot, plot_sankey, sankey_converter
from skimage import io
import matplotlib.gridspec as gridspec
import matplotlib.patches as patches

from loguru import logger
logger.info('Import OK')

output_folder = 'figures/'

if not os.path.exists(output_folder):
    os.makedirs(output_folder)

palette = {
    'Control': '#04A777',
    'short_Control': '#D7EEE7',
    'long_Control': '#7DD0B7',
    'SOD1': '#820263',
    'short_SOD1': '#EBD0E4',
    'long_SOD1': '#C26FAD',
    'TDP-43': '#FB8B24',
    'short_TDP-43': '#FFDAB7',
    'long_TDP-43': '#FBA04B',
}
region_map = {
    'Cerebellum': 'Cerebellum',
    'FCx': 'Frontal cortex',
    'Fcx': 'Frontal cortex',
}

# ================Panel : measurements================

input_path = 'data/microscopy-tissue/SR_measurement_summary.csv'

# Read in dataset
properties = pd.read_csv(input_path)
properties.drop([col for col in properties.columns.tolist() if 'Unnamed: ' in col], axis=1, inplace=True)

# Filter out only cluster props for tissue samples
mean_properties = properties[
    (properties['Sample_ID'] != 'Image controls') &
    (properties['prop_type'] == 'cluster') &
    (properties['channel'] == 638) &
    (properties['sample_type'] == 'soak')
    ].copy()

mean_properties['density'] = mean_properties['density'] * 10000
# ================Panel G: proportions================
input_path = 'data/microscopy-tissue/SR_flows.csv'

# Read in dataset
flows = pd.read_csv(input_path)
flows.drop([col for col in flows.columns.tolist() if 'Unnamed: ' in col], axis=1, inplace=True)
        
flow_data = flows[
    (flows['sample_type'] == 'soak') &
    (flows['prop_type'] == 'cluster')
].copy()

for (sample_type, region), df in flow_data.groupby(['sample_type', 'Region']):
    plot_sankey(flow_data=df, palette=palette, output_path=f'data/microscopy-tissue/sankey_{sample_type}_{region}.svg', h=500, w=750)

flow_data['cat_eccentricity'] = flow_data['cat_eccentricity'].str.capitalize().str.replace('Fibril', 'Fibrillar')
flow_data['cat_smoothed_length'] = flow_data['cat_smoothed_length'].str.capitalize()
flow_data['cat_density'] = flow_data['cat_density'].str.capitalize()

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


fig = plt.figure(figsize=(18.8*cm, 6.1*3*cm))

gs = fig.add_gridspec(nrows=3, ncols=3, wspace=0.43, hspace=0.5)

gsA = gridspec.GridSpecFromSubplotSpec(1, 5, subplot_spec=gs[0:1, 0:3], wspace=0)
gsB = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs[1:2, 0:1], wspace=0.45)
gsC = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs[1:2, 1:2], wspace=0.57)
gsD = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs[1:2, 2:3], wspace=0.45)

gsE = gridspec.GridSpecFromSubplotSpec(1, 6, subplot_spec=gs[2:3, 0:3], wspace=0.45)

axA = fig.add_subplot(gsA[0:1, 0:5])
axA1 = fig.add_subplot(gsA[0:1, 0:1])
axA3 = fig.add_subplot(gsA[0:1, 1:2])
axA2 = fig.add_subplot(gsA[0:1, 2:3])
axA4 = fig.add_subplot(gsA[0:1, 3:4])
axA5 = fig.add_subplot(gsA[0:1, 4:5])

axBi = fig.add_subplot(gsB[0:1])
axBii = fig.add_subplot(gsB[1:2])
axCi = fig.add_subplot(gsC[0:1])
axCii = fig.add_subplot(gsC[1:2])
axDi = fig.add_subplot(gsD[0:1])
axDii = fig.add_subplot(gsD[1:2])

axE = fig.add_subplot(gsE[0:3])
axEi = fig.add_subplot(gsE[0:3])
axF = fig.add_subplot(gsE[3:6])
axFi = fig.add_subplot(gsE[3:6])

for ax, label, (xpos, ypos) in zip([axA, axBi, axCi, axDi, axEi, axFi], ['A', 'B', 'C', 'D', 'E', 'F'], [(-25, 2), (-25, -8), (-25, -8), (-25, -8), (-25, 4), (-25, 4)]):
# label physical distance to the left and up:
    trans = mtransforms.ScaledTranslation(xpos/72, ypos/72, fig.dpi_scale_trans)
    ax.text(0.0, 1.1, label, transform=ax.transAxes + trans,
            fontsize=12, va='bottom', fontweight='bold')


axA.axis('off')

# -----------Barplots-----------
for label, props, axes in zip(['Length (nm)', 'Eccentricity', 'Localisation density\n(×10$^{-4}$ per nm$^2$)'], ['smoothed_length', 'eccentricity', 'density'], [(axBi, axBii), (axCi, axCii), (axDi, axDii)]):
    for (ax, (region, df)) in zip(axes, mean_properties.groupby('Region')):
        bardotplot(
            df, xcol='Category', ycol=props, order=['Control', 'SOD1', 'TDP-43'], 
            hue=None, hue_order=None, scat_hue=None, scat_hue_order=None, 
            palette=palette, 
            xlabel='', ylabel=label, 
            pairs=[('Control', 'TDP-43')], correction='holm-bonferroni', 
            xticks=None, 
            groups=None, group_label_y=-0.18, group_line_y=-0.05, 
            ax=ax, legend='', 
            dot_size=5, cap_size=0.2, cap_width=1
        )

        ax.set_title(region_map[region.capitalize()], fontsize=8)
        ax.set_xticks(np.arange(0, 3), labels=['CRL   ', 'SOD', '   TDP'], fontsize=8)
    if props == 'eccentricity':
        pad = 1
    else:
        pad = -1.2
    axes[0].set_ylabel(label, labelpad=pad)
    axes[1].set_ylabel('')

# -------------Proportions-------------
axes = [axE, axF, axEi, axFi]
for ax in axes:
    ax.axis('off')

for i, ((sample_type, region), df) in enumerate(flow_data.groupby(['sample_type', 'Region'])):
    ax = axes[i]
    img = sankey_converter(input_path=f'data/microscopy-tissue/sankey_{sample_type}_{region}.svg', w=625)
    ax.imshow(img)
    
    ax.set_position([(0.43*(i)+0.11), 0.08, 0.34, 0.32])
    
    label_pos = {}
    for col in ['cat_eccentricity', 'cat_smoothed_length', 'cat_density']:
        pos = df.groupby(['Category', col]).sum()['proportion'].reset_index()
        pos['key'] = [f'{label}_{cat}' for label, cat in pos[[col, 'Category']].values]
        label_pos.update(dict(pos[['key', 'proportion']].values))
    
    for i, (cat, title_pos) in enumerate(zip(['TDP-43', 'SOD1', 'Control'], [215, 800, 1350])):
        ax.annotate(cat, xy=(-50, title_pos), xycoords='data', rotation=90, ha='center', va='center', annotation_clip=False, fontsize=8)
        
        labels = ['Long', 'Short', 'Fibrillar', 'Globular', 'Dense', 'Sparse']
        for j, (label, xpos) in enumerate(zip(
            labels, 
            [1200, 1200, 2200, 2200, 3250, 3250]
            )):
            if j % 2 == 0:
                ypos = label_pos[f'{label}_{cat}']/2*550+i*550+20
            else:
                ypos = (label_pos[f'{labels[j-1]}_{cat}']*550) + (label_pos[f'{label}_{cat}']*550)/2+i*550+20
                

            ax.annotate(label, xy=(xpos, ypos), xycoords='data', ha='left', va='center', annotation_clip=False, fontsize=5)
    ax.annotate(f'{region_map[region]}', xy=(1650, -150), xycoords='data', ha='center', va='center', annotation_clip=False, fontsize=8)
    
    for label, xpos in zip(['Length', 'Eccentricity', 'Density'], [1150, 2150, 3200]):
        ax.annotate(label, xy=(xpos, 1700), xycoords='data', ha='center', va='center', annotation_clip=False, fontsize=8)

# Add dotted boxes to highlight populations of interest

from matplotlib import patches
axes = [axE, axF]
for ax in axes:
    # Create a Rectangle patch
    rects = [
        patches.Rectangle((2135, 185), 1050, 150, linewidth=0.5, edgecolor='grey', facecolor='none', linestyle='--'),
        patches.Rectangle((2135, 760), 1050, 140, linewidth=0.5, edgecolor='grey', facecolor='none', linestyle='--'),
        patches.Rectangle((2135, 1375), 1050, 75, linewidth=0.5, edgecolor='grey', facecolor='none', linestyle='--'),
    ]

    # Add the patch to the Axes
    for rect in rects:
        ax.add_patch(rect)

# ------------Figure admin------------

for ax in [axBi, axBii, axCi, axCii, axDi, axDii]:
    ax.spines[['right', 'top']].set_visible(False)
plt.savefig(f'{output_folder}figure_3.svg', transparent=True)
plt.savefig(f'{output_folder}figure_3.png', transparent=True, dpi=1000)
plt.show()
