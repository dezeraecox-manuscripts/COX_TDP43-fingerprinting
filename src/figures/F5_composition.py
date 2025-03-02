import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from src.utilities.visualise import bardotplot, plot_volcano
from skimage import io
import subprocess
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches

from loguru import logger

logger.info('Import OK')

output_folder = 'figures/'

if not os.path.exists(output_folder):
    os.makedirs(output_folder)
       
# =================Panel B: flTDP coloc %=================
input_path = 'data/microscopy-tissue/flTDP_pSER-summary.csv'

fl_coloc = pd.read_csv(input_path)
fl_coloc.drop([col for col in fl_coloc.columns.tolist() if 'Unnamed: ' in col], axis=1, inplace=True)

# samples
fl_coloc = fl_coloc[
    (fl_coloc['channel'] == 638) &
    (fl_coloc['Region'] == 'FCx') &
    # (fl_coloc['sample_type'] == 'soaked') &
    (fl_coloc['protein_detected'] == 'flTDP43')
    ].copy().dropna(subset='Category')

fl_coloc = fl_coloc.groupby(['Sample_ID','EXP-number', 'sample_type',	'Region',	'Category',	'Subcategory', 'protein_detected']).mean().reset_index().groupby(['Sample_ID', 'sample_type', 'Region',	'Category',	'Subcategory', 'protein_detected']).mean().reset_index()
       
# ================Panel C: Volcano ================
input_path = 'data/proteomics/normalised-ratios.csv'

ms_ratios = pd.read_csv(input_path)
ms_ratios.drop([col for col in ms_ratios.columns.tolist() if 'Unnamed: ' in col], axis=1, inplace=True)
ms_ratios['-log10_pval'] = ms_ratios['-log10_pval'].replace(np.inf, np.nan).replace(-np.inf, np.nan)

ms_ratios = ms_ratios[(ms_ratios['Category'] == 'TDP-43') & (ms_ratios['sample_type'] == 'Beads') & (ms_ratios['norm_type'] == 'cnorm')].copy().dropna(subset=['log2_ratio', '-log10_pval'])

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

# ================Panel D: Colabelling %================
input_path = 'data/microscopy-tissue/flTDP_colabelling-summary.csv'
colocalised = pd.read_csv(input_path)
colocalised.drop([col for col in colocalised.columns.tolist() if 'Unnamed: ' in col], axis=1, inplace=True)

proteo_coloc = colocalised[
    (colocalised['Category'] == 'TDP-43') &
    (colocalised['Region'] == 'FCx') &
    (colocalised['channel'] == 638)
].copy()

proteo_coloc['sample_type'] = ['h' if stype == 'Homogenate' else ('s' if stype == 'Soaked' else stype) for stype in proteo_coloc['sample_type']]

input_stats = 'data/microscopy-tissue/flTDP_tstats_summary.csv'
tstats = pd.read_csv(input_stats)
tstats.drop([col for col in tstats.columns.tolist() if 'Unnamed: ' in col], axis=1, inplace=True)

tstats = tstats[
    (tstats['category'] == 'TDP-43') &
    (tstats['Region'] == 'FCx')
].copy()

# ================Panel E-G: Binned coloc properties================
input_properties = f'data/microscopy-tissue/SR_coloc-dists.csv'
sr_coloc_spots = pd.read_csv(input_properties)
sr_coloc_spots.drop([col for col in sr_coloc_spots.columns.tolist() if 'Unnamed: ' in col], axis=1, inplace=True)

input_stats = f'data/microscopy-tissue/SR_coloc-dists-stats.csv'
sr_coloc_stats = pd.read_csv(input_stats)
sr_coloc_stats.drop([col for col in sr_coloc_stats.columns.tolist() if 'Unnamed: ' in col], axis=1, inplace=True)

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
    'up': '#FB8B24',
    'down': '#FB8B24',
    'no change': '#FEC086',
    'highlight': '#291720',
    # Coloc
    'Control_col': '#036346',
    'SOD1_col': '#50013D',
    'TDP-43_col': '#FB8B24',
    'homog': '#e08c3f',
    'h': '#FB8B24',
    'soak': '#e08c3f',
    's': '#FB8B24',
    'Control_hatch': '#025038',
    'SOD1_hatch': '#50013D',
    'TDP-43_hatch': '#B45904',
}
    
region_map = {
    'Cerebellum': 'Cerebellum',
    'FCx': 'Frontal cortex',
    'Fcx': 'Frontal cortex',
    'H': 'Embedded',
    'Homogenate': 'Embedded',
    'S': 'Diffusible',
    'Soaked': 'Diffusible',
}

prop_map = {
    'smoothed_length': 'Length (nm)', 
    'eccentricity': 'Eccentricity', 
    'density': 'Localisation density\n(×10$^{-4}$ per nm$^2$)'
}


fig = plt.figure(figsize=(18.8*cm, 6.1*3*cm))
gs = fig.add_gridspec(nrows=3, ncols=1, hspace=0.25)
gs1 = gridspec.GridSpecFromSubplotSpec(1, 10, subplot_spec=gs[0:1, 0:1], wspace=1.5)
gs2 = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=gs[1:2, 0:1], wspace=0.57)
gs3 = gridspec.GridSpecFromSubplotSpec(4, 3, subplot_spec=gs[2:3, 0:1], wspace=0.57)

axA = fig.add_subplot(gs1[0:2])
axAi = fig.add_subplot(gs1[0:2])
axBi = fig.add_subplot(gs1[2:4])
axBii = fig.add_subplot(gs1[4:6])
axC = fig.add_subplot(gs1[6:10])
axD = fig.add_subplot(gs2[:])

axE = fig.add_subplot(gs3[0:4, 0:1])
axE1 = fig.add_subplot(gs3[0:1, 0:1])
axE2 = fig.add_subplot(gs3[1:2, 0:1])
axE3 = fig.add_subplot(gs3[2:3, 0:1])
axE4 = fig.add_subplot(gs3[3:4, 0:1])

axF = fig.add_subplot(gs3[0:4, 1:2])
axF1 = fig.add_subplot(gs3[0:1, 1:2])
axF2 = fig.add_subplot(gs3[1:2, 1:2])
axF3 = fig.add_subplot(gs3[2:3, 1:2])
axF4 = fig.add_subplot(gs3[3:4, 1:2])

axG = fig.add_subplot(gs3[0:4, 2:3])
axG1 = fig.add_subplot(gs3[0:1, 2:3])
axG2 = fig.add_subplot(gs3[1:2, 2:3])
axG3 = fig.add_subplot(gs3[2:3, 2:3])
axG4 = fig.add_subplot(gs3[3:4, 2:3])

for ax, label, xpos in zip([axAi, axBi, axC, axD, axE1, axF1, axG1], ['A', 'B', 'C', 'D', 'E', 'F', 'G'], [-35, -30, -30, -35, -35, -30, -30]):
# label physical distance to the left and up:
    trans = mtransforms.ScaledTranslation(xpos/72, -8/72, fig.dpi_scale_trans)
    ax.text(0.0, 1.1, label, transform=ax.transAxes + trans,
            fontsize=12, va='bottom', fontweight='bold')

# -------flTDP coloc-------
# The Holm–Bonferroni method also controls the FWER at \alpha , but with a lower increase of type II error risk than the classical Bonferroni method.
for (ax, stype) in zip([axBi, axBii], ['soaked', 'homogenate']):
    bardotplot(
        fl_coloc[fl_coloc['sample_type'] == stype].copy(), xcol='Category', ycol='proportion_coloc', order=['Control', 'SOD1', 'TDP-43'], 
        hue=None, hue_order=None, scat_hue=None, scat_hue_order=None, 
        palette=palette, 
        xlabel='', ylabel='Colabelled (%)', 
        pairs=[('Control', 'TDP-43')],
        correction='holm-bonferroni', 
        xticks=None, 
        groups=None, group_label_y=-0.18, group_line_y=-0.05, 
        ax=ax, legend='', 
        dot_size=5, cap_size=0.2, cap_width=1
    )

    ax.axhline(fl_coloc[fl_coloc['sample_type'] == stype].copy()['chance_proportion_coloc'].mean(), linestyle='--', color='lightgrey')
    ax.set_title(region_map[stype.capitalize()], fontsize=8)
    ax.set_xticks(np.arange(0, 3), labels=['CRL  ', 'SOD', '   TDP'], fontsize=8)
    
    if stype == 'homogenate':
        for category, patch in zip(['Control', 'SOD1', 'TDP-43'], ax.patches):
            patch.set_edgecolor(palette[f'{category}_hatch'])
            patch.set_hatch('///')
            patch.set_linewidth(0)
        
axBi.set_ylabel('Colabelled puncta (%)', labelpad=-1, fontsize=8)
axBii.set_ylabel('', labelpad=-1)

# ------------Volcano plot------------
plot_volcano(ms_ratios[ms_ratios['thresh_cat'] == 'no change'], xcol='log2_ratio', ycol='-log10_pval', palette=palette, hue_col='thresh_cat', ax=axC, ythresh=1.3, xthresh=[-1, 1], alpha=0.7)
plot_volcano(ms_ratios[ms_ratios['thresh_cat'].isin(['up', 'down'])], xcol='log2_ratio', ycol='-log10_pval', palette=palette, hue_col='thresh_cat', ax=axC, ythresh=1.3, xthresh=[-1, 1], alpha=0.9)
plot_volcano(ms_ratios[ms_ratios['thresh_cat'] == 'highlight'], xcol='log2_ratio', ycol='-log10_pval', palette=palette, hue_col='thresh_cat', ax=axC, ythresh=1.3, xthresh=[-1, 1])
axC.legend('', frameon=False)
axC.set_xlabel('Log$_2$ Corrected ratio', labelpad=-1)
axC.set(ylabel='-Log$_{10}$ p-value', xlim=(-3.6, 3.6))

# Add annotations for highlighted spots

positions = ms_ratios[ms_ratios['thresh_cat'] == 'highlight'][['Majority protein IDs', 'log2_ratio', '-log10_pval']].copy()
positions = {key: (round(pos1, 3), round(pos2, 3)) for key, pos1, pos2 in positions.values}

label_positions = {
    # 'Q13148': (0, 0), #TDP
    'P35637': (2.0, 3.5), #FUS
    'P21926': (-2.2, 4.6), #CD9
    'P0DMV9;P0DMV8' : (-1.8, 1.5),  #HSPs
    'O60763' : (3, 1.75),  #USO1
    'P27824': (-1, 0.5), #CANX
    'P06748': (-2.5, 3), # Nucleophosmin
}

for label in label_positions:
    label_pos = label_positions[label]
    arrow_pos = positions[label]
    ha = 'left' if label_pos[0] > 0 else 'right'
    offset = 0.1 if label_pos[0] < 0 else -0.1
    axC.annotate(highlights[label], xy=label_positions[label], fontsize=5, ha=ha, va='center', zorder=100000)
    axC.plot([label_pos[0]+offset, arrow_pos[0]], [label_pos[1], arrow_pos[1]], color=palette['highlight'], linewidth=0.5, zorder=100000)
    
# ------------Colabelling ------------
order = ['USO1', 'NPM1', 'HSP70', 'FUS', 'CD9', 'Calnexin', 'CD63', 'TSG101']
bardotplot(
    proteo_coloc, xcol='protein_detected', ycol='change_coloc', order=order, 
    hue='sample_type', hue_order=['s', 'h'], scat_hue='sample_type', scat_hue_order=['s', 'h'], 
    palette=palette, 
    xlabel='Target', ylabel='Colabelled (%)', 
    pairs=None,
    correction=None, 
    xticks=None, 
    groups=None, group_label_y=-0.18, group_line_y=-0.05, 
    ax=axD, legend='', 
    dot_size=5, cap_size=0.2, cap_width=1
)

# Add annotations
formatted_pvalues = [tstats[(tstats['protein_detected'] == protein) & (tstats['sample_type'] == stype)]['sign'].tolist()[0] for protein in ['USO1', 'NPM1', 'HSP70', 'FUS', 'CD9', 'Calnexin', 'CD63', 'TSG101'] for stype in ['s', 'h']] 
for x in range(8):
    axD.annotate(formatted_pvalues[x*2], xy=(x-0.2, 1.4), fontsize=8, ha='center')
    axD.annotate(formatted_pvalues[x*2+1], xy=(x+0.2, 1.4), fontsize=8, ha='center')
    

for stype, patch in zip(['s']*len(order) + ['h']*len(order), axD.patches):
    if stype == 'h':
        patch.set_facecolor(palette['homog'])
        patch.set_edgecolor(palette['TDP-43_hatch'])
        patch.set_linewidth(0)
        patch.set_hatch('///')
    else:
        patch.set_facecolor(palette['homog'])
        patch.set_edgecolor('#fff')
        patch.set_linewidth(1)
        
s_patch = mpatches.Patch(facecolor=palette['homog'], label='Diffusible', edgecolor='#fff', linewidth=1)
h_patch = mpatches.Patch(facecolor=palette['homog'], label='Embedded', edgecolor=palette['TDP-43_hatch'], linewidth=0, hatch='///')
# axD.legend(handles=[s_patch, h_patch], frameon=False, ncols=2)

axD.axhline(1, linestyle='--', color='lightgrey')
axD.set_ylabel('$\Delta$ Colabelled', labelpad=-1, fontsize=8)
axD.set_xlabel('')
axD.set_yticks([0.75, 1.0, 1.25, 1.5])
axD.set_ylim(0.75, 1.5)
axD.set_xlim(-0.6, 7.6)
handles, labels = axD.get_legend_handles_labels()
by_label = dict(zip(reversed(labels), reversed(handles)))
axD.legend([s_patch, h_patch], ['Diffusible', 'Embedded'], frameon=False, loc='upper left', handlelength=1., handletextpad=0.25, bbox_to_anchor=(-0.005, 1.12), ncol=2)

# ------------Coloc v Noncoloc distributions------------
for ax in [axE, axF, axG]:
    ax.axis('off')
    ax.annotate('Density', xy=(-0.1, 0.5), xycoords='axes fraction', rotation=90, ha='center', va='center')

category = 'TDP-43'
color = palette['homog']
for prop, axes in zip(['smoothed_length', 'eccentricity', 'density'], [(axE1, axE2, axE3, axE4), (axF1, axF2, axF3, axF4), (axG1, axG2, axG3, axG4)]):
    
    if prop == 'density':
        xlims = (-0.005, 0.05)
        label_pos = 1
        label_orient = 'right'
    if prop == 'eccentricity':
        xlims = (-0.1, 1.1)
        label_pos = 0.05
        label_orient = 'left'
    if prop == 'smoothed_length':
        xlims = (-10, 500)
        label_pos = 1
        label_orient = 'right'
        
    for j, ((detected), data) in enumerate(sr_coloc_spots.groupby(['protein_detected'])):
            
        data = data[(data[prop] <= xlims[1]) & (data['Category'] == category) & (data['coloc?'] == 1)].copy()
        noncoloc_data = sr_coloc_spots[(sr_coloc_spots['Category'] == category) & (sr_coloc_spots['coloc?'] == 0) & (sr_coloc_spots[prop] <= xlims[1])].copy()
        
        sns.kdeplot(
            data=data,
            x=prop,
            ax=axes[j],
            color='black',
            linewidth=1
        )
        sns.kdeplot(
            data=noncoloc_data,
            x=prop,
            ax=axes[j],
            color=color,
            linewidth=1,
            fill=True,
            alpha=0.5
        )
        
        axes[j].set_xlim(*xlims)
        axes[j].legend('', frameon=False)
        axes[j].set_xlabel('')
        axes[j].set_ylabel('')
        if j < 3:
            axes[j].set_xticklabels([])
        axes[j].set_yticklabels([])
        axes[j].set_yticks([])
        axes[j].annotate(f'+{detected}', xy=(label_pos, 0.8), xycoords='axes fraction', color='black', fontsize=8, va='center', ha=label_orient)
        
        ks_stat = sr_coloc_stats[(sr_coloc_stats['Category'] == category) & (sr_coloc_stats['parameter'] == prop) & (sr_coloc_stats['protein_detected'] == detected)]['sign'].tolist()[0]
        axes[j].annotate(ks_stat, xy=(1.1, 0.5), xycoords='axes fraction', color='black', fontsize=8, va='center', ha='center')
        
    axes[0].annotate('TDP-43', xy=(label_pos, 1.1), xycoords='axes fraction', color=color, fontsize=8, va='center', ha=label_orient)
    axes[3].set_xlabel(prop_map[prop])

# ------------Figure admin------------
for ax in [axBi, axBii, axC, axD, axE1, axE2, axE3, axE4, axF1, axF2, axF3, axF4, axG1, axG2, axG3, axG4]:
    ax.spines[['right', 'top']].set_visible(False)
plt.savefig(f'{output_folder}figure_5.svg', transparent=True, dpi=1000)
plt.savefig(f'{output_folder}figure_5.png', transparent=True, dpi=1000)
plt.show()
