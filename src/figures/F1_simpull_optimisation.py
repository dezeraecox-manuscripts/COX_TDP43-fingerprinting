
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from src.utilities.visualise import bardotplot, plot_interpolated_ecdf
from src.utilities.analyze import fitting_ecdfs
from loguru import logger
from skimage import io, exposure
import subprocess
from scipy import stats

logger.info('Import OK')

output_folder = 'figures/'

if not os.path.exists(output_folder):
    os.makedirs(output_folder)
    


# ================Panel C: ThT Kinetics================
input_kinetics = 'data/microscopy-RRMs/complete_fluorescence.csv'

# Read in data
fluorescence = pd.read_csv(input_kinetics)
fluorescence.drop([col for col in fluorescence.columns.tolist() if 'Unnamed: ' in col], axis=1, inplace=True)

fluorescence['norm_fluorescence'] = (fluorescence['fluorescence'] - fluorescence[fluorescence['sample'].str.contains('Tht')]['fluorescence'].min())/ fluorescence['fluorescence'].max() 
fluorescence['time_h'] = fluorescence['time_combined'] / 60
fluorescence = fluorescence[fluorescence['protein'] == 'RRM1+2'].copy()

# ================Panel D: SiMPull diff limited================
input_count_tht = 'data/microscopy-RRMs/summary_488.csv'

# Read in summary FOV data
spots_tht = pd.read_csv(f'{input_count_tht}')
spots_tht.drop([col for col in spots_tht.columns.tolist() if 'Unnamed: ' in col], axis=1, inplace=True)
spots_tht = spots_tht[['replicate', 'HSA+ThT-t0-10000', 'RRM1+2-t0-5000', 'RRM1+2-t72-5000']].copy()
spots_tht = pd.melt(spots_tht, id_vars='replicate',  var_name='sample', value_name='tht_count')

input_count_ab = 'data/microscopy-RRMs/summary_638.csv'

# Read in summary FOV data
spots_ab = pd.read_csv(f'{input_count_ab}')
spots_ab.drop([col for col in spots_ab.columns.tolist() if 'Unnamed: ' in col], axis=1, inplace=True)
spots_ab = spots_ab[['replicate', 'HSA+ThT-t0-10000', 'RRM1+2-t0-5000', 'RRM1+2-t72-5000']].copy()
spots_ab = pd.melt(spots_ab, id_vars='replicate', var_name=['sample'], value_name='flTDP_count')

spots = pd.merge(spots_tht, spots_ab, on=['replicate', 'sample'], how='outer')
spots['label'] = spots['sample'].map({'HSA+ThT-t0-10000': 'Control', 'RRM1+2-t0-5000': 'Before', 'RRM1+2-t72-5000': 'After'})

# ===============Panel E: SiMPull super-res size===============
input_lengths= 'data/microscopy-RRMs/branch_lengths.csv'

lengths = pd.read_csv(input_lengths)
lengths.drop([col for col in lengths.columns.tolist() if 'Unnamed: ' in col], axis=1, inplace=True)

lengths = lengths[
    (lengths['protein'] == 'RRM1+2') &
    (lengths['concentration'] == 5000) &
    (lengths['branch-distance'] > 0)
    ].copy()
lengths['replicate'] = lengths['well_info'].str[4:8]
lengths['sample'] = lengths['timepoint'].map({'t0': 'Before', 't72': 'After'})
lengths['log_length'] = np.log10(lengths['branch-distance'])

ks_stat, ks_pval = stats.ks_2samp(
            lengths[lengths['sample'] == 'Before']['branch-distance'].values,
            lengths[lengths['sample'] == 'After']['branch-distance'].values
        )
        
lengths_ecdf = fitting_ecdfs(lengths, group_cols=['sample', 'replicate', 'timepoint'], val_col='branch-distance').dropna(subset=['branch-distance'])

# ================Compile figure================
import matplotlib
import matplotlib.transforms as mtransforms
font = {'family' : 'arial',
'weight' : 'normal',
'size'   : 8 }
matplotlib.rc('font', **font)
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['figure.dpi'] = 300
cm = 1/2.54  # centimeters in inches


# Visualise average fluorescence
palette = {
    'RRM1+2-10': '#33658A',
    'RRM1+2-5': '#668CA8',
    'RRM1+2-1': '#81A0B7',
    'Before': '#234E69',
    'After': '#234E69',
}

fig = plt.figure(figsize=(18.8*cm, 6.1*2*cm))

gs = fig.add_gridspec(nrows=2, ncols=6, wspace=0.85, hspace=0.35)
axA = fig.add_subplot(gs[0:1, 0:2])
axAi = fig.add_subplot(gs[0:1, 0:2])
axBi = fig.add_subplot(gs[0:1, 2:4])
axBii = fig.add_subplot(gs[0:1, 4:6])
axC = fig.add_subplot(gs[1:2, 0:2])
axDi = fig.add_subplot(gs[1:2, 2:3])
axDii = fig.add_subplot(gs[1:2, 3:4])
axE = fig.add_subplot(gs[1:2, 4:6])

# -------SiMPull schematic-------
axA.axis('off')
axAi.axis('off')
axAi.set_position([0.065, 0.5, 0.3, 0.4])

# --------Example images--------
import matplotlib.patches as patches

recti = patches.Rectangle((425, 475), 46.7, 2, linewidth=1, edgecolor='white', facecolor='white')

axBi.annotate('Before', xy=(0.05, 0.9), xycoords='axes fraction', fontsize=8, color='white')
axBi.add_patch(recti)
axBi.axis('off')

rectii = patches.Rectangle((425, 475), 46.7, 2, linewidth=1, edgecolor='white', facecolor='white')

axBii.annotate('After', xy=(0.05, 0.9), xycoords='axes fraction', fontsize=8, color='white')
axBii.add_patch(rectii)
axBii.axis('off')

# ------------Kinetics------------
sns.lineplot(
    data=fluorescence,
    x='time_h',
    y='norm_fluorescence',
    hue='sample',
    palette=palette,
    ci='sd',
    hue_order=[f'RRM1+2-10', f'RRM1+2-5', f'RRM1+2-1'],
    ax=axC
)
axC.set_ylim(0, 1)
axC.set_xlabel('Time (h)')
axC.set_ylabel('ThT Fluorescence (A.U.)')
handles, labels = axC.get_legend_handles_labels()
by_label = dict(zip([f'{label.split("-")[1]} µM' for label in labels], handles))
axC.legend(by_label.values(), by_label.keys(), frameon=False, handlelength=1)

# ------------Spot count------------

buffer_control = dict(spots[spots['label'] == 'Control'].copy().mean().T.reset_index().values)
spots = spots[spots['sample'] != 'Control'].copy()

bardotplot(
    spots, xcol='label', ycol='tht_count', order=['Before', 'After'], 
    hue=None, hue_order=None, scat_hue=None, scat_hue_order=None, 
    palette=palette, 
    xlabel='', ylabel='Number of puncta', 
    pairs=[('Before', 'After')], correction=None, 
    xticks=None, 
    groups=None, group_label_y=-0.18, group_line_y=-0.05, 
    ax=axDi, legend='', 
    dot_size=5, cap_size=0.2, cap_width=1
    )
axDi.axhline(buffer_control['tht_count'], linestyle='--', color='lightgrey', linewidth=1)
axDi.set_title('ThT')
axDi.set_ylabel('Number of puncta', labelpad=1)

bardotplot(
    spots, xcol='label', ycol='flTDP_count', order=['Before', 'After'], 
    hue=None, hue_order=None, scat_hue=None, scat_hue_order=None, 
    palette=palette, 
    xlabel='', ylabel='Number of puncta', 
    pairs=[('Before', 'After')], correction=None, 
    xticks=None, 
    groups=None, group_label_y=-0.18, group_line_y=-0.05, 
    ax=axDii, legend='', 
    dot_size=5, cap_size=0.2, cap_width=2
    )
axDii.axhline(buffer_control['flTDP_count'], linestyle='--', color='lightgrey', linewidth=1)
axDii.set_title('flTDP')
axDii.set_ylabel('', labelpad=0.15)

# ------------Size------------
plot_interpolated_ecdf(
    fitted_ecdfs=lengths_ecdf[(lengths_ecdf['sample'] == 'Before') & (lengths_ecdf['ecdf'] < 1.0)].copy().dropna(), ycol='branch-distance',
    huecol='sample', palette=palette, ax=axE, orientation='h',
    linestyle='--')
plot_interpolated_ecdf(
    fitted_ecdfs=lengths_ecdf[(lengths_ecdf['sample'] == 'After')  & (lengths_ecdf['ecdf'] < 1.0)].copy().dropna(), ycol='branch-distance',
    huecol='sample', palette=palette, ax=axE, orientation='h',
    linestyle='-'
    )
axE.set_xlim(-10, 1000)
axE.set_xlabel('Length (nm)')
axE.set_ylabel('Proportion')
axE.set_title('flTDP')
axE.set_xticks(np.arange(0, 1001, 200))
axE.legend(frameon=False, handlelength=1.5)

lengths.groupby('sample').mean()

# ------------Figure admin------------
for ax in [axA, axBi, axBii, axC, axDi, axDii, axE]:
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    
for ax, label in zip([axA, axBi, axC, axDi, axE], ['A', 'B', 'C', 'D', 'E']):
# label physical distance to the left and up:
    trans = mtransforms.ScaledTranslation(-35/72, -11/72, fig.dpi_scale_trans)
    ax.text(0.0, 1.1, label, transform=ax.transAxes + trans,
            fontsize=12, va='bottom', fontweight='bold')

plt.tight_layout()
plt.savefig(f'{output_folder}figure_1.svg', transparent=True)
plt.savefig(f'{output_folder}figure_1.png', transparent=True)
plt.show()
    