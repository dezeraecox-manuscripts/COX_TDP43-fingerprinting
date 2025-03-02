"""General purpose visualisation functions"""

import os
import re
import subprocess
from math import pi
from textwrap import wrap

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import seaborn as sns
from loguru import logger
from matplotlib.patches import Circle, RegularPolygon
from matplotlib.path import Path
from matplotlib.projections import register_projection
from matplotlib.projections.polar import PolarAxes
from matplotlib.spines import Spine
from matplotlib.transforms import Affine2D
from scipy.interpolate import interp1d
from skimage import io
from statannotations.Annotator import Annotator

logger.info('Import OK')

def bardotplot(data, xcol, ycol, order, hue=None, hue_order=None, scat_hue=None, scat_hue_order=None, palette=False, xlabel='', ylabel=False, pairs=False, correction=None, xticks=None, groups=None, group_label_y=-0.18, group_line_y=-0.05, ax=None, legend='', dot_size=5, cap_size=0.2, cap_width=2):
    if ax is None:
        fig, ax = plt.subplots()
    sns.barplot(
        data=data,
        x=xcol,
        y=ycol,
        hue=hue,
        palette=palette,
        capsize=cap_size,
        errwidth=cap_width,
        ax=ax,
        dodge=True,
        order=order,
        hue_order=hue_order,
        edgecolor='white'
    )
    sns.stripplot(
        data=data,
        x=xcol,
        y=ycol,
        hue=scat_hue,
        palette=palette,
        ax=ax,
        edgecolor='#fff',
        linewidth=1,
        s=dot_size,
        order=order,
        hue_order=scat_hue_order,
        dodge=True,
    )

    if pairs:
        annotator = Annotator(
            ax=ax, pairs=pairs, data=data, x=xcol, y=ycol, order=order, hue=hue, hue_order=hue_order)
        annotator.configure(test='t-test_ind', text_format='star',
                            loc='inside', comparisons_correction=correction, line_width=0.5)
        annotator.apply_and_annotate()

    ax.set(ylabel=ylabel)

    ax.set_xlabel(xlabel)
    if xticks:
        ax.set_xticks(xticks)
        ax.set_xticklabels(hue_order*len(order))
    if groups:
        for group_label, (x0, x1, x2) in groups.items():
            ax.annotate(group_label, xy=(x0, group_label_y),
                        xycoords='data', ha='center', annotation_clip=False)
            trans = ax.get_xaxis_transform()
            ax.plot([x1, x2], [group_line_y, group_line_y],
                    color="black", transform=trans, clip_on=False)

    if legend == '':
        ax.legend('', frameon=False)
    else:    
        handles, labels = plt.gca().get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        ax.legend(by_label.values(), by_label.keys())
    
    return ax

# interpolate ecdf for combined visualisation

def plot_interpolated_ecdf(fitted_ecdfs, ycol, huecol, palette, ax=None, orientation=None, linestyle=None):
    """Generate fitted cumulative distribution

    Args:
        fitted_ecdfs (dataframe): Dataframe containing interpolated cumulative distributions
        ycol (str): Value of the cumulative distribution
        huecol (str): Column on which the hue is based
        palette (dictionary): Dictionary containing the palette
        ax (str, optional): Subpanel number in multipanel figures. Defaults to None.
        orientation (str, optional): Orientation of the plot, typically h. Defaults to None.

    Returns:
        _type_: _description_
    """

    if not ax:
        fig, ax = plt.subplots()

    if orientation == 'h':
        means = fitted_ecdfs[fitted_ecdfs['type'] == 'interpolated'].groupby(
            [huecol, 'ecdf']).agg(['mean', 'std'])[ycol].reset_index()
        means.columns = [huecol, 'ecdf', 'mean', 'std']
        means['pos_err'] = means['mean'] + means['std']
        means['neg_err'] = means['mean'] - means['std']

        for hue, data in means.groupby([huecol]):

            ax.plot(
                data['mean'],
                data['ecdf'],
                color=palette[hue],
                label=hue,
                linestyle=linestyle,
            )
            ax.fill_betweenx(
                y=data['ecdf'].tolist(),
                x1=(data['neg_err']).tolist(),
                x2=(data['pos_err']).tolist(),
                color=palette[hue],
                alpha=0.3
            )

    else:
        sns.lineplot(
            data=fitted_ecdfs,
            y=ycol,
            x='ecdf',
            hue=huecol,
            palette=palette,
            ci='sd',
            ax=ax,
            linestyle=linestyle
        )

    return fitted_ecdfs, ax


def plot_sankey(flow_data, palette, title=None, output_path=None, h=500, w=500):
    data = []
    y=0
    for category, df in flow_data.groupby('Category'):
        indices = [category, "long", "short", "fibril", 'globular', 'dense', 'sparse']
        indices = dict(zip(indices, np.arange(len(indices))))
    
        df['target_A'] = df['Category'].map(indices)
        df['target_B'] = df['cat_smoothed_length'].map(indices)
        df['target_C'] = df['cat_eccentricity'].map(indices)
        df['target_D'] = df['cat_density'].map(indices)
        df['color'] = [palette[f'{cat}_{category}'] for cat in df['cat_smoothed_length']]

        links = []
        targets = ['target_A', 'target_B', 'target_C', 'target_D', ]
        for x in range(len(targets)-1):
            links.extend(df[[targets[x], targets[x+1], 'proportion', 'color']].values)
        links = pd.DataFrame(links, columns=['source', 'target', 'proportion', 'color'])
        links = links.groupby(['color', 'source', 'target',]).sum().reset_index()

        ypos = dict(links.groupby('target').sum()['proportion'].reset_index().values)

        data.append(go.Sankey(
            domain={
                'x': [0, 1],
                'y': [y/3, (1+y)/3-0.05],
            },
            node = dict(
            pad = 15,
            thickness = 20,
            line = dict(color = "black", width = 0.5),
            # label = [category, "Long", "Short", "Fibrillar", 'Globular', 'Dense', 'Sparse'],
            color = palette[category],
            x = [0.001, 0.33, 0.33, 0.66, 0.66, 0.999, 0.999],
            y = [0.5, ypos[1]/2, (ypos[1]+ypos[2]/2), ypos[3]/2, (ypos[3]+ypos[4]/2), ypos[5]/2, (ypos[5]+ypos[6]/2)]
            ),
            link = dict(
            source = links['source'].tolist(), 
            target = links['target'].tolist(),
            value = links['proportion'].tolist(),
            color = links['color'].tolist(),
        )))
    
        y+=1 
    
    fig = go.Figure(
    data=data,
    )

    fig.update_layout(title=title, font_size=8, height=h, width=w,)
    if output_path:
        fig.write_image(f'{output_path}')
    fig.show()


def sankey_converter(input_path, x0=60.750, w=380.500, y0=110.750, h=314.500):
    subprocess.run(f"inkscape {input_path} --export-filename={input_path.replace('.svg', '.png')} --export-dpi=500 --export-area={x0-0.5}:{y0-0.5}:{x0+w+0.5}:{y0+h+0.5}")
    return io.imread(input_path.replace('.svg', '.png'))


# Visualise volcano plot
def plot_volcano(df, xcol, ycol, palette, hue_col=False, ax=False, xthresh=False, ythresh=False, xlim=False, ylim=False, alpha=1):
    if not ax:
        fig, ax = plt.subplots()
    sns.scatterplot(data=df, x=xcol,
                    y=ycol, palette=palette, hue=hue_col, ax=ax, zorder=1000, alpha=alpha)
    if xlim:
        ax.set_ylim(*ylim)
    if ylim:
        ax.set_xlim(*xlim)
    if xthresh:
        for thresh in xthresh:
            ax.axvline(thresh, linestyle='--', color='lightgrey', zorder=-10)
    if ythresh:
        ax.axhline(1.3, linestyle='--', color='lightgrey',  zorder=-100)
    return ax


def plot_lda_scatter(df, hue, palette, style, ax=None, s=300):
    if not ax:
        fig, ax = plt.subplots(figsize=(6, 5.5))
    sns.scatterplot(
        data=df,
        x='dim1',
        y='dim2',
        hue=hue,
        style=style,
        palette=palette,
        s=s,
        ax=ax
    )

    # ax.legend(bbox_to_anchor=(1.0, 1.0), ncol=3)
    ax.set_xlabel('Dimension 1')
    ax.set_ylabel('Dimension 2')
    
    return ax



def radar_factory(num_vars, frame='circle'):
    """
    Create a radar chart with `num_vars` axes.

    This function creates a RadarAxes projection and registers it.

    Parameters
    ----------
    num_vars : int
        Number of variables for radar chart.
    frame : {'circle', 'polygon'}
        Shape of frame surrounding axes.

    """
    # calculate evenly-spaced axis angles
    theta = np.linspace(0, 2*np.pi, num_vars, endpoint=False)

    class RadarTransform(PolarAxes.PolarTransform):

        def transform_path_non_affine(self, path):
            # Paths with non-unit interpolation steps correspond to gridlines,
            # in which case we force interpolation (to defeat PolarTransform's
            # autoconversion to circular arcs).
            if path._interpolation_steps > 1:
                path = path.interpolated(num_vars)
            return Path(self.transform(path.vertices), path.codes)

    class RadarAxes(PolarAxes):

        name = 'radar'
        PolarTransform = RadarTransform

        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            # rotate plot such that the first axis is at the top
            self.set_theta_zero_location('N')

        def fill(self, *args, closed=True, **kwargs):
            """Override fill so that line is closed by default"""
            return super().fill(closed=closed, *args, **kwargs)

        def plot(self, *args, **kwargs):
            """Override plot so that line is closed by default"""
            lines = super().plot(*args, **kwargs)
            for line in lines:
                self._close_line(line)

        def _close_line(self, line):
            x, y = line.get_data()
            # FIXME: markers at x[0], y[0] get doubled-up
            if x[0] != x[-1]:
                x = np.append(x, x[0])
                y = np.append(y, y[0])
                line.set_data(x, y)

        def set_varlabels(self, labels):
            self.set_thetagrids(np.degrees(theta), labels)

        def _gen_axes_patch(self):
            # The Axes patch must be centered at (0.5, 0.5) and of radius 0.5
            # in axes coordinates.
            if frame == 'circle':
                return Circle((0.5, 0.5), 0.5)
            elif frame == 'polygon':
                return RegularPolygon((0.5, 0.5), num_vars,
                                      radius=.5, edgecolor="k")
            else:
                raise ValueError("Unknown value for 'frame': %s" % frame)

        def _gen_axes_spines(self):
            if frame == 'circle':
                return super()._gen_axes_spines()
            elif frame == 'polygon':
                # spine_type must be 'left'/'right'/'top'/'bottom'/'circle'.
                spine = Spine(axes=self,
                              spine_type='circle',
                              path=Path.unit_regular_polygon(num_vars))
                # unit_regular_polygon gives a polygon of radius 1 centered at
                # (0, 0) but we want a polygon of radius 0.5 centered at (0.5,
                # 0.5) in axes coordinates.
                spine.set_transform(Affine2D().scale(.5).translate(.5, .5)
                                    + self.transAxes)
                return {'polar': spine}
            else:
                raise ValueError("Unknown value for 'frame': %s" % frame)

    register_projection(RadarAxes)
    return theta


def plot_radar(dataframe, label_col, value_col, hue_col, palette, err_bars=True, label_dict=None, ax_labels=None, ax=None):
    
    categories = dataframe[label_col].unique().tolist()
    N = len(categories)
    theta = radar_factory(N, frame='polygon')

    if not ax:
        fig, ax = plt.subplots(figsize=(9, 9),subplot_kw=dict(projection='radar'))

    
    for category, df in dataframe.groupby(hue_col):
        
        color = palette[category]
        yvals = dict(df.groupby(label_col).mean()[value_col].reset_index().values)
        err = dict(df.groupby(label_col).std()[value_col].reset_index().values)
        vals = np.array(list(yvals.values()))
        ax.plot(theta, vals, color=color, label=category)
        
        if err_bars:
        
            y_err = np.array([err[key] for key in yvals])
            err_up =  vals + y_err
            err_down =  vals - y_err
            
            ax.fill_between( [*theta, theta[0]], [*err_up, err_up[0]], [*err_down, err_down[0]], facecolor=color, alpha=0.25)
    if ax_labels:
        ax.set_rgrids(ax_labels, fontsize=5) # axis labels
    if label_dict:
        categories = [label_dict[label] for label in yvals]
    else:
        categories = list(yvals.keys())
    ax.set_varlabels(categories)
        
    return ax    

