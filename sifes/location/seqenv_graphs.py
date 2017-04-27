# Futures #
from __future__ import division

# Built-in modules #
import warnings

# Internal modules #
from sifes.metadata.correspondence import reverse_corr

# First party modules #
from plumbing.graphs import Graph

# Third party modules #
import pandas, numpy
from matplotlib import pyplot
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap

# Modules using deprecated interfaces #
with warnings.catch_warnings():
    warnings.filterwarnings("ignore")
    import seaborn

# Constants #
__all__ = ['SeqenvHeatmap']

# Colors #
row_colors = [(0.510530904690042, 0.6614299289084904, 0.1930849118538962),
              (0.775731904186273, 0.5784925270759935, 0.1947556653855187),
              (0.216629789230736, 0.6676586160122123, 0.7318695594345369),
              (0.958705008049440, 0.3662259565791745, 0.9231469575614251)]

################################################################################
class SeqenvHeatmap(Graph):
    """A cool clustered heatmap."""

    # Parameters #
    short_name = 'seqenv_heatmap'
    width  = 14.0
    height = 10.0
    left   = 0.00
    right  = 0.90
    bottom = 0.2

    # Options #
    custom_metadata = None
    max_columns     = 24

    def plot(self, **kwargs):
        # Data #
        self.df = self.parent.samples_to_names.copy()
        # Filter #
        best_columns = self.df.sum(axis=1).sort_values(ascending=False).index[0:self.max_columns]
        self.df = self.df.reindex_axis(best_columns)
        # Transform #
        self.df = self.df.transpose()
        self.df = self.df.applymap(numpy.sqrt)
        # Row colors #
        custom_groups   = list(set(s.grouping for s in self.parent.samples))
        name_to_color   = lambda n: row_colors[custom_groups.index(self.parent.cluster[n].grouping)]
        self.row_colors = pandas.Series({n: name_to_color(n) for n,r in self.df.iterrows()})
        self.patches    = [mpatches.Patch(color=row_colors[i], label=g) for i,g in enumerate(custom_groups)]
        # Row right labels #
        if self.custom_metadata:
            name_to_md        = lambda n: self.parent.cluster[n].info.get(self.custom_metadata)
            self.custom_md_fn = lambda n: "" if pandas.isnull(name_to_md(n)) else name_to_md(n)
        # Custom color map #
        start_color   = (1.0,   1.0,   0.15)
        sentinel_cl   = (0.95,  0.5,   0.15)
        mid_color     = (1.0,   0.0,   0.0)
        end_colors    = (0.09,  0.09,  1.0)
        custom_colors = [end_colors, mid_color, sentinel_cl, start_color]
        # Different color maps options #
        self.color_map = pyplot.get_cmap("viridis")
        self.color_map = LinearSegmentedColormap.from_list('r_cmap', colors=custom_colors)
        # Plot #
        self.clustergrid = seaborn.clustermap(self.df, cmap        = self.color_map,
                                                       linewidths  = .03,
                                                       row_colors  = self.row_colors,
                                                       row_cluster = False)
        # Objects #
        self.fig  = pyplot.gcf()
        self.axes = self.clustergrid.ax_heatmap
        # Rotate labels #
        pyplot.setp(self.axes.yaxis.get_majorticklabels(), rotation=0)
        pyplot.setp(self.axes.xaxis.get_majorticklabels(), rotation=90)
        # Copy left labels to row colors axis #
        self.clustergrid.ax_row_colors.set_yticks(self.axes.get_yticks())
        # Rename ticks to the left (samples names) #
        labels = reversed([n for n,r in self.df.iterrows()])
        self.clustergrid.ax_row_colors.set_yticklabels(labels)
        # Rename ticks to the right (custom metadata param) #
        if self.custom_metadata:
            labels = reversed([self.custom_md_fn(n) for n,r in self.df.iterrows()])
            self.axes.set_yticklabels(labels)
        # Y labels names #
        self.clustergrid.ax_row_colors.set_ylabel("Sample names")
        if self.custom_metadata: self.axes.set_ylabel(reverse_corr.get(self.custom_metadata, self.custom_metadata))
        # Remove the 'None' x label #
        self.clustergrid.ax_row_colors.set_xticklabels([])
        # Add legend #
        pyplot.legend(handles=self.patches, bbox_to_anchor=(1, 1), bbox_transform=self.fig.transFigure)
        # Change font of sample names #
        ticks = self.clustergrid.ax_row_colors.get_yticklabels()
        self.clustergrid.ax_row_colors.set_yticklabels(ticks, fontname='Menlo', fontweight='bold', size='medium')
        # Save it #
        self.save_plot(self.fig, self.axes, **kwargs)