# Font files directory #
# /home/lucas/.local/lib/python2.7/site-packages/matplotlib/mpl-data/fonts/ttf/

# Built-in modules #
from collections import OrderedDict

# Internal modules #
from plumbing.graphs import Graph
from plumbing.common import split_thousands

# Third party modules #
import numpy, matplotlib
from matplotlib import pyplot

# Constants #
__all__ = ['MistagHeatmap']

################################################################################
class MistagHeatmap(Graph):
    """A heatmap showing mistags."""
    short_name = 'mistag_heatmap'
    width  = 8.0
    height = 8.0

    top    = 0.8
    bottom = 0.05
    left   = 0.2
    right  = 0.72

    x_label = "Forward barcodes"
    y_label = "Reverse barcodes"

    def plot(self, **kwargs):
        # Attributes #
        samples = self.parent.samples
        # Names #
        fwd_barcodes = OrderedDict((s.info['forward_mid'], s.info['forward_num']) for s in samples)
        rev_barcodes = OrderedDict((s.info['reverse_mid'], s.info['reverse_num']) for s in samples)
        # Data #
        df = self.parent.read_counts
        # Sample counts #
        df = df.replace({s.short_name: len(s) for s in samples})
        df = df.astype(int)
        # Masks #
        mistag_mask = (df == df) # All True
        sample_mask = (df != df) # All False
        for s in samples: mistag_mask[s.info['forward_mid']][s.info['reverse_mid']] = False
        for s in samples: sample_mask[s.info['forward_mid']][s.info['reverse_mid']] = True
        # Colors #
        mistag_vals = df.where(mistag_mask)
        sample_vals = df.where(sample_mask)
        # Figure #
        fig = pyplot.figure()
        axes = pyplot.gca()
        # Plot #
        mistag_mesh = axes.matshow(mistag_vals, cmap='winter' )
        sample_mesh = axes.matshow(sample_vals, cmap='spring')
        # Titles #
        axes.set_xlabel(self.x_label)
        axes.set_ylabel(self.y_label)
        axes.xaxis.set_label_position('top')
        # Labels #
        axes.set(xticks = numpy.arange(df.shape[1]), # fwd = df.columns
                 yticks = numpy.arange(df.shape[0])) # rev = df.rows
        axes.set_xticklabels(("%s (%s)" % (v[0],k) for k,v in fwd_barcodes.items()),
                             rotation=90, fontname='Helvetica Neue', fontweight='normal')
        axes.set_yticklabels(("(%s) %s" % (k,v[0]) for k,v in fwd_barcodes.items()),
                             rotation=0,  fontname='Helvetica Neue', fontweight='normal')
        axes.xaxis.tick_top()
        axes.yaxis.tick_left()
        axes.tick_params(direction='out')
        # Add dual colorbars #
        mistags_cbar = fig.colorbar(mistag_mesh, cax=fig.add_axes([0.75, 0.05, 0.04, 0.83]))
        samples_cbar = fig.colorbar(sample_mesh, cax=fig.add_axes([0.86, 0.05, 0.04, 0.83]))
        # Split thousands #
        separate = lambda x,pos: split_thousands(x)
        mistags_cbar.ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(separate))
        samples_cbar.ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(separate))
        # Smaller labels #
        #for tick in mistags_cbar.ax.yaxis.get_major_ticks(): tick.label.set_fontsize(5)
        #for tick in samples_cbar.ax.yaxis.get_major_ticks(): tick.label.set_fontsize(5)
        mistags_cbar.ax.set_yticklabels(mistags_cbar.ax.get_yticklabels(), fontsize='smaller')
        samples_cbar.ax.set_yticklabels(samples_cbar.ax.get_yticklabels(), fontsize='smaller')
        # Add legend #
        mistags_cbar.ax.text(0.55, 0.8, 'Mistags', rotation=90, ha='center', va='center',
                     transform=mistags_cbar.ax.transAxes, color='black')
        samples_cbar.ax.text(0.55, 0.8, 'Real samples', rotation=90, ha='center', va='center',
                     transform=samples_cbar.ax.transAxes, color='black')
        # Save it #
        self.save_plot(fig, axes, **kwargs)
        pyplot.close(fig)
        # For convenience #
        return self
