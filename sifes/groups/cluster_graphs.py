# Futures #
from __future__ import division

# Built-in modules #
from functools import partial

# Internal modules #
from sifes.location.map_figure     import MapFigure
from sifes.metadata.correspondence import reverse_corr

# First party modules #
from plumbing.graphs import Graph
from plumbing.common import uniquify_list

# Third party modules #
import matplotlib, numpy
from matplotlib import pyplot
from mpl_toolkits.axes_grid.anchored_artists import AnchoredText
from skbio.diversity import alpha_diversity  as alphadiv
from skbio.stats     import subsample_counts as subsample
from scipy.stats     import linregress

# Constants #
__all__ = ['Chao1', 'Ace', 'Shannon', 'Simpson']

###############################################################################
class ClusterLocationMap(MapFigure):
    """Map of several sample locations."""
    short_name = 'cluster_location_map'

    def __init__(self, group_name, samples, *args, **kwargs):
        self.group_name = group_name
        self.samples    = samples
        super(ClusterLocationMap, self).__init__(*args, **kwargs)

    def plot(self, **kwargs):
        # Colors based on custom attributes #
        colors = list(set(s.attribute for s in self.samples))
        # Loop #
        for s in self.samples:
            self.add_marker(s.latitude, s.longitude, color=colors.index(s.attribute))
        # Do it #
        self.save_map(**kwargs)

###############################################################################
class DiversityRegression(Graph):
    """Group samples and plot their alpha diversity along a metric."""

    # Parameters #
    bottom = 0.1

    # Options #
    custom_metadata  = 'depth'

    def plot(self, **kwargs):
        # The custom sample grouping #
        groups = uniquify_list([s.grouping for s in self.parent.samples])
        # Make as many subplots as there are groups #
        self.fig, self.axes = pyplot.subplots(1, len(groups), sharex=True, sharey=True)
        # Make each subplot #
        for axe, group in zip(self.axes, groups): self.subplot(axe, group)
        # Y axis label #
        self.axes[0].set_ylabel(self.description)
        # X axis label #
        for a in self.axes: a.set_xlabel(reverse_corr.get(self.custom_metadata, self.custom_metadata))
        # Save it #
        self.save_plot(self.fig, self.axes, **kwargs)

    def subplot(self, axe, group):
        # Basic data #
        samples = [s for s in self.parent.samples if s.grouping == group]
        x = numpy.array([float(s.info.get(self.custom_metadata))    for s in samples])
        y = numpy.array([float(self.div_fn(s, self.parent.down_to)) for s in samples])
        axe.plot(x, y, 'ko')
        axe.title.set_text(group)
        # Regression #
        slope, intercept, r_value, p_value, std_err = linregress(x,y)
        prediction_y       = intercept + slope * x
        prediction_error   = y - prediction_y
        degrees_of_freedom = len(x) - 2
        residual_std_error = numpy.sqrt(numpy.sum(prediction_error**2) / degrees_of_freedom)
        # Plot line #
        axe.plot(x, prediction_y, 'r-')
        # Add regression result #
        matplotlib.rc('text', usetex=True)
        text  = "\\centerline{\\textbf{Least-squares reg.}}\n"
        text += "{\\raggedright \\textit{Residual stderr}: \\texttt{%.2f}}\n"
        text += "{\\raggedright \\textit{R\\textsuperscript{2}:} \\texttt{%.2f}}\n"
        text += "{\\raggedright \\textit{Slope}: \\texttt{%.2E}}"
        text  = text % (residual_std_error, r_value**2, slope)
        anchor = AnchoredText(text, prop=dict(size=11), frameon=True, loc=4)
        anchor.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
        axe.add_artist(anchor)

    def div_fn(self, s, k):
        return alphadiv(self.diversity_metric, subsample(s.otu_counts, k))

###############################################################################
class Chao1(DiversityRegression):
    diversity_metric = 'chao1'
    short_name       = 'diversity_reg_chao1'
    description      = "Chao1 (bias-corrected version)"

class Ace(DiversityRegression):
    diversity_metric = 'ace'
    short_name       = 'diversity_reg_ace'
    description      = "Ace (Abundance-based Coverage Estimator)"

class Shannon(DiversityRegression):
    diversity_metric = 'shannon'
    short_name       = 'diversity_reg_shannon'
    description      = "Shannon (entropy of counts H, in bits)"

class Simpson(DiversityRegression):
    diversity_metric = 'simpson'
    short_name       = 'diversity_reg_simpson'
    description      = "Simpson (1 - dominance)"
