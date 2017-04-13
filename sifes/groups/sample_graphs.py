"""
Check here for the list of estimators:
http://scikit-bio.org/docs/latest/generated/skbio.diversity.alpha.html
"""

# Built-in modules #
from functools import partial

# Internal modules #
from sifes.location.map_figure import MapFigure

# First party modules #
from plumbing.graphs import Graph

# Third party modules #
from skbio.diversity import alpha_diversity
from skbio.stats     import subsample_counts as subsample
from matplotlib      import pyplot
from numpy           import linspace

# Constants #
__all__ = ['Chao1',
           'Ace',
           'Shannon',
           'Simpson',
           'LocationMap']

###############################################################################
class AlphaDiversityGraph(Graph):
    bottom = 0.10
    right  = 0.96
    sep    = 'x'
    y_grid = True
    x_min  = 0

    @property
    def x(self):
        return map(int, linspace(0, sum(self.parent.otu_counts), 600))

    @property
    def y(self):
        try: return [self.div_fn(subsample(self.parent.otu_counts, k)) for k in self.x]
        except ValueError: return [0 for k in self.x]

    def plot(self, **kwargs):
        # Plot #
        fig  = pyplot.figure()
        axes = fig.add_subplot(111)
        axes.plot(self.x, self.y, 'ro')
        # Labels #
        axes.set_title("Rarefaction curve of the '" + self.title + "' diversity estimate")
        axes.set_xlabel('Sequences rarefied down to this many')
        axes.set_ylabel(self.title + " diversity estimate")
        # Save #
        self.save_plot(fig, axes, **kwargs)
        pyplot.close(fig)

###############################################################################
class Chao1(AlphaDiversityGraph):
    """Will graph the Chao1 rarefaction curve."""
    short_name = 'chao1'
    title = "Chao1 (bias-corrected version)"

    @property
    def div_fn(self): return partial(alpha_diversity, 'chao1')

###############################################################################
class Ace(AlphaDiversityGraph):
    """Will graph the Ace rarefaction curve."""
    short_name = 'ace'
    title = "Ace (Abundance-based Coverage Estimator)"

    @property
    def div_fn(self): return partial(alpha_diversity, 'ace')

    @property
    def x(self):
        total = sum(self.parent.otu_counts)
        return map(int,linspace(total/2, total, 600))

###############################################################################
class Shannon(AlphaDiversityGraph):
    """Will graph the Shannon rarefaction curve."""
    short_name = 'shannon'
    title = "Shannon (entropy of counts H, in bits)"

    @property
    def div_fn(self): return partial(alpha_diversity, 'shannon')

###############################################################################
class Simpson(AlphaDiversityGraph):
    """Will graph the Simpson rarefaction curve."""
    short_name = 'simpson'
    title = "Simpson (1 - dominance)"

    @property
    def div_fn(self): return partial(alpha_diversity, 'simpson')

###############################################################################
class LocationMap(MapFigure):
    """Map of sample location."""
    short_name = 'location_map'

    def plot(self, **kwargs):
        self.add_marker(self.parent.latitude, self.parent.longitude)
        self.save_map(**kwargs)