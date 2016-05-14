"""
Check here for the list of estimators:
http://scikit-bio.org/docs/latest/generated/skbio.diversity.alpha.html
"""

# Built-in modules #
from functools import partial

# Internal modules #
from plumbing.graphs import Graph

# Third party modules #
from skbio.diversity import alpha_diversity
from skbio.stats     import subsample_counts as subsample
from matplotlib      import pyplot
from numpy           import linspace

###############################################################################
class AlphaDiversity(object):
    """All the alpha diversity estimates to be calculated on an OTU table."""

    def __init__(self, parent):
        self.sample, self.parent = parent, parent
        self.chao1   = Chao1(self, base_dir=self.sample.p.graphs_dir)
        self.ace     = Ace(self, base_dir=self.sample.p.graphs_dir)
        self.shannon = Shannon(self, base_dir=self.sample.p.graphs_dir)
        self.simpson = Simpson(self, base_dir=self.sample.p.graphs_dir)
        self.graphs  = [self.chao1, self.ace, self.shannon, self.simpson]

    def plot_all_graphs(self):
        for graph in self.graphs: graph.plot()

###############################################################################
class AlphaDiversityGraph(Graph):
    bottom = 0.10
    right  = 0.96
    sep    = 'x'

    @property
    def x(self):
        return map(int, linspace(0, sum(self.parent.sample.counts), 600))

    @property
    def y(self):
        try: return [self.div_fn(subsample(self.parent.sample.counts, k)) for k in self.x]
        except ValueError: return [0 for k in self.x]

    def plot(self, **kwargs):
        # Plot #
        fig  = pyplot.figure()
        axes = fig.add_subplot(111)
        axes.plot(self.x, self.y, 'ro')
        # Labels #
        axes.set_title("Rarefaction curve of the " + self.title + " diversity estimate")
        axes.set_xlabel('Sequences rarefied down to this many')
        axes.set_ylabel(self.title + " diversity estimate")
        # Options #
        axes.yaxis.grid(True)
        axes.set_xlim(0, axes.get_xlim()[1])
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
        total = sum(self.parent.sample.counts)
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