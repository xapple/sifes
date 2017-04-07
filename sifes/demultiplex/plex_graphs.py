# Built-in modules #

# Internal modules #
from plumbing.graphs import Graph

# Third party modules #
from matplotlib import pyplot

# Constants #
__all__ = ['ReadCountsHeatmap']

################################################################################
class ReadCountsHeatmap(Graph):
    """A heatmap showing mistags."""
    short_name = 'mistag_heatmap'

    def plot(self, **kwargs):
        # Data #
        x = 1
        y = 1
        # Plot #
        fig = pyplot.figure()
        pyplot.plot(x, y)
        axes = pyplot.gca()
        # Save it #
        self.save_plot(fig, axes, **kwargs)
        pyplot.close(fig)
        # For convenience #
        return self
