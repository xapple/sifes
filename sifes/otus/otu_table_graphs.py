# Futures #
from __future__ import division

# Built-in modules #

# Internal modules #
from plumbing.graphs import Graph
from plumbing.common import split_thousands

# Third party modules #
import matplotlib
from matplotlib import pyplot

# Constants #
__all__ = ['OtuSizesDist', 'OtuSumsPerSample', 'SampleSumsPerOtu', 'CumulativePresence']

################################################################################
class OtuSizesDist(Graph):
    """Distribution of OTU cluster sizes in log-log."""

    short_name = 'otu_sizes_dist'
    x_grid     = True
    y_grid     = True
    x_scale    = 'symlog'
    y_scale    = 'log'
    x_label    = 'Number of sequences in an OTU'
    y_label    = 'Number of OTUs with that many sequences in them'

    def plot(self, **kwargs):
        # Sum by column and count frequencies #
        distrib = self.parent.otu_table.sum().value_counts()
        x = distrib.keys()
        y = distrib.values
        # Make scatter #
        fig = pyplot.figure()
        axes = fig.add_subplot(111)
        axes.plot(x, y, 'ro')
        axes.set_title('Distribution of sizes for %s OTUs' % split_thousands(sum(y)))
        # Add annotations #
        for i in range(min(5,len(x))):
            pyplot.annotate("%i: %s" % (x[i], split_thousands(y[i])), size=13, xy = (x[i], y[i]), xytext = (10, 0),
                            textcoords = 'offset points', ha = 'left', va = 'center',
                            bbox = dict(boxstyle = 'round,pad=0.2', fc = 'yellow', alpha = 0.3))
        # Save it #
        self.save_plot(fig, axes, **kwargs)
        pyplot.close(fig)

################################################################################
class OtuSumsPerSample(Graph):

    short_name = 'otu_sums_graph'
    title      = 'Histogram of OTU appearance sums per sample'
    x_label    = 'Number of OTUs present (non-null) in a sample'
    y_label    = 'Number of samples with that many OTUs in them'

    def plot(self, **kwargs):
        # Sum by row and count frequencies #
        self.frame = self.parent.otu_table.astype(bool).sum(axis=1)
        # Make hist #
        fig = pyplot.figure()
        axes = self.frame.hist(color='gray', bins=40)
        # Save it #
        self.save_plot(fig, axes, **kwargs)
        pyplot.close(fig)

################################################################################
class SampleSumsPerOtu(Graph):

    short_name = 'sample_sums_graph'
    title      = 'Histogram of OTU appearance sums per OTU'
    y_label    = 'Number of OTUs that appear in these many samples'

    def plot(self, **kwargs):
        # Sum by row and count frequencies #
        self.frame = self.parent.otu_table.astype(bool).sum(axis=0).value_counts().sort_index()
        # Make hist #
        fig = pyplot.figure()
        axes = self.frame.plot(kind='bar', color='gray')
        axes.set_xlabel('Number of samples an OTU appears in (max. %i)' % self.parent.otu_table.shape[0])
        # Save it #
        self.save_plot(fig, axes, **kwargs)
        pyplot.close(fig)

################################################################################
class CumulativePresence(Graph):
    """Cumulative graph of cluster presence in samples. This means something such as:
    - 0% of OTUs appear in 100% of the samples,
    - 10% of OTUs appear in 90% of the samples,
    - 90% of OTUs appear in 1% of the samples."""

    short_name = 'cumulative_presence'

    def plot(self, **kwargs):
        # Number of samples #
        num_of_otus    = self.parent.otu_table.shape[1]
        num_of_samples = self.parent.otu_table.shape[0]
        samples_index  = list(reversed(range(1, num_of_samples+1)))
        # Get value frequencies #
        counts = self.parent.otu_table.astype(bool).sum(axis=0).value_counts()
        # Add missing values #
        for n in samples_index:
            if n not in counts:
                counts = counts.set_value(n,0)
        # Sort it #
        counts = counts.sort_index(ascending=False)
        # Cumulative sum #
        self.y = list(counts.cumsum())
        # Percentage of samples #
        self.x = [n/num_of_samples for n in samples_index]
        # Make step plot #
        fig = pyplot.figure()
        axes = fig.add_subplot(111)
        axes.step(self.x, self.y, fillstyle='bottom')
        axes.set_title('Cumulative graph of OTU presence in samples for %s OTUs' % num_of_otus)
        axes.set_xlabel('Fraction of samples (100%% equates %i samples)' % num_of_samples)
        axes.set_ylabel('Number of OTUs that appear in that fraction of samples or more')
        axes.invert_xaxis()
        axes.set_xticks([min(self.x) + (max(self.x)-min(self.x))* n / 9 for n in range(10)])
        axes.set_yscale('log')
        axes.xaxis.grid(True)
        # Set percentage #
        percentage = lambda x, pos: '%1.0f%%' % (x*100.0)
        axes.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(percentage))
        # Save it #
        self.save_plot(fig, axes, **kwargs)
        pyplot.close(fig)
