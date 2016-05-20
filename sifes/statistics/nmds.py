# Internal modules #
from plumbing.dataframes import r_matrix_to_dataframe
from plumbing.graphs     import Graph

# Third party modules #
from matplotlib import pyplot

################################################################################
class GraphNMDS(Graph):
    """Non-metric dimensional scaling plot."""

    short_name = 'nmds_horn'
    bottom     = 0.03
    top        = 0.97

    def plot(self, **kwargs):
        # The otu table path #
        self.tsv = self.parent.otu_table.p.otu_table_flat
        # Run via R #
        self.run_via_R()
        # Data #
        x     = self.coords['NMDS1'].values
        y     = self.coords['NMDS2'].values
        names = self.coords['NMDS1'].keys()
        # Make scatter #
        fig  = pyplot.figure()
        axes = fig.add_subplot(111)
        axes.plot(x, y, 'ro')
        axes.set_title('Non-Metric Multidimensional scaling')
        axes.set_xlabel('Dimension 1')
        axes.set_ylabel('Dimension 2')
        # Add annotations #
        for i in range(len(names)):
            pyplot.annotate(names[i], size=9, xy = (x[i], y[i]), xytext = (10, 0),
                            textcoords = 'offset points', ha = 'left', va = 'center',
                            bbox = dict(boxstyle = 'round,pad=0.2', fc = 'yellow', alpha = 0.3))
        # Save it #
        self.save_plot(fig, axes, **kwargs)
        pyplot.close(fig)

    def run_via_R(self):
        # Module on demand #
        from rpy2 import robjects as ro
        # Load dataframe #
        ro.r("library(vegan)")
        ro.r("table = read.table('%s', sep='\t', header=TRUE, row.names='X')" % (self.tsv))
        # Run computation #
        ro.r("nmds = metaMDS(table, distance='horn', trymax=200)")
        # Extract result #
        ro.r("coord = scores(nmds)")
        ro.r("loadings = nmds$species")
        # Retrieve values #
        self.coords = r_matrix_to_dataframe(ro.r.coord)