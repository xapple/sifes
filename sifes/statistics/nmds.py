# Internal modules #
from plumbing.dataframes import r_matrix_to_dataframe
from plumbing.graphs     import Graph

# Third party modules #
from matplotlib import pyplot

################################################################################
class GraphNMDS(Graph):
    """Non-metric dimensional scaling plot."""

    short_name = 'nmds_horn'

    def plot(self, **kwargs):
        # The otu table path #
        self.tsv = self.parent.otu_table.p.otu_table_flat
        # Run via R #
        self.run_via_R()
        # Data #
        x      = self.coords['NMDS1'].values
        y      = self.coords['NMDS2'].values
        names  = self.coords['NMDS1'].keys()
        labels = [self.parent[k].info.get('short_label', '').replace(' ','') for k in names]
        # Make scatter #
        fig  = pyplot.figure()
        axes = fig.add_subplot(111)
        axes.plot(x, y, 'ro')
        axes.set_title('Non-Metric Multidimensional scaling')
        axes.set_xlabel('Dimension 1')
        axes.set_ylabel('Dimension 2')
        # Add annotations #
        for i in range(len(names)):
            pyplot.annotate(labels[i] or names[i], size=9, xy = (x[i], y[i]), xytext = (10, 0),
                            textcoords = 'offset points', ha = 'left', va = 'center',
                            bbox = dict(boxstyle = 'round,pad=0.2', fc = 'yellow', alpha = 0.3))
        # Save it #
        self.save_plot(fig, axes, **kwargs)
        pyplot.close(fig)

    def run_via_R(self):
        # Import needed rpy2 objects #
        from rpy2 import robjects   as ro
        from rpy2 import rinterface as ri
        # Record output and don't display it #
        self.stdout = []
        self.stderr = []
        def add_to_stdout(line): self.stdout.append(line)
        def add_to_stderr(line): self.stderr.append(line)
        ri.set_writeconsole_regular(add_to_stdout)
        ri.set_writeconsole_warnerror(add_to_stderr)
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

    @property
    def stress_value(self): return self.stdout[-1]
