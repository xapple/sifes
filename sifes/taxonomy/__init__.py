# Futures #
from __future__ import division

# Built-in modules #

# Internal modules #
from plumbing.common import prepend_to_file
from plumbing.cache import property_cached
from fasta import FASTA

# Third party modules #
import pandas, numpy

###############################################################################
class Taxonomy(object):
    """Can assign taxonomy to a FASTA file of 16S sequences."""

    # Parameters #
    unwanted = ['Plastid', 'Mitochondrion', 'Thaumarchaeota', 'Crenarchaeota', 'Euryarchaeota']

    def __repr__(self): return '<%s object of %s>' % (self.__class__.__name__, self.parent)

    def at_level(self, level):
        return dict((k,v[level]) for k,v in self.assignments.items() if len(v) > level)

    def make_plots(self):
        for graph in self.graphs: graph.plot()

    def make_otu_table(self):
        """Ask the counts from the OTU class and do some modifications.
        OTUs are columns and sample names are rows."""
        # Remove unwanted #
        result = self.otu.cluster_counts_table.copy()
        for otu_name in result:
            species = self.assignments.get(otu_name, None)
            if len(species) > 2 and species[2] in self.unwanted: result = result.drop(otu_name, 1)
        # Convert to CSV #
        result.to_csv(str(self.otu_csv), sep='\t', encoding='utf-8')
        prepend_to_file(self.otu_csv, 'X')

    @property_cached
    def otu_table(self):
        """OTUs are columns and sample names are rows."""
        return pandas.io.parsers.read_csv(self.otu_csv, sep='\t', index_col=0, encoding='utf-8')

    @property
    def otu_table_norm(self):
        """The same thing as otu_table but normalized so that the sum of a sample is always one"""
        return self.otu_table.apply(lambda x: x/x.sum(), axis=1).replace(numpy.inf, 0.0)

    def make_otu_table_norm(self):
        """Convert to CSV. OTUs are columns and sample names are rows."""
        self.otu_table_norm.to_csv(self.otu_csv_norm.path, sep='\t', float_format='%.5g', encoding='utf-8')
        prepend_to_file(self.otu_csv_norm, 'X')

    def make_filtered_centers(self):
        """Regenerate the centers file with only the OTUs that haven't been
        filtered out previously."""
        self.otus_to_keep = [otu for otu in self.otu_table]
        def filter_otus(f):
            for seq in f:
                if seq.id in self.otus_to_keep: yield seq
        self.centers.write(filter_otus(self.otu.centers))