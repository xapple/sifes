# Futures #
from __future__ import division

# Built-in modules #

# Internal modules #
from sifes.otus import otu_table_graphs

# First party modules #
from fasta import FASTA
from plumbing.autopaths import AutoPaths
from plumbing.common    import prepend_to_file
from plumbing.cache     import property_cached

# Third party modules #
import pandas, numpy

# Constants #
class Dummy(object): pass

###############################################################################
class OtuTable(object):
    """Takes the centers and the assignments and makes an OTU table."""

    # Attributes #
    short_name = 'otu_table'

    # Properties #
    unwanted_taxa = ['Plastid', 'Mitochondrion'] #, 'Thaumarchaeota', 'Crenarchaeota', 'Euryarchaeota']

    all_paths = """
    /otu_table_flat.tsv
    /otu_table_norm.tsv
    /centers_filtered.fasta
    /graphs/
    """

    def __repr__(self): return '<%s object on %s>' % (self.__class__.__name__, self.centering.results.centers)
    def __nonzero__(self): return bool(self.p.flat)

    def __init__(self, centering, taxonomy, result_dir):
        # Attributes #
        self.centering  = centering
        self.taxonomy   = taxonomy
        self.result_dir = result_dir
        # Auto paths #
        self.base_dir = self.result_dir
        self.p = AutoPaths(self.base_dir, self.all_paths)

    def run(self):
        # Message #
        print "Making OTU table with '%s'" % self.centering.results.centers
        # Do it #
        self.make_otu_table()
        self.make_otu_table_norm()
        self.make_filtered_centers()

    def make_otu_table(self):
        """
        Ask the counts from the OTU making class and do some modifications.
        OTUs are columns and sample names are rows.
        """
        # Main objects #
        cluster_table = self.centering.results.cluster_counts_table.copy()
        assignments   = self.taxonomy.results.assignments
        # Remove unwanted #
        for otu_name in cluster_table:
            species = assignments[otu_name]
            if any(bad in species for bad in self.unwanted_taxa):
                cluster_table = cluster_table.drop(otu_name, 1)
        # Convert to CSV #
        cluster_table.to_csv(self.p.flat.path, sep='\t', encoding='utf-8')
        # Symbolic X in the upper left corner #
        prepend_to_file(self.p.flat, 'X')

    def make_otu_table_norm(self):
        """Convert to CSV. OTUs are columns and sample names are rows."""
        normed = self.results.otu_table.apply(lambda x: x/x.sum(), axis=1).replace(numpy.inf, 0.0)
        normed.to_csv(self.p.norm.path, sep='\t', float_format='%.5g', encoding='utf-8')
        # Symbolic X in the corner #
        prepend_to_file(self.p.norm, 'X')

    def make_filtered_centers(self):
        """
        Regenerate the centers file with only the OTUs that haven't been
        filtered out previously.
        """
        self.otus_to_keep = [otu for otu in self.results.otu_table]
        def filter_otus(seqs):
            for seq in seqs:
                if seq.id in self.otus_to_keep: yield seq
        self.results.centers.write(filter_otus(self.centering.results.centers))

    @property_cached
    def results(self):
        results = OtuTableResults(self)
        message = "You can't access results from the otu table before running the algorithm."
        if not results: raise Exception(message)
        return results

###############################################################################
class OtuTableResults(object):

    def __nonzero__(self): return bool(self.p.flat)

    def __init__(self, table):
        # Attributes #
        self.table = table
        self.p     = table.p
        # Auto paths #
        self.centers = FASTA(self.p.centers)

    @property_cached
    def otu_table(self):
        """OTUs are columns and sample names are rows."""
        return pandas.io.parsers.read_csv(self.p.flat, sep='\t', index_col=0, encoding='utf-8')

    @property
    def otu_table_norm(self):
        """The same thing as otu_table but normalized so that the sum of a sample is always one"""
        return pandas.io.parsers.read_csv(self.p.norm, sep='\t', index_col=0, encoding='utf-8')

    @property_cached
    def graphs(self):
        """Sorry for the black magic. The result is an object whose attributes
        are all the graphs found in otu_table_graphs.py initialized with this
        instance as only argument."""
        result = Dummy()
        for graph in otu_table_graphs.__all__:
            cls = getattr(otu_table_graphs, graph)
            setattr(result, cls.short_name, cls(self))
        return result
