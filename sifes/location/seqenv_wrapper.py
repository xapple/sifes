# Futures #
from __future__ import division

# Built-in modules #
import shutil

# Internal modules #
from sifes.location import seqenv_graphs

# First party modules #
import seqenv
from seqsearch.databases.nt import nt
from plumbing.autopaths     import AutoPaths
from plumbing.cache         import property_cached
from plumbing.csv_tables    import TSVTable

# Third party modules #
import pandas

# Constants #
class Dummy(object): pass

###############################################################################
class Seqenv(object):
    """Base class for Seqenv results processing."""

    # Parameters #
    short_name = 'seqenv'
    long_name  = 'Seqenv version ' + seqenv.__version__
    article    = "https://peerj.com/articles/2690/"
    base_url   = "https://github.com/xapple/seqenv"
    database   = nt.long_name

    # Paths #
    all_paths = """
    /abundances.tsv
    /output/samples_to_names.tsv
    /graphs/
    """

    # Options #
    default_N         = 2000
    default_threshold = 3.0
    default_threads   = 1

    def __nonzero__(self): return bool(self.p.samples_to_names)

    def __init__(self, parent, base_dir=None):
        # Parent #
        self.cluster, self.parent = parent, parent
        # Other #
        self.N         = self.default_N
        self.threshold = self.default_threshold
        self.cpus      = self.default_threads
        # Directory #
        if base_dir is None: self.base_dir = self.parent.p.seqenv_dir
        else: self.base_dir = base_dir
        self.p = AutoPaths(self.base_dir, self.all_paths)
        # Files #
        self.abundances = TSVTable(self.p.abundances)

    @property_cached
    def analysis(self):
        return seqenv.Analysis(self.parent.otu_table.results.centers,
                               out_dir      = self.p.output_dir,
                               abundances   = self.abundances,
                               N            = self.N,
                               num_threads  = self.cpus,
                               min_identity = (100 - self.threshold) / 100)

    def run(self, cleanup=True):
        # Clean up #
        if cleanup: shutil.rmtree(self.p.output_dir)
        # Make the abundances file #
        self.parent.otu_table.results.p.norm.transpose(self.abundances, d='\t')
        # Do it #
        self.analysis.run()

    @property_cached
    def results(self):
        results = SeqenvResults(self)
        message = "You can't access results from seqenv before running the algorithm."
        if not results: raise Exception(message)
        return results

###############################################################################
class SeqenvResults(object):

    def __nonzero__(self): return bool(self.p.samples_to_names)

    def __init__(self, wrapper):
        # Attributes #
        self.wrapper  = wrapper
        self.p        = wrapper.p
        self.analysis = wrapper.analysis
        self.cluster  = wrapper.cluster
        self.samples  = wrapper.cluster.samples

    @property_cached
    def samples_to_names(self):
        return pandas.io.parsers.read_csv(self.p.samples_to_names, sep='\t', index_col=0, encoding='utf-8')

    @property_cached
    def graphs(self):
        """Sorry for the black magic. The result is an object whose attributes
        are all the graphs found in otu_table_graphs.py initialized with this
        instance as only argument."""
        result = Dummy()
        for graph in seqenv_graphs.__all__:
            cls = getattr(seqenv_graphs, graph)
            setattr(result, cls.short_name, cls(self))
        return result
