# Futures #
from __future__ import division

# Built-in modules #
import shutil, multiprocessing

# Internal modules #
from plumbing.autopaths import AutoPaths
from plumbing.cache import property_cached
from plumbing.csv_tables import TSVTable

# Third party modules #
import seqenv

###############################################################################
class Seqenv(object):
    """Base class for Seqenv results processing."""

    all_paths = """
    /abundances.tsv
    /output/
    """

    default_N         = 2000
    default_threshold = 3.0
    default_threads   = 1

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
