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

###############################################################################
class TaxaTable(object):
    """Takes centers and assignments and makes a taxa table."""

    # Attributes #
    short_name = 'taxa_table'

    all_paths = """
    /taxa_table.csv
    """

    def __repr__(self): return '<%s object on %s>' % (self.__class__.__name__, self.centering.results.centers)
    def __nonzero__(self): return bool(self.p.flat)

    def __init__(self, centering, taxonomy, result_dir):
        # Attributes #
        self.centering  = centering
        self.taxonomy   = taxonomy
        self.result_dir = result_dir
        # Auto paths #
        self.base_dir = self.result_dir + '/'
        self.p = AutoPaths(self.base_dir, self.all_paths)

    def run(self):
        # Message #
        print "Making taxa table with '%s'" % self.centering.results.centers
        # Do it #
        pass

    @property_cached
    def results(self):
        results = TaxaTableResults(self)
        message = "You can't access results from the taxa table before running the algorithm."
        if not results: raise Exception(message)
        return results

###############################################################################
class TaxaTableResults(object):

    def __nonzero__(self): return bool(self.p.flat)

    def __init__(self, table):
        # Attributes #
        self.table = table
        self.p     = table.p
