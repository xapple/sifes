# Futures #
from __future__ import division

# Built-in modules #
from collections import defaultdict

# Internal modules #
from sifes.composition  import taxa_table_graphs

# First party modules #
from plumbing.autopaths import AutoPaths, FilePath, DirectoryPath
from plumbing.cache     import property_cached
from plumbing.common    import prepend_to_file

# Third party modules #
import pandas

###############################################################################
class TaxaTable(object):
    """Takes the otu table and assignments to make taxa tables."""

    # Attributes #
    short_name = 'taxa_table'

    all_paths = """
    /taxa_table_domain.tsv
    /taxa_table_phylum.tsv
    /taxa_table_class.tsv
    /taxa_table_order.tsv
    /taxa_table_family.tsv
    /taxa_table_genus.tsv
    /taxa_table_species.tsv
    /graphs/
    """

    def __repr__(self): return '<%s object on %s>' % (self.__class__.__name__, self.centering.results.centers)
    def __nonzero__(self): return bool(self.p.taxa_table_phylum)

    def __init__(self, otu_table, taxonomy, result_dir):
        # Attributes #
        self.otu_table  = otu_table
        self.taxonomy   = taxonomy
        self.result_dir = result_dir
        # Auto paths #
        self.base_dir = self.result_dir
        self.p = AutoPaths(self.base_dir, self.all_paths)

    def run(self):
        # Message #
        print "Making taxa tables with '%s'" % self.otu_table
        # Do it #
        for i, rank_name in enumerate(self.taxonomy.results.rank_names):
            table = self.taxa_table_at_rank(i)
            path  = self.base_dir + 'taxa_table_' + rank_name.lower() + '.tsv'
            # Make a TSV #
            table.to_csv(path.path, sep='\t', encoding='utf-8')
            # The prepend will fail if the line containing column names is omitted #
            prepend_to_file(path, 'X')

    def taxa_table_at_rank(self, rank):
        # Build a new frame #
        result = defaultdict(lambda: defaultdict(int))
        for sample_name, column in self.otu_table.results.otu_table.iterrows():
            for otu_name, count in column.iteritems():
                assignment = self.taxonomy.results.assignments[otu_name]
                taxa_term = assignment[rank]
                result[taxa_term][sample_name] += count
        # Fill the holes #
        result = pandas.DataFrame(result)
        result = result.fillna(0)
        result = result.astype(int)
        # Sort the table by sum #
        sums = result.sum()
        sums = sums.sort_values(ascending=False)
        result = result.reindex_axis(sums.keys(), axis=1)
        # Return #
        return result

    @property_cached
    def results(self):
        results = TaxaTableResults(self)
        message = "You can't access results from the taxa table before running the algorithm."
        if not results: raise Exception(message)
        return results

###############################################################################
class TaxaTableResults(object):

    def __nonzero__(self): return bool(self.p.taxa_table_phylum)

    def __init__(self, table):
        # Attributes #
        self.table    = table
        self.p        = table.p
        self.taxonomy = table.taxonomy

    def load_table(self, path):
        return pandas.io.parsers.read_csv(path, sep='\t', index_col=0, encoding='utf-8')

    @property_cached
    def taxa_table_domain(self):  return self.load_table(self.p.domain)
    @property_cached
    def taxa_table_phylum(self):  return self.load_table(self.p.phylum)
    @property_cached
    def taxa_table_class(self):   return self.load_table(self.p('class'))
    @property_cached
    def taxa_table_order(self):   return self.load_table(self.p.order)
    @property_cached
    def taxa_table_family(self):  return self.load_table(self.p.family)
    @property_cached
    def taxa_table_genus(self):   return self.load_table(self.p.genus)
    @property_cached
    def taxa_table_species(self): return self.load_table(self.p.species)

    @property_cached
    def taxa_tables_by_rank(self):
        return [getattr(self, 'taxa_table_' + n.lower()) for n in self.taxonomy.results.rank_names]

    @property_cached
    def graphs(self):
        """Sorry for the black magic. The result is an object whose attributes
        are all the graphs found in taxa_table_graphs.py initialized with this
        instance as only argument."""
        class Graphs(object): pass
        result = Graphs()
        for graph in taxa_table_graphs.__all__:
            cls = getattr(taxa_table_graphs, graph)
            setattr(result, cls.short_name, cls(self))
        return result

