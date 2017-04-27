# Futures #
from __future__ import division

# Built-in modules #
from collections import defaultdict

# Internal modules #
from sifes.composition.taxa_table_graphs import TaxaBarstack

# First party modules #
from plumbing.autopaths import AutoPaths, FilePath, DirectoryPath
from plumbing.cache     import property_cached
from plumbing.common    import prepend_to_file

# Third party modules #
import pandas

# Graphs #
class Dummy(object): pass

###############################################################################
class TaxaTable(object):
    """Takes the otu table and assignments to make taxa tables."""

    # Attributes #
    short_name = 'taxa_table'

    all_paths = """
    /taxa_table_domain.tsv
    /taxa_table_kingdom.tsv
    /taxa_table_phylum.tsv
    /taxa_table_class.tsv
    /taxa_table_order.tsv
    /taxa_table_family.tsv
    /taxa_table_tribe.tsv
    /taxa_table_genus.tsv
    /taxa_table_species.tsv
    /graphs/
    """

    def __repr__(self): return '<%s object on %s>' % (self.__class__.__name__, self.result_dir)
    def __nonzero__(self): return bool(self.p.taxa_table_phylum)

    def __init__(self, otu_table, taxonomy, result_dir):
        # Attributes #
        self.otu_table  = otu_table
        self.taxonomy   = taxonomy
        self.result_dir = result_dir
        # Short cuts #
        self.rank_names  = self.taxonomy.results.rank_names
        self.otu_df      = self.otu_table.results.otu_table
        self.assignments = self.taxonomy.results.assignments
        # Auto paths #
        self.base_dir = self.result_dir
        self.p        = AutoPaths(self.base_dir, self.all_paths)

    def run(self, verbose=False):
        # Message #
        if verbose: print "Making all taxa tables in '%s'" % self.base_dir
        # Make directory #
        self.base_dir.create_if_not_exists()
        # Do it #
        for i, rank_name in enumerate(self.rank_names):
            table = self.taxa_table_at_rank(i)
            path  = self.base_dir + 'taxa_table_' + rank_name.lower() + '.tsv'
            # Make a TSV #
            table.to_csv(path.path, sep='\t', encoding='utf-8')
            # The prepend will fail if the line containing column names is omitted #
            prepend_to_file(path, 'X')

    def taxa_table_at_rank(self, rank):
        # Build a new frame #
        result = defaultdict(lambda: defaultdict(int))
        for sample_name, column in self.otu_df.iterrows():
            for otu_name, count in column.iteritems():
                assignment = self.assignments[otu_name]
                if rank >= len(assignment): taxa_term = "Unassigned"
                else:                       taxa_term = assignment[rank]
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
        self.table      = table
        self.p          = table.p
        self.base_dir   = table.base_dir
        self.rank_names = table.rank_names

    def load_table(self, path):
        return pandas.io.parsers.read_csv(path, sep='\t', index_col=0, encoding='utf-8')

    @property_cached
    def taxa_tables_by_rank(self):
        return [self.load_table(self.base_dir + 'taxa_table_' + n.lower() + '.tsv') for n in self.rank_names]

    @property_cached
    def graphs(self):
        """The result is an object whose attributes are all the graphs initialized with
        this instance as only argument. The graphs are also in a list."""
        result = Dummy()
        result.by_rank = []
        for i, rank_name in enumerate(self.rank_names):
            attributes = dict(base_rank  = i,
                              label      = rank_name,
                              short_name ='taxa_barstack_' + rank_name.lower())
            graph = type("Composition" + rank_name, (TaxaBarstack,), attributes)(self)
            setattr(result, graph.short_name, graph)
            result.by_rank.append(graph)
        return result
