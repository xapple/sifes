# Built-in modules #
from collections import defaultdict

# Internal modules #
import illumitag
from plumbing.autopaths import AutoPaths
from plumbing.csv_tables import CSVTable
from plumbing.cache import property_cached
from illumitag.clustering.statistics import StatsOnComposition

# Third party modules #
import pandas

###############################################################################
class Composition(object):
    """Base class for taxonomic compositing"""

    all_paths = """
    /taxa_table.csv
    /graphs/
    /stats/
    """

    def __len__(self): return len(set(self.taxonomy.assignments[o] for o in self.taxonomy.otu_table))

    def __init__(self, parent, base_dir=None):
        # Parent #
        self.taxonomy, self.parent = parent, parent
        # Inherited #
        self.samples = self.parent.samples
        # Dir #
        if base_dir is None: self.base_dir = self.parent.p.composition_dir
        else: self.base_dir = base_dir
        self.p = AutoPaths(self.base_dir, self.all_paths)
        # Taxa table #
        self.taxa_csv = CSVTable(self.p.taxa_csv)
        # Stats #
        self.stats = StatsOnComposition(self)

    def make_plots(self, **kwargs):
        for graph in self.graphs: graph.plot(**kwargs)

    def make_taxa_table(self):
        """Convert to CSV"""
        self.taxa_table.to_csv(str(self.taxa_csv), sep='\t', encoding='utf-8')

###############################################################################
class SimpleComposition(Composition):
    def __init__(self, parent, base_dir):
        self.parent, self.taxonomy = parent, parent
        self.base_dir = base_dir
        self.p = AutoPaths(self.base_dir, self.all_paths)
        # Simple graph #
        self.graph = illumitag.clustering.composition.plots.TaxaBarstack(self)
        self.graph.bottom = 0.40
        self.graph.legend_anchor = -0.3
        self.formats = ('pdf',)

    @property_cached
    def taxa_table(self):
        result = defaultdict(int)
        for taxa in self.parent.assignments.values():
            phyla = taxa[2] if len(taxa) > 2 else taxa[-1]
            result[phyla] += 1
        return pandas.DataFrame(result, index=["simple"])