# Built-in modules #

# Internal modules #
import sifes
from sifes.location.seqenv_wrapper    import Seqenv
from sifes.report.clusters            import ClusterReport
from sifes.groups.aggregate           import Aggregate
from sifes.groups                     import cluster_graphs
from sifes.groups.cluster_graphs      import ClusterLocationMap
from sifes.centering.uparse           import Uparse
from sifes.taxonomy.crest             import Crest
from sifes.taxonomy.rdp               import Rdp
from sifes.taxonomy.mothur_classify   import MothurClassify
from sifes.taxonomy.qiime_classify    import QiimeClassify
from sifes.otus.otu_table             import OtuTable
from sifes.composition.taxa_table     import TaxaTable
from sifes.composition.sub_taxa_table import SubTaxaTable
from sifes.statistics.nmds            import GraphNMDS
from sifes.distances.unifrac          import UniFrac

# First party modules #
from fasta import FASTA
from plumbing.cache import property_cached

# Third party modules #
import numpy
from shell_command import shell_output

# Constants #
class Dummy(object): pass

###############################################################################
class Cluster(Aggregate):
    """Analyzes a group of samples."""

    # Parameters #
    default_taxonomy = "mothur"
    read_count_cutoff_factor = 0.01

    # Options #
    sub_taxa = []

    # Paths #
    all_paths = """
    /logs/
    /reads/all_reads.fasta
    /centers/
    /otu_table/
    /taxa_table/
    /taxonomy/
    /sub_taxa/
    /distances/
    /seqenv/
    /graphs/
    /report/report.pdf
    /report/replicates.pdf
    """

    def __nonzero__(self):
        """When we havn't run anything yet, return False."""
        return bool(self.reads)

    def __init__(self, name, samples, out_dir=None):
        # Directory #
        if out_dir is None: out_dir = sifes.clusters_dir
        # Compute cutoff for throwing away samples #
        read_counts            = [s.clean.count for s in samples]
        self.read_count_cutoff = numpy.mean(read_counts) * self.read_count_cutoff_factor
        # Auto-filter low read count samples #
        self.good_samples = [s for s in samples if len(s.clean) >= self.read_count_cutoff]
        self.bad_samples  = [s for s in samples if len(s.clean) < self.read_count_cutoff]
        self.num_dropped_samples = len(self.bad_samples)
        # Super #
        super(self.__class__,self).__init__(name, self.good_samples, out_dir)
        # Figure out if it's a subset of a project #
        self.project = None
        if set(self.samples) <= set(self.first.project.samples):
            self.project = self.first.project
        # FASTA #
        self.reads = FASTA(self.p.all_reads)

    def combine_reads(self):
        """This is the first method you should call. It will combine all the
        reads of all the samples of this cluster into one big FASTA file."""
        print "Combining all reads for cluster '%s' (%i samples)" % (self.name, len(self.samples))
        paths = [sample.filter.results.clean for sample in self]
        shell_output('cat %s > %s' % (' '.join(paths), self.reads))
        assert sum([len(s.filter.results.clean) for s in self]) == self.reads.count
        return self.reads

    @property_cached
    def centering(self):
        """Will cluster the sequences at 97% and pick centers."""
        return Uparse(self.reads, self.p.centers_dir)

    @property_cached
    def taxonomy(self):
        """Will predict the taxonomy."""
        choices = {'mothur': (MothurClassify, (self.centering.results.centers, self.p.taxonomy_dir)),
                   'qiime':  (QiimeClassify,  (self.centering.results.centers, self.p.taxonomy_dir)),
                   'rdp':    (Rdp,            (self.centering.results.centers, self.p.taxonomy_dir)),
                   'crest':  (Crest,          (self.centering.results.centers, self.p.taxonomy_dir))}
        cls, params = choices.get(self.default_taxonomy)
        return cls(*params)

    @property_cached
    def otu_table(self):
        """Will produce the OTU table."""
        return OtuTable(self.centering, self.taxonomy, self.p.otu_table_dir)

    @property_cached
    def taxa_table(self):
        """Will produce the taxonomy-based tables."""
        return TaxaTable(self.otu_table, self.taxonomy, self, self.p.taxa_table_dir)

    @property_cached
    def sub_taxa_tables(self):
        """Will produce the sub taxa tables."""
        return [SubTaxaTable(self.taxa_table, rank, taxa, self, self.p.sub_taxa_dir) for rank,taxa in self.sub_taxa]

    @property_cached
    def down_to(self):
        """The number that we need to rarefy each sample so that diversity is comparable."""
        return min(sum(s.otu_counts) for s in self.samples)

    @property_cached
    def seqenv(self):
        """Will produce the isolation source linear combination predictions."""
        return Seqenv(self, self.p.seqenv_dir)

    @property_cached
    def unifrac_matrix(self):
        """Will produce the UniFrac distance matrix."""
        return UniFrac(self.otu_table, self.p.distances_dir)

    @property_cached
    def nmds_graph(self):
        """Non-metric multidimensional scaling. Using the information in the OTU table and a
         distance metric such as the Horn 1966 (adapted from Morisita 1959)" one,
         will place every sample on an ordination plot."""
        return GraphNMDS(self, self.p.graphs_dir)

    @property_cached
    def redundancy(self):
        """Will produce the Redundancy Analysis (RDA)."""
        return RedundancyAnalysis(self)

    @property_cached
    def locations_maps(self):
        """Make as many maps as there are custom groupings."""
        maps   = []
        groups = set(s.grouping for s in self.samples)
        for g in groups:
            samples = [s for s in self.samples if s.grouping == g]
            short_name = 'location_map_' + g.lower()
            maps.append(ClusterLocationMap(g, samples, self, short_name=short_name))
        return maps

    @property_cached
    def graphs(self):
        """Sorry for the black magic. The result is an object whose attributes
        are all the graphs found in cluster_graphs.py initialized with this
        instance as only argument."""
        result = Dummy()
        for graph in cluster_graphs.__all__:
            cls = getattr(cluster_graphs, graph)
            setattr(result, cls.short_name, cls(self))
        return result

    @property_cached
    def report(self):
        """The PDF report."""
        return ClusterReport(self)