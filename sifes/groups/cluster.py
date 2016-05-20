# Built-in modules #

# Internal modules #
import sifes
from sifes.report.clusters          import ClusterReport
from sifes.groups.aggregate         import Aggregate
from sifes.centering.uparse         import Uparse
from sifes.taxonomy.crest           import Crest
from sifes.taxonomy.rdp             import Rdp
from sifes.taxonomy.mothur_classify import MothurClassify
from sifes.otus.otu_table           import OtuTable
from sifes.composition.taxa_table   import TaxaTable
from sifes.statistics.nmds          import GraphNMDS

# Composition #
#from sifes.clustering.composition.custom_rank import CompositionPhyla, CompositionOrder, CompositionClass
#from sifes.clustering.composition.custom_rank import CompositionFamily, CompositionGenus
#from sifes.clustering.composition.tips        import CompositionTips
#from sifes.clustering.statistics import StatsOnTaxonomy

# First party modules #
from fasta import FASTA
from plumbing.cache import property_cached

# Third party modules #
import numpy
from shell_command import shell_output

###############################################################################
class Cluster(Aggregate):
    """Analyzes a group of samples."""

    # Parameters #
    read_count_cutoff_percentile = 1

    all_paths = """
    /logs/
    /reads/all_reads.fasta
    /centers/
    /otu_table/
    /taxa_table/
    /taxonomy/
    /graphs/
    /report/report.pdf
    """

    def __init__(self, name, samples, out_dir=None):
        # Directory #
        if out_dir is None: out_dir = sifes.clusters_dir
        # Compute cutoff for throwing away samples #
        read_counts            = [s.clean.count for s in samples]
        self.read_count_cutoff = numpy.percentile(read_counts, self.read_count_cutoff_percentile)
        # Auto-filter low read count samples #
        self.good_samples = [s for s in samples if len(s.clean) >= self.read_count_cutoff]
        self.bad_samples  = [s for s in samples if len(s.clean) < self.read_count_cutoff]
        self.num_dropped_samples = len(self.bad_samples)
        # Super #
        super(self.__class__,self).__init__(name, self.good_samples, out_dir)
        # Figure out if it's a subset of a project #
        self.project = None
        if set(self.samples) < set(self.first.project.samples):
            self.project = self.first.project
        # FASTA #
        self.reads = FASTA(self.p.all_reads)
        # Composition tables #
        #self.comp_phyla  = CompositionPhyla(self,  self.p.comp_phyla)
        #self.comp_order  = CompositionOrder(self,  self.p.comp_order)
        #self.comp_class  = CompositionClass(self,  self.p.comp_class)
        #self.comp_family = CompositionFamily(self, self.p.comp_order)
        #self.comp_genus  = CompositionGenus(self,  self.p.comp_class)
        # The complete one #
        #self.comp_tips  = CompositionTips(self,  self.p.comp_tips)
        # Stats #
        #self.stats = StatsOnTaxonomy(self)
        # Source tracking #
        #self.seqenv = Seqenv(self)
        # Report #
        self.report = ClusterReport(self)

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
        return MothurClassify(self.centering.results.centers, 'silva123', self.p.taxonomy_dir)
        return Rdp(           self.centering.results.centers, 'internal', self.p.taxonomy_dir)
        return Crest(         self.centering.results.centers, 'silva123', self.p.taxonomy_dir)

    @property_cached
    def otu_table(self):
        """Will produce the OTU table."""
        return OtuTable(self.centering, self.taxonomy, self.p.otu_table_dir)

    @property_cached
    def taxa_table(self):
        """Will produce the taxonomy-based tables."""
        return TaxaTable(self.otu_table, self.taxonomy, self.p.taxa_table_dir)

    @property_cached
    def nmds_graph(self):
        """Non-metric multidimensional scaling. Using the information in the OTU table and a
         distance metric such as the Horn 1966 (adapted from Morisita 1959)" one,
         will place every sample on an ordination plot."""
        return GraphNMDS(self, self.p.graphs_dir)
