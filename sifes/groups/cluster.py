# Built-in modules #

# Internal modules #
import sifes
from sifes.report.clusters  import ClusterReport
from sifes.groups.aggregate import Aggregate
from sifes.otu.uparse       import Uparse
from sifes.taxonomy.crest   import Crest
from sifes.taxonomy.rdp     import Rdp

# Composition #
#from sifes.clustering.composition.custom_rank import CompositionPhyla, CompositionOrder, CompositionClass
#from sifes.clustering.composition.custom_rank import CompositionFamily, CompositionGenus
#from sifes.clustering.composition.tips        import CompositionTips
#from sifes.clustering.statistics import StatsOnTaxonomy

# First party modules #
from fasta import FASTA
from plumbing.cache import property_cached

# Third party modules #
from shell_command import shell_output

###############################################################################
class Cluster(Aggregate):
    """Analyzes a group of samples."""

    all_paths = """
    /logs/
    /reads/all_reads.fasta
    /otus/
    /taxonomy/
    /report/report.pdf
    """

    def __init__(self, name, samples, out_dir=None):
        # Directory #
        if out_dir is None: out_dir = sifes.clusters_dir
        # Super #
        super(self.__class__,self).__init__(name, samples, out_dir)
        # Figure out if it's a project #
        self.project = None
        if set(self.samples) == set(self.first.project.samples):
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
        """This is the first function you should call. It will combine all the
        reads of all the samples of this cluster into one big FASTA file."""
        paths = [sample.filter.results.clean for sample in self]
        shell_output('cat %s > %s' % (' '.join(paths), self.reads))
        assert sum([len(s.filter.results.clean) for s in self]) == self.reads.count
        return self.reads

    @property_cached
    def otus(self):
        """Will put cluster the sequences at 97%."""
        return Uparse(self.reads, self.p.otus_dir)

    @property_cached
    def taxonomy(self):
        """Will predict the taxonomy."""
        return Rdp(self.otus.results.centers,   'internal', self.p.taxonomy_dir)
        return Crest(self.otus.results.centers, 'silva123', self.p.taxonomy_dir)
