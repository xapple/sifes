# Built-in modules #

# Internal modules #
import sifes
from sifes.report.clusters       import ClusterReport
from sifes.groups.aggregate      import Aggregate
from sifes.clustering.otu.uparse import UparseOTUs

# First party modules #
from plumbing.autopaths import AutoPaths
from fasta import FASTA

# Third party modules #
from shell_command import shell_output

###############################################################################
class Cluster(Aggregate):
    """Analyzes a group of samples."""

    all_paths = """
    /logs/
    /reads/all_reads.fasta
    /otus/
    /report/report.pdf
    """

    def __init__(self, name, samples, base_dir=None):
        # Directory #
        if base_dir is None: base_dir = sifes.clusters_dir
        # Super #
        super(self.__class__,self).__init__(name, samples, base_dir)
        # Figure out if it's a project #
        self.project = None
        if set(self.samples) == set(self.first.project.samples):
            self.project = self.first.pool.project
        # FASTA #
        self.reads = FASTA(self.p.all_reads)
        # OTU picking #
        self.otus_uparse = UparseOTUs(self)
        self.otus        = self.otu_uparse
        # Composing #
        self.report = ClusterReport(self)

    @property
    def first(self): return self.children[0]

    @property
    def count_seq(self):
        return sum([len(sample) for sample in self])

    def combine_reads(self):
        """This is the first function you should call. It will combine all the
        reads of all the samples of this cluster into one big FASTA file."""
        paths = [sample.filter.results.clean for sample in self]
        shell_output('cat %s > %s' % (' '.join(paths), self.reads))
        return self.reads