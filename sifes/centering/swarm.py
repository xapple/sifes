# Futures #
from __future__ import division

# Built-in modules #
import multiprocessing

# Internal modules #

# First party modules #
from plumbing.autopaths import AutoPaths, FilePath
from plumbing.cache     import property_cached, LazyString
from fasta              import FASTA, SizesFASTA

# Third party modules #
import sh

# Constants #

###############################################################################
class Swarm(object):
    """Will use SWARM to create OTU from a given FASTA file."""

    # Attributes #
    short_name = 'swarm'
    long_name  = 'SWARM single-linkage clustering method'
    executable = 'swarm'
    url        = 'https://github.com/torognes/swarm'
    article    = 'https://peerj.com/articles/593/'
    version    = '2.1.13'

    all_paths = """
    /centers.fasta
    /details.txt
    /statistics.txt
    /stdout.txt
    /stderr.txt
    /graphs/
    """

    def __repr__(self): return '<%s object on %s>' % (self.__class__.__name__, self.reads)
    def __nonzero__(self): return bool(self.p.readmap)

    def __init__(self, reads, out_dir):
        # Attributes #
        self.reads = reads
        # Paths #
        self.base_dir = out_dir + self.short_name + '/'
        self.p = AutoPaths(self.base_dir, self.all_paths)
        # Output file #
        self.centers = FASTA(self.p.centers)

    def run(self, cpus=1, verbose=True):
        # Message #
        if verbose: print "Making OTUs on '%s' with '%s'" % (self.reads, self.short_name)
        # Number of cores #
        if cpus is None: cpus = min(multiprocessing.cpu_count(), 32)
        # Check version #
        assert "Swarm " + self.version in sh.swarm('-v').stderr
        # Launch #
        sh.swarm('--output-file',     self.p.details,
                 '--seeds',           self.centers,
                 '--threads',         cpus,
                 '--statistics-file', self.p.statistics,
                 '--usearch-abundance',
                 self.reads,
                 _out=self.p.stdout, _err=self.p.stderr)
        # Rename the centers #
        #self.centers.rename_with_num('OTU-')
        # Map the reads back to the centers #
        #pass
        # Checks #
        #assert len(self.reads) == len(self.xxx)
        #assert len(self.reads) == len(self.xxx)

    @property_cached
    def results(self):
        results = SwarmResults(self)
        message = "You can't access results from SWARM before running the algorithm."
        if not results: raise Exception(message)
        return results

###############################################################################
class SwarmResults(object):

    def __nonzero__(self): return bool(self.p)

    def __init__(self, parent):
        self.parent  = parent
        self.p       = parent.p
        self.centers = parent.centers

    @property_cached
    def cluster_counts_table(self):
        pass