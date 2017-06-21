# Futures #
from __future__ import division

# Built-in modules #
import multiprocessing, re
from collections import defaultdict

# Internal modules #

# First party modules #
from plumbing.common    import natural_sort
from plumbing.autopaths import AutoPaths, FilePath
from plumbing.cache     import property_cached, LazyString
from fasta              import FASTA, SizesFASTA

# Third party modules #
import sh, pandas

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
    /derep.fasta
    /sorted.fasta
    /centers.fasta
    /clusters.fasta
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
        # Files #
        self.derep    = SizesFASTA(self.p.derep)
        self.sorted   = SizesFASTA(self.p.sorted)
        self.clusters = SizesFASTA(self.p.clusters)
        self.centers  = FASTA(self.p.centers)

    def run(self, cpus=None, verbose=True):
        # Message #
        if verbose: print "Making OTUs on '%s' with '%s'" % (self.reads, self.short_name)
        # Number of cores #
        if cpus is None: cpus = min(multiprocessing.cpu_count(), 32)
        # Check version #
        assert "Swarm " + self.version in sh.swarm('-v').stderr
        # Dereplicate #
        sh.usearch8("--derep_fulllength", self.reads,
                    '-fastaout',          self.derep,
                    '-sizeout',
                    '-threads', cpus)
        # Order by size and kill singletons (likely chimeras) #
        sh.usearch8("--sortbysize", self.derep,
                    '-fastaout',    self.sorted,
                    '-minsize',     2,
                    '-threads',     cpus)
        # Launch #
        sh.swarm('--output-file',     self.p.details,
                 '--seeds',           self.clusters,
                 '--threads',         cpus,
                 '--statistics-file', self.p.statistics,
                 '--fastidious',      # link nearby low-abundance swarms
                 '--usearch-abundance',
                 self.sorted,
                 _out=self.p.stdout, _err=self.p.stderr)
        # Rename the centers #
        self.clusters.rename_with_num('OTU-', self.centers)

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
        """We just need to read the file 'details.txt' to build this table."""
        # Put results in a dict of dicts #
        result = defaultdict(lambda: defaultdict(int))
        # Loop #
        for i, line in enumerate(self.p.details):
            target = 'OTU-%i' % i
            for item in line.split():
                sample_name, num, size = re.findall("\A(\w+):(\d+);size=(\d+);\Z", item)
                result[target][sample_name] += int(size)
        # Return #
        result = pandas.DataFrame(result)
        result = result.fillna(0)
        result = result.astype(int)
        result = result.reindex_axis(sorted(result.columns, key=natural_sort), axis=1)
        return result
