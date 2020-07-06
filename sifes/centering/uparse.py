# Futures #
from __future__ import division

# Built-in modules #
import multiprocessing

# Internal modules #
from sifes.centering import UClusterFile

# First party modules #
from plumbing.common      import natural_sort
from autopaths.auto_paths import AutoPaths
from plumbing.cache       import property_cached, LazyString
from fasta                import FASTA, SizesFASTA

# Third party modules #
import sh, pandas

# Constants #
uparse_version = LazyString(lambda: sh.usearch8('-version').stdout[8:].strip('\n'))

###############################################################################
class Uparse(object):
    """Will use UPARSE to create OTU clusters from a given FASTA file."""

    # Attributes #
    short_name = 'uparse'
    long_name  = 'UPARSE denovo picking'
    executable = 'usearch8'
    article    = "http://www.nature.com/doifinder/10.1038/nmeth.2604"
    version    = uparse_version

    # Parameters #
    threshold = 3.0

    all_paths = """
    /derep.fasta
    /sorted.fasta
    /centers.fasta
    /readmap.uc
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
        self.derep   = SizesFASTA(self.p.derep)
        self.sorted  = SizesFASTA(self.p.sorted)
        self.centers = FASTA(self.p.centers)
        self.readmap = UClusterFile(self.p.readmap)

    def run(self, threshold=None, cpus=1):
        # Message #
        print("Making OTUs on '%s' with '%s'" % (self.reads, self.short_name))
        # Optional threshold #
        if threshold is None: threshold = self.threshold
        # Number of cores #
        if cpus is None: cpus = min(multiprocessing.cpu_count(), 32)
        # Identity cutoff in percent #
        identity = (100 - threshold) / 100
        # Dereplicate (uparse version 32-bit version can run out of memory) #
        sh.usearch8("--derep_fulllength", self.reads,
                    '-fastaout',          self.derep,
                    '-sizeout',
                    '-threads', cpus)
        # Order by size and kill singletons (likely chimeras) #
        sh.usearch8("--sortbysize", self.derep,
                    '-fastaout',    self.sorted,
                    '-minsize',     2,
                    '-threads',     cpus)
        # Compute the centers #
        sh.usearch8("--cluster_otus",  self.sorted,
                    '-otus',           self.centers,
                    '-otu_radius_pct', threshold,
                    '-threads',        cpus)
        # Rename the centers #
        self.centers.rename_with_num('OTU-')
        # Check that we don't have a file size problem at this point #
        assert self.reads.count_bytes < 2**32
        # Map the reads back to the centers #
        sh.usearch8("-usearch_global", self.reads,
                    '-db',             self.centers,
                    '-strand',         'plus',
                    '-id',             identity,
                    '-uc',             self.readmap,
                    '-threads',        cpus)
        # Checks #
        assert len(self.reads) == len(self.derep)
        assert len(self.reads) == len(self.readmap)

    @property_cached
    def results(self):
        results = UparseResults(self)
        message = "You can't access results from UPARSE before running the algorithm."
        if not results: raise Exception(message)
        return results

###############################################################################
class UparseResults(object):

    def __nonzero__(self): return bool(self.p.readmap)

    def __init__(self, parent):
        self.parent  = parent
        self.p       = parent.p
        self.centers = parent.centers

    @property_cached
    def cluster_counts_table(self):
        """Parse that custom output for creating the unfiltered OTU table."""
        result = pandas.DataFrame(self.parent.readmap.otu_sample_counts)
        result = result.fillna(0)
        result = result.astype(int)
        result = result.reindex_axis(sorted(result.columns, key=natural_sort), axis=1)
        return result
