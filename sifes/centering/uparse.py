# Futures #
from __future__ import division

# Built-in modules #
import os, multiprocessing
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
        print "Making OTUs on '%s' with '%s'" % (self.reads, self.short_name)
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

###############################################################################
class UClusterFile(FilePath):
    """A special format outputted by UPARSE
    An example line:
    H       1474    422     97.6    +       0       0       422M    run4_pool1_sample1_read2        OTU-1474

    The corresponding legend:
    # 1=Type, 2=ClusterNr, 3=SeqLength or ClusterSize, 4=PctId, 5=Strand, 6=QueryStart, 7=SeedStart, 8=Alignment, 9=QueryLabel, 10=TargetLabel
    """

    def __len__(self):     return self.count
    def __repr__(self):    return '<%s object on "%s">' % (self.__class__.__name__, self.path)
    def __nonzero__(self): return os.path.getsize(self.path) != 0

    @property
    def count(self):
        return int(sh.wc('-l', self.path).stdout.split()[0])

    @property
    def count_mapped(self):
        """Count only the reads that mapped to a center"""
        return sum([1 for line in self if line.split()[-1] != '*'])

    @property
    def otu_sample_counts(self):
        # Put results in a dict of dicts #
        result = defaultdict(lambda: defaultdict(int))
        # Loop #
        for line in self:
            # Parse the line (the query is the name of the read) #
            kind, num, length, similarity, strand, start, seed, alignment, query, target = line.split()
            # Skip no hits (the target is the OTU name) #
            if target == '*': continue
            # Remove read number to get sample name #
            sample_name = query.split(':')[0]
            # Count #
            result[target][sample_name] += 1
        # Return #
        return result
