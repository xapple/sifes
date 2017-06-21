# Futures #
from __future__ import division

# Built-in modules #
import multiprocessing, re, base64, hashlib
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
    /counts.tsv
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
        self.clusters = FASTA(self.p.clusters)
        self.centers  = FASTA(self.p.centers)

    def run(self, cpus=None, verbose=True):
        # Message #
        if verbose: print "Making OTUs on '%s' with '%s'" % (self.reads, self.short_name)
        # Number of cores #
        if cpus is None: cpus = min(multiprocessing.cpu_count(), 32)
        # Check version #
        assert "Swarm " + self.version in sh.swarm('-v').stderr
        # Dereplicate sequences, but remember everything along the way #
        unique_seqs = defaultdict(list)
        for read in self.reads: unique_seqs[str(read.seq)] += [read.id]
        # Write the result to a file #
        self.derep.create()
        hash = lambda x: base64.b32encode(hashlib.sha1(x).digest())[0:10]
        for seq, samples in unique_seqs.items():
            if len(samples) == 1: continue
            self.derep.add_str(seq, hash(seq), len(samples))
        self.derep.close()
        # Launch #
        sh.swarm('--output-file',     self.p.details,
                 '--seeds',           self.clusters,
                 '--threads',         cpus,
                 '--statistics-file', self.p.statistics,
                 '--fastidious',      # link nearby low-abundance swarms
                 '--usearch-abundance',
                 self.derep,
                 _out=self.p.stdout, _err=self.p.stderr)
        # Rename the centers #
        self.clusters.rename_with_num('OTU-', self.centers)
        # Use the information in 'details.txt' to get the final counts #
        result = defaultdict(lambda: defaultdict(int))
        # Loop #
        for i, line in enumerate(self.p.details):
            target = 'OTU-%i' % i
            for item in line.split():
                hash, size = re.findall("\A(\w+);size=(\d+);\Z", item)[0]
                read_ids   = unique_seqs[hash]
                print read_ids #TODO empty
                for read_id in read_ids:
                    result[target][read_id.split(':')[0]] += 1
        # Return #
        result = pandas.DataFrame(result)
        result = result.fillna(0)
        result = result.astype(int)
        result = result.reindex_axis(sorted(result.columns, key=natural_sort), axis=1)
        result.to_csv(self.p.counts.path, sep='\t', encoding='utf-8')

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
        return pandas.io.parsers.read_csv(self.p.counts.path, sep='\t', index_col=0, encoding='utf-8')