"""
Special module to demultiplex samples.
"""

# Built-in modules #
import re

# Internal modules #

# First party modules #
from plumbing.autopaths import DirectoryPath
from plumbing.cache     import property_cached
from fasta              import PairedFASTQ
from fasta.primers      import iupac

# Third party modules #
import regex
from shell_command import shell_output
from tqdm import tqdm
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq

# Constants #

###############################################################################
class Demultiplexer(object):

    # Parameters #
    barcode_mismatches = 0

    def __repr__(self): return '<%s object with %i samples>' % (self.__class__.__name__, len(self.proj))

    def __init__(self, plexed, proj):
        # Attributes #
        self.plexed = plexed
        self.proj   = proj
        # Check #
        assert all([s.info.get('multiplex_group') for s in self.plexed])
        assert all([s.info.get('multiplexed_in')  for s in self.proj])

    def run(self):
        # Merge lanes together #
        for lane in self.lane_pools: lane.merge()
        # Search for barcodes #
        for s in tqdm(self.proj):
            lane = [p for p in self.lane_pools if p.name == s.info.get('multiplexed_in')][0]
            s.pair.create()
            barcode = s.barcode
            for f,r in lane.pair:
                start_pos, end_pos = barcode.search(f)
                if True: s.pair.add(f,r)
                else: continue
            s.pair.close()

    @property_cached
    def lane_pools(self):
        """We got a few runs that had several lanes in the same sample directory.
        We want to `cat` all these to files called fwd.fastq.gz and rev.fastq.gz"""
        pools = sorted(list(set([s.info['multiplex_group'] for s in self.plexed])))
        pools = [LanePool(name, [s for s in self.plexed if s.info['multiplex_group']==name]) for name in pools]
        return pools

    @property_cached
    def results(self):
        results = DemultiplexerResults(self)
        message = "You can't access results from demultiplexing before running the algorithm."
        if not results: raise Exception(message)
        return results

###############################################################################
class DemultiplexerResults(object):

    def __nonzero__(self): return bool(False)

    def __init__(self, parent):
        self.parent = parent

###############################################################################
class LanePool(object):
    """Several lanes that need to be merged"""

    def __repr__(self): return '<%s object with %i samples>' % (self.__class__.__name__, len(self.samples))
    def __nonzero__(self): return bool(self.pair)

    def __init__(self, name, samples):
        # Check #
        assert all([s.info.get('multiplex_group') == name for s in samples])
        # Attributes #
        self.name      = name
        self.samples   = samples
        self.project   = samples[0].project
        self.base_dir  = DirectoryPath(self.project.p.lane_cat_dir + self.name + '/')
        self.fwd_path  = self.base_dir + 'fwd.fastq.gz'
        self.rev_path  = self.base_dir + 'rev.fastq.gz'
        self.pair      = PairedFASTQ(self.fwd_path, self.rev_path)
        self.fwd_files = [s.pair.fwd for s in self.samples]
        self.rev_files = [s.pair.rev for s in self.samples]

    def merge_lanes(self):
        # Do it #
        shell_output("zcat %s |gzip > %s" % (' '.join(self.fwd_files), self.pair.fwd))
        shell_output("zcat %s |gzip > %s" % (' '.join(self.rev_files), self.pair.rev))
        # Check #
        assert self.pair.fwd.first.id == self.pair.rev.first.id

###############################################################################
class SingleBarcode(object):
    """Useful when samples are actually combined lanes and you need to find the
    barcodes yourself."""

    def __len__(self): return len(self.string)

    def __init__(self, string, mismatches=0):
        self.string     = string
        self.mismatches = mismatches
        self.seq        = Seq(self.string, IUPAC.ambiguous_dna)
        self.pattern    = ''.join(['[' + iupac[char] + ']' for char in self.fwd_seq])
        self.regex      = regex.compile("(%s){s<=%i}" % (self.pattern, self.mismatches))

    def search(self, read):
        match     = self.regex.search(str(read.seq))
        start_pos = self.match.start() if match else None
        end_pos   = self.match.end()   if match else None
        return start_pos, end_pos