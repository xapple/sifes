"""
Special module to demultiplex samples.
"""

# Built-in modules #
import multiprocessing

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

    def __init__(self, plexed, samples):
        # Attributes #
        self.plexed  = plexed
        self.samples = samples
        # Check #
        assert all([s.info.get('multiplex_group') for s in self.plexed])
        assert all([s.info.get('multiplexed_in')  for s in self.samples])
        # Lane pool #
        get_samples = lambda name: [s for s in self.plexed if s.info['multiplex_group'] == name]
        self.lane_pools = sorted(list(set([s.info['multiplex_group'] for s in self.plexed])))
        self.lane_pools = [LanePool(name, get_samples(name)) for name in self.lane_pools]

    def run(self, cpus=None):
        # Number of cores #
        if cpus is None: cpus = min(multiprocessing.cpu_count(), 32)
        # Merge lanes together #
        #thread_pool = multiprocessing.Pool(processes=cpus)
        #thread_pool.map(lambda lane: lane.merge(), self.lane_pools)
        # Search for barcodes #
        for i,s in enumerate(self.samples):
            print "Demultiplexing sample %i out of %i" % (i+1, len(self.samples))
            lane = [p for p in self.lane_pools if p.name == s.info.get('multiplexed_in')][0]
            s.pair.create()
            barcode = SingleBarcode(s.info['custom_barcode'])
            for f,r in tqdm(lane.pair.parse_primers(s.primers, 2)):
                if f.fwd_start_pos:
                    start_pos, end_pos = barcode.search(f.read)
                    if end_pos == f.fwd_start_pos:
                        s.pair.add(f.read, r.read)
            s.pair.close()

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

    def __init__(self, name, input):
        # Check #
        assert all([s.info.get('multiplex_group') == name for s in input])
        # Attributes #
        self.name      = name
        self.input     = input
        self.project   = input[0].project
        self.base_dir  = DirectoryPath(self.project.p.lane_cat_dir + self.name + '/')
        self.fwd_path  = self.base_dir + 'fwd.fastq.gz'
        self.rev_path  = self.base_dir + 'rev.fastq.gz'
        self.pair      = PairedFASTQ(self.fwd_path, self.rev_path)
        self.fwd_files = [s.pair.fwd for s in self.input]
        self.rev_files = [s.pair.rev for s in self.input]
        # Link the input samples back #
        for s in input: s.lane_pool = self

    def merge(self):
        """We got a few runs that had several lanes in the same sample directory.
        We want to `cat` all these to files called fwd.fastq.gz and rev.fastq.gz"""
        # Make output directory #
        self.pair.fwd.directory.create()
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
        self.pattern    = ''.join(['[' + iupac[char] + ']' for char in self.string])
        self.regex      = regex.compile("(%s){s<=%i}" % (self.pattern, self.mismatches))

    def search(self, read):
        match     = self.regex.search(str(read.seq))
        start_pos = match.start() if match else None
        end_pos   = match.end()   if match else None
        return start_pos, end_pos