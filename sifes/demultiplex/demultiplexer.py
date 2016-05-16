"""
Special module to demultiplex samples.
"""

# Built-in modules #

# Internal modules #

# First party modules #
from plumbing.autopaths import DirectoryPath
from plumbing.cache     import property_cached
from fasta              import PairedFASTQ

# Third party modules #
from shell_command import shell_output

# Constants #

###############################################################################
class Demultiplexer(object):

    def __repr__(self): return '<%s object with %i samples>' % (self.__class__.__name__, len(self.proj))

    def __init__(self, plexed, proj):
        # Attributes #
        self.plexed = plexed
        self.proj   = proj

    def run(self):
        # Merge lanes together #
        for lane in self.lane_pools: lane.merge()
        # Search for barcodes #
        pass

    @property_cached
    def lane_pools(self):
        """We got a few runs that had several lanes in the same sample directory.
        We want to `cat` all these to files called fwd.fastq.gz and rev.fastq.gz"""
        assert all([s.info.get('multiplex_group') for s in self.plexed])
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

    def __repr__(self): return '<%s object with %i samples>' % (self.__class__.__name__, len(self.proj))

    def __init__(self, name, samples):
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
