# Built-in modules #

# Internal modules #

# First party modules #
from plumbing.autopaths import AutoPaths
from plumbing.cache     import property_cached

###############################################################################
class SeqFilter(object):
    """Lorem."""

    short_name = 'seq_filter'
    primer_mismatches = "TODO" #TODO

    all_paths = """
    /x.fasta
    /x.fastq
    """

    def __repr__(self): return '<%s object on %s>' % (self.__class__.__name__, self.input)
    def __nonzero__(self): return bool(self.p.x)

    def __init__(self, input, result_dir):
        # Save attributes #
        self.input      = pair
        self.result_dir = result_dir
        # Auto paths #
        self.base_dir = self.result_dir + self.short_name + '/'
        self.p = AutoPaths(self.base_dir, self.all_paths)

    def run(self):
        pass

    @property_cached
    def results(self):
        results = PandaseqResults(self)
        message = "You can't access results from filtering before running the algorithm."
        if not results: raise Exception(message)
        return results

###############################################################################
class SeqFilterResults(object):

    def __nonzero__(self): return bool(self.x)

    def __init__(self, x):
        self.x    = x
