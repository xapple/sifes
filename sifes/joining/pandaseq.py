# Built-in modules #

# Internal modules #

# First party modules #
from plumbing.autopaths import AutoPaths
from plumbing.cache     import property_cached

# Third party modules #
from shell_command import shell_call

###############################################################################
class Pandaseq(object):
    """Use Pandaseq to join read pairs together."""

    short_name = 'pandaseq'
    long_name  = 'PANDAseq Version 2.8.1'
    executable = 'pandaseq28'
    url        = 'https://github.com/neufeld/pandaseq/'
    license    = 'GPLv3'
    dependencies = []

    all_paths = """
    /assembled.fasta
    /unassembled.fastq
    """

    def __repr__(self): return '<%s object on %s>' % (self.__class__.__name__, self.pair)
    def __nonzero__(self): return bool(self.p.assembled)

    def __init__(self, pair, result_dir):
        # Save attributes #
        self.pair       = pair
        self.result_dir = result_dir
        # Auto paths #
        self.base_dir = self.result_dir + self.short_name + '/'
        self.p = AutoPaths(self.base_dir, self.all_paths)

    def run(self):
        command = 'pandaseq28 -T 1 -f %s -r %s -u %s -F 1> %s 2> %s'
        command = command % (self.fwd, self.rev,
                             self.p.unassembled,
                             self.p.assembled,
                             self.p.stderr)
        # Call it #
        shell_call(command)
        # Check #
        assert self.assembled.path.exists

    @property_cached
    def results(self):
        results = PandaseqResults(self)
        message = "You can't access results from PANDAseq before running the algorithm."
        if not results: raise Exception(message)
        return results

###############################################################################
class PandaseqResults(object):

    def __nonzero__(self): return bool(self.pandaseq.p.assembled)

    def __init__(self, pandaseq):
        self.pandaseq = pandaseq

###############################################################################
#sh.pandaseq28('-T', 1,        # direct output to file <f>, not stdout
#              '-f', self.fwd, # parsable table of per-sequence hits
#              '-r', self.rev, # set RNG seed to <n>
#              '-u', self.unassembled.path # unlimited ASCII text output line width
#              self.assembled, # prefer accessions over names in output
