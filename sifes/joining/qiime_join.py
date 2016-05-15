# Built-in modules #

# Internal modules #
import sifes

# First party modules #
from fasta              import FASTA
from plumbing.autopaths import AutoPaths
from plumbing.cache     import property_cached

# Third party modules #
import sh

###############################################################################
class QiimeJoin(object):
    """A wrapper to join_paired_ends.py"""

    short_name = 'qiime_join'
    long_name  = 'QIIME Version 1.9.1'
    executable = 'join_paired_ends.py'

    all_paths = """
    /assembled.fasta
    /stderr.out
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
        command = sh.Command(self.executable)
        env = {'HOME': sifes.home + 'no_home', 'PYTHONPATH':''}
        command('-f', self.pair.fwd,
                '-r', self.pair.rev,
                '-o', self.p.assembled,
                '--pe_join_method', 'SeqPrep',
                _env = env, # Requires outdated scikit-bio
                _err = self.p.stderr.path)
        # Check #
        assert self.p.assembled

    def clean(self):
        self.p.stderr.remove()
        self.p.assembled.remove()

    @property_cached
    def results(self):
        results = QiimeJoinResults(self)
        message = "You can't access results from QiimeJoin before running the algorithm."
        if not results: raise Exception(message)
        return results

###############################################################################
class QiimeJoinResults(object):

    def __nonzero__(self): return bool(self.p.assembled)
    def __len__(self):     return len(self.p.assembled)

    def __init__(self, parent):
        # Attributes #
        self.parent      = parent
        self.p           = parent.p
        self.pair        = parent.pair
        self.assembled   = FASTA(self.p.assembled)
