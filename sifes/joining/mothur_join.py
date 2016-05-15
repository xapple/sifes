# Built-in modules #

# Internal modules #
import os, shutil

# First party modules #
from plumbing.autopaths import AutoPaths
from plumbing.cache     import property_cached
from plumbing.tmpstuff import new_temp_dir
from plumbing.autopaths import FilePath
from fasta import FASTA, FASTQ
from fasta.common import add_dummy_scores

# Third party modules #
import sh

###############################################################################
class MothurJoin(object):
    """A wrapper to mothur make contigs"""

    short_name = 'mothur_join'
    long_name  = 'Mothur Version 1.37.4'
    executable = 'mothur'

    all_paths = """
    /fwd.fastq.gz
    /rev.fastq.gz
    /stdout.txt
    /stderr.txt
    /assembled.fasta"""
   #/fwd.trim.contigs.fasta
   #/fwd.contigs.report
   #/fwd.scrap.contigs.fasta

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
        """The make contigs command.
        http://www.mothur.org/wiki/Make.contigs"""
        # Check input #
        assert hasattr(self.pair, 'fwd') and hasattr(self.pair, 'rev')
        # Prepare #
        self.p.fwd.link_from(self.pair.fwd)
        self.p.rev.link_from(self.pair.rev)
        current_dir = os.getcwd()
        os.chdir(self.base_dir)
        # Run #
        sh.mothur("#make.contigs(ffastq=%s, rfastq=%s);" % (self.p.fwd, self.p.rev),
                  _out=self.p.stdout.path,_err=self.p.stderr.path)
        # Check output #
        if "ERROR" in self.p.stdout.contents:
            raise Exception("Mothur make contigs didn't run correctly.")
        # Back #
        os.chdir(current_dir)
        # Check #
        self.p.assembled.link_from(self.results.joined)
        assert self.results.joined

    def clean(self):
        self.p.stderr.remove()
        self.p.assembled.remove()

    @property_cached
    def results(self):
        results = MothurJoinResults(self)
        message = "You can't access results from MothurJoin before running the algorithm."
        if not results: raise Exception(message)
        return results

###############################################################################
class MothurJoinResults(object):

    def __nonzero__(self): return bool(self.p.assembled)
    def __len__(self):     return len(self.p.assembled)

    def __init__(self, parent):
        # Attributes #
        self.parent      = parent
        self.p           = parent.p
        self.pair        = parent.pair
        self.joined      = FASTA(   self.p.fwd.prefix_path + '.trim.contigs.fasta')
        self.report      = FilePath(self.p.fwd.prefix_path + '.contigs.report')
        self.scrap       = FASTA(   self.p.fwd.prefix_path + '.scrap.contigs.fasta')
        self.outputs     = (self.joined, self.report, self.scrap)
        self.assembled   = FASTA(self.p.assembled)
