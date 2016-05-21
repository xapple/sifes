# Built-in modules #

# Internal modules #
import os, multiprocessing

# First party modules #
from plumbing.autopaths import AutoPaths
from plumbing.cache     import property_cached
from plumbing.autopaths import FilePath
from fasta import FASTA

# Third party modules #
import sh

###############################################################################
class MothurJoin(object):
    """A wrapper to mothur make contigs"""

    short_name = 'mothur_join'
    long_name  = 'Mothur Version 1.37.4'
    executable = 'mothur'

    all_paths = """
    /fwd.fastq
    /rev.fastq
    /stdout.txt
    /stderr.txt
    /assembled.fasta
    /unassembled.fasta"""
   #/fwd.trim.contigs.fasta
   #/fwd.contigs.report
   #/fwd.scrap.contigs.fasta

    def __repr__(self): return '<%s object on %s>' % (self.__class__.__name__, self.pair)
    def __nonzero__(self): return bool(self.p.assembled)

    def __init__(self, pair, result_dir, sample_name):
        # Save attributes #
        self.pair        = pair
        self.result_dir  = result_dir
        self.sample_name = sample_name
        # Auto paths #
        self.base_dir = self.result_dir + self.short_name + '/'
        self.p = AutoPaths(self.base_dir, self.all_paths)

    def run(self, cpus=1):
        """The make contigs command.
        http://www.mothur.org/wiki/Make.contigs"""
        # Message #
        print "Joining sample '%s'" % self.sample_name
        # Number of cores #
        if cpus is None: cpus = min(multiprocessing.cpu_count(), 32)
        # Check input #
        assert hasattr(self.pair, 'fwd') and hasattr(self.pair, 'rev')
        # Prepare #
        self.p.fwd.link_from(self.pair.fwd)
        self.p.rev.link_from(self.pair.rev)
        current_dir = os.getcwd()
        os.chdir(self.base_dir)
        # Run #
        command = "#make.contigs(ffastq=%s, rfastq=%s, processors=%s);"
        sh.mothur(command % (self.p.fwd, self.p.rev, cpus),
                  _out=self.p.stdout.path,_err=self.p.stderr.path)
        # Check output #
        if "ERROR" in self.p.stdout.contents:
            raise Exception("Mothur make contigs didn't run correctly.")
        # Back #
        os.chdir(current_dir)
        # Check #
        assert self.results.joined
        # Link #
        self.p.assembled.link_from(self.results.joined)
        self.p.unassembled.link_from(self.results.scrap)
        # Return #
        return self.results.assembled

    @property_cached
    def results(self):
        results = MothurJoinResults(self)
        message = "You can't access results from MothurJoin before running the algorithm."
        if not results: raise Exception(message)
        return results

###############################################################################
class MothurJoinResults(object):

    def __nonzero__(self): return bool(self.p.stdout)
    def __len__(self):     return len(self.p.assembled)

    def __init__(self, parent):
        # Attributes #
        self.parent      = parent
        self.base_dir    = parent.base_dir
        self.p           = parent.p
        self.joined      = FASTA(   self.base_dir + 'fwd.trim.contigs.fasta')
        self.report      = FilePath(self.base_dir + 'fwd.contigs.report')
        self.scrap       = FASTA(   self.base_dir + 'fwd.scrap.contigs.fasta')
        self.outputs     = (self.joined, self.report, self.scrap)
        self.assembled   = FASTA(self.p.assembled)
        self.unassembled = FASTA(self.p.unassembled)

    @property
    def unassembled_count(self):
        return len(self.parent.pair) - len(self.assembled)

    @property
    def unassembled_percent(self):
        percent = (self.unassembled_count / len(self.parent.pair)) * 100.0
        return "%.1f%%" % percent
