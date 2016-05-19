# Built-in modules #
import os, multiprocessing

# Internal modules #
import sifes

# First party modules #
from plumbing.autopaths import AutoPaths
from plumbing.cache     import property_cached
from plumbing.autopaths import FilePath
from fasta import FASTA

# Third party modules #
import sh

# Constants #
silva_path    = sifes.home + 'databases/silva_mothur_v123/silva.nr_v123.fasta'
taxonomy_path = sifes.home + 'databases/silva_mothur_v123/silva.nr_v123.tax'

###############################################################################
class MothurClassify(object):
    """A wrapper to mothur classify"""

    short_name = 'mothur_classify'
    long_name  = 'Mothur Version 1.37.4'
    executable = 'mothur'
    doc        = 'http://www.mothur.org/wiki/Classify.seqs'

    all_paths = """
    /stdout.txt
    /stderr.txt
    /centers.fasta
    /assignment.txt
    """

    def __repr__(self): return '<%s object on %s>' % (self.__class__.__name__, self.pair)
    def __nonzero__(self): return bool(self.p.assignment)

    def __init__(self, centers, database, result_dir):
        # Attributes #
        self.centers    = centers
        self.database   = database
        self.result_dir = result_dir
        # Auto paths #
        self.base_dir = self.result_dir + self.short_name + '/'
        self.p = AutoPaths(self.base_dir, self.all_paths)

    def run(self, cpus=None):
        # Message #
        print "Classifying file '%s'" % self.centers
        # Number of cores #
        if cpus is None: cpus = min(multiprocessing.cpu_count(), 32)
        # Prepare #
        self.p.centers.link_from(self.centers)
        current_dir = os.getcwd()
        os.chdir(self.base_dir)
        # Run #
        command = "#classify.seqs(fasta=%s, template=%s, taxonomy=%s, processors=%s);"
        sh.mothur(command % (self.p.centers, silva_path, taxonomy_path, cpus),
                  _out=self.p.stdout.path,_err=self.p.stderr.path)
        # Check output #
        if "ERROR" in self.p.stdout.contents:
            raise Exception("Mothur classify didn't run correctly.")
        # Back #
        os.chdir(current_dir)
        # Check #
        pass

    def clean(self):
        self.p.stderr.remove()
        self.p.assembled.remove()

    @property_cached
    def results(self):
        results = MothurClassifyResults(self)
        message = "You can't access results from MothurClassify before running the algorithm."
        if not results: raise Exception(message)
        return results

###############################################################################
class MothurClassifyResults(object):

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
