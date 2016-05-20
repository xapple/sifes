# Built-in modules #
import multiprocessing

# Internal modules #
import sifes

# First party modules #
from plumbing.autopaths import AutoPaths
from plumbing.cache import property_cached

# Third party modules #
import sh

###############################################################################
class Rdp(object):

    short_name = 'rdp'
    long_name  = 'RDP Classifier version 2.11'
    executable = 'rdp_multiclassifier'
    url        = "https://github.com/rdpstaff/classifier"
    article    = "http://www.ncbi.nlm.nih.gov/pmc/articles/PMC1950982/"
    exec_path  = sifes.home + "programs/rdp_classifier/dist/classifier.jar"

    all_paths = """
    /assignments.txt
    /composition.txt
    /rdp_stdout.txt
    /rdp_stderr.txt
    """

    def __nonzero__(self): return bool(self.p.assignments)

    def __init__(self, centers, database, result_dir):
        # Attributes #
        self.centers    = centers
        self.database   = database
        self.result_dir = result_dir
        # Auto paths #
        self.base_dir = self.result_dir + self.short_name + '/'
        self.p = AutoPaths(self.base_dir, self.all_paths)

    def run(self, cpus=None):
        # Number of cores #
        if cpus is None: cpus = min(multiprocessing.cpu_count(), 32)
        # Run #
        sh.java('-Xmx10g',
                '-jar',             self.exec_path,
                'classify',
                '--conf='           + '0.5',
                '--format='         + 'fixrank',
                '--hier_outfile='   + self.p.composition,
                '--outputFile='     + self.p.assignments,
                self.centers,
                _out=self.p.rdp_stdout.path, _err=self.p.rdp_stderr.path)
        # Check #
        if self.p.rdp_stdout.contents.startswith('Command Error'):
            raise Exception("RDP command error. See stdout.")
        if not self.p.assignments:
            raise Exception("Assignments file empty. The RDP process probably didn't run.")

    @property_cached
    def results(self):
        results = RdpResults(self)
        message = "You can't access results from RDP before running the algorithm."
        if not results: raise Exception(message)
        return results

###############################################################################
class RdpResults(object):

    def __nonzero__(self): return bool(self.p.assignments)
    def __len__(self):     return len(self.p.assembled)

    def __init__(self, rdp):
        self.rdp = rdp
        self.p   = rdp.p

    @property_cached
    def assignments(self):
        result = {}
        with open(self.p.assignments, 'r') as handle:
            for line in handle:
                line = line.strip('\n').split('\t')
                code = line.pop(0)
                species = tuple(x.strip('"') for x in line[1::3])
                result[code] = species
        return result

    @property
    def count_assigned(self):
        """How many got a position"""
        return len([s for s in self.assignments.values() if s != ('No hits',)])