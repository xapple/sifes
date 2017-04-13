# Built-in modules #
import os, multiprocessing, glob

# Internal modules #

# First party modules #
from plumbing.autopaths import AutoPaths
from plumbing.cache     import property_cached
from plumbing.autopaths import FilePath

# Different databases #
from seqsearch.databases.silva_mothur import silva_mothur
from seqsearch.databases.foraminifera import foraminifera

# Third party modules #
import sh

###############################################################################
class MothurClassify(object):
    """A wrapper to mothur classify"""

    # Attributes #
    short_name    = 'mothur_classify'
    long_name     = 'Mothur Version 1.37.4'
    executable    = 'mothur'
    doc           = 'http://www.mothur.org/wiki/Classify.seqs'
    database_name = 'the non-redundant, no-gaps Silva version 123 database'

    # Parameters #
    default_database = 'silva'
    bootstrap_cutoff = 70

    all_paths = """
    /stdout.txt
    /stderr.txt
    /centers.fasta
    /assignments.txt
    /summary.txt
    /flipped.txt
    """

    def __repr__(self): return '<%s object on %s>' % (self.__class__.__name__, self.pair)
    def __nonzero__(self): return bool(self.p.assignments)

    def __init__(self, centers, result_dir, database=None):
        # Attributes #
        self.centers    = centers
        self.result_dir = result_dir
        # Auto paths #
        self.base_dir = self.result_dir + self.short_name + '/'
        self.p = AutoPaths(self.base_dir, self.all_paths)
        # Default or user specified #
        if database is None: self.database = self.default_database
        else:                self.database = database
        # The database to use #
        if   self.database == 'silva':        self.database = silva_mothur
        elif self.database == 'foraminifera': self.database = foraminifera
        else: raise NotImplemented()

    def run(self, cpus=None, bootstrap_cutoff=None):
        # Message #
        print "Classifying file '%s'" % self.centers
        # Number of cores #
        if cpus is None: cpus = min(multiprocessing.cpu_count(), 32)
        # Bootstrap cutoff #
        cutoff = self.bootstrap_cutoff if bootstrap_cutoff is None else bootstrap_cutoff
        # Prepare #
        self.p.centers.link_from(self.centers)
        current_dir = os.getcwd()
        os.chdir(self.base_dir)
        # Clean #
        self.p.assignments.remove()
        self.p.summary.remove()
        self.p.flipped.remove()
        for p in glob.glob('mothur.*.logfile'): os.remove(p)
        # Run #
        command = ("#classify.seqs(", # The command
                   " fasta=%s,"     , # The input file
                   " reference=%s," , # The database (was also called 'template')
                   " taxonomy=%s,"  , # The taxonomy file
                   " processors=%s,", # The number of threads
                   " cutoff=%s,"    , # At 0 you get a full taxonomy for every sequence.
                   " probs=F);")      # Disable the output of bootstrap values with your taxonomy
        command = ''.join(command)
        sh.mothur(command % (self.p.centers, self.database.alignment, self.database.taxonomy, cpus, cutoff),
                             _out=self.p.stdout.path, _err=self.p.stderr.path)
        # Check output #
        if "ERROR" in self.p.stdout.contents:
            raise Exception("Mothur classify didn't run correctly.")
        # Back #
        os.chdir(current_dir)
        # Outputs #
        self.assignments = FilePath(self.base_dir + "centers.%s.wang.taxonomy"    % self.database.nickname)
        self.summary     = FilePath(self.base_dir + "centers.%s.wang.tax.summary" % self.database.nickname)
        self.flipped     = FilePath(self.base_dir + "centers.%s.wang.flip.accnos" % self.database.nickname)
        # Move #
        self.assignments.move_to(self.p.assignments)
        self.summary.move_to(self.p.summary)
        if self.flipped: self.flipped.move_to(self.p.flipped)
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
    """Here are examples outputs in assignments:
         'OTU-3576': ('Bacteria',
                      'Bacteria_unclassified',
                      'Bacteria_unclassified',
                      'Bacteria_unclassified',
                      'Bacteria_unclassified',
                      'Bacteria_unclassified',
                      '')
         'OTU-3577': ('Bacteria',
                      'Firmicutes',
                      'Clostridia',
                      'Clostridiales',
                      'Clostridiaceae_1',
                      'Clostridiaceae_1_unclassified',
                      '')
         'OTU-4086': ('unknown',
                      'unknown_unclassified',
                      'unknown_unclassified',
                      'unknown_unclassified',
                      'unknown_unclassified',
                      '')
        '"""

    def __nonzero__(self): return bool(self.p.stdout)
    def __len__(self):     return len(self.p.assignments)

    def __init__(self, mothur):
        # Attributes #
        self.mothur = mothur
        self.p      = mothur.p

    @property_cached
    def assignments(self):
        result = {}
        with open(self.p.assignments, 'r') as handle:
            for line in handle:
                otu_name, species = line.split('\t')
                species           = [i.strip('\n') for i in species.split(';')]
                result[otu_name]  = tuple(species)
        return result

    @property
    def rank_names(self):
        """The names of the ranks. Mothur skips Kingdom."""
        return ['Domain',
                'Phylum',
                'Class',
                'Order',
                'Family',
                'Genus',
                'Species']

    @property_cached
    def count_unassigned(self):
        """How many did not get a prediction at each level"""
        return [sum((1 for x in self.assignments.values() if x[0] == 'unknown')),      # Domain
                sum((1 for x in self.assignments.values() if 'unclassified' in x[1])), # Phylum
                sum((1 for x in self.assignments.values() if 'unclassified' in x[2])), # Class
                sum((1 for x in self.assignments.values() if 'unclassified' in x[3])), # Order
                sum((1 for x in self.assignments.values() if 'unclassified' in x[4])), # Family
                sum((1 for x in self.assignments.values() if 'unclassified' in x[5])), # Genus
                sum((1 for x in self.assignments.values() if x[6] == '')),             # Species
                ]

    @property_cached
    def count_assigned(self):
        """How many got a prediction at each level"""
        minus_unassigned = lambda u: len(self.assignments) - u
        return map(minus_unassigned, self.count_unassigned)

