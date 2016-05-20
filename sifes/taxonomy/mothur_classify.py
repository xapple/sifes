# Built-in modules #
import os, multiprocessing

# Internal modules #
import sifes

# First party modules #
from plumbing.autopaths import AutoPaths
from plumbing.cache     import property_cached
from plumbing.autopaths import FilePath

# Third party modules #
import sh

# Constants #
silva_path    = sifes.home + 'databases/silva_mothur_v123/silva.nr_v123.fasta'
taxonomy_path = sifes.home + 'databases/silva_mothur_v123/silva.nr_v123.tax'

###############################################################################
class MothurClassify(object):
    """A wrapper to mothur classify"""

    # Attributes #
    short_name = 'mothur_classify'
    long_name  = 'Mothur Version 1.37.4'
    executable = 'mothur'
    doc        = 'http://www.mothur.org/wiki/Classify.seqs'
    database   = 'the non-redundant, no-gaps Silva v123 database'

    # Parameters #
    bootstrap_cutoff = 60

    all_paths = """
    /stdout.txt
    /stderr.txt
    /centers.fasta
    /assignments.txt
    /summary.txt
    """

    def __repr__(self): return '<%s object on %s>' % (self.__class__.__name__, self.pair)
    def __nonzero__(self): return bool(self.p.assignments)

    def __init__(self, centers, database, result_dir):
        # Attributes #
        self.centers    = centers
        self.database   = database
        self.result_dir = result_dir
        # Auto paths #
        self.base_dir = self.result_dir + self.short_name + '/'
        self.p = AutoPaths(self.base_dir, self.all_paths)

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
        # Run #
        command = "#classify.seqs(fasta=%s, template=%s, taxonomy=%s, processors=%s, cutoff=%s, probs=F);"
        sh.mothur(command % (self.p.centers, silva_path, taxonomy_path, cpus, cutoff),
                  _out=self.p.stdout.path,_err=self.p.stderr.path)
        # Check output #
        if "ERROR" in self.p.stdout.contents:
            raise Exception("Mothur classify didn't run correctly.")
        # Back #
        os.chdir(current_dir)
        # Outputs #
        self.assignments = FilePath(self.base_dir + "centers.nr_v123.wang.taxonomy")
        self.summary     = FilePath(self.base_dir + "centers.nr_v123.wang.tax.summary")
        # Move #
        self.p.assignments.remove()
        self.p.summary.remove()
        self.assignments.move_to(self.p.assignments)
        self.summary.move_to(self.p.summary)
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

