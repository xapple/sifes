# Built-in modules #
import os
from collections import Counter

# Internal modules #
from sifes.taxonomy import Classify

# First party modules #
from plumbing.cache import property_cached

# Third party modules #
import sh

# Constants #
home = os.environ.get('HOME', '~') + '/'

###############################################################################
class QiimeClassify(Classify):
    """A wrapper to qiime assign_taxonomy.py
    You have to install qiime like this:
    
        $ cd ~/programs
        $ virtualenv qiime
        $ source qiime/bin/activate
        $ pip install numpy
        $ pip install qiime
        $ deactivate"""

    # Attributes #
    short_name    = 'qiime_classify'
    long_name     = 'Qiime Version 1.9.1'
    executable    = (home + 'programs/qiime/bin/python',
                     home + 'programs/qiime/bin/assign_taxonomy.py')
    doc           = 'http://qiime.org/scripts/assign_taxonomy.html'

    all_paths = """
    /stdout.txt
    /stderr.txt
    /centers_tax_assignments.log
    /centers_tax_assignments.txt
    """

    def __repr__(self): return '<%s object on %s>' % (self.__class__.__name__, self.centers)
    def __nonzero__(self): return bool(self.p.assignments_txt)

    def run(self, cpus=None):
        # Message #
        print self.message
        # It's a virtual env #
        qiime = sh.Command(self.executable[0])
        assert "1.9.1" in qiime(self.executable[1], '--version')
        # Output dir must exist #
        self.base_dir.create_if_not_exists()
        # Run #
        command = ("-i", self.centers,            # Path to the input fasta file
                   "-m", "uclust",                # Taxon assignment method
                   "-t", self.database.taxonomy,  # TSV mapping sequences to assigned taxonomy
                   "-r", self.database.alignment, # Path to reference sequences
                   "-o", self.base_dir)           # Directory to store result file
        qiime(self.executable[1], *command, _out=self.p.stdout.path, _err=self.p.stderr.path)

    @property_cached
    def results(self):
        results = QiimeClassifyResults(self)
        message = "You can't access results from QiimeClassify before running the algorithm."
        if not results: raise Exception(message)
        return results

###############################################################################
class QiimeClassifyResults(object):

    def __nonzero__(self): return bool(self.p.assignments_txt)
    def __len__(self):     return len(self.p.assignments_txt)

    @property
    def rank_names(self): return self.database.rank_names

    def __init__(self, qiime):
        # Attributes #
        self.qiime    = qiime
        self.database = qiime.database
        self.p        = qiime.p

    @property_cached
    def assignments(self):
        result = {}
        with open(self.p.assignments_txt, 'r') as handle:
            for line in handle:
                line = line.strip('\n')
                otu_name, species, confidence, hits = line.split('\t')
                result[otu_name] = tuple(i for i in species.split(';') if i)
        return result

    @property_cached
    def count_unassigned(self):
        """How many did not get a prediction at each level. Must be monotonically increasing."""
        result = [sum(1 for x in self.assignments.values() if x[0] == 'Unassigned')]
        for i, rank_name in enumerate(self.rank_names):
            if i == 0: continue
            result.append(sum(1 for x in self.assignments.values() if len(x) < i+1))
        return result

    @property_cached
    def count_assigned(self):
        """How many got a prediction at each level."""
        minus_unassigned = lambda u: len(self.assignments) - u
        return map(minus_unassigned, self.count_unassigned)

    @property
    def assignments_file(self):
        return self.p.assignments_txt