# Built-in modules #
import os, shutil, multiprocessing

# Internal modules #

# First party modules #
from plumbing.autopaths import AutoPaths
from plumbing.common import which
from plumbing.cache import property_cached
from plumbing.csv_tables import TSVTable

# Third party modules #
import sh

###############################################################################
class Crest(object):

    short_name = 'crest'
    long_name  = 'LCAClassifier/CREST version 2.0.4'
    executable = 'classify'
    url        = "https://github.com/lanzen/CREST"
    article    = "http://dx.plos.org/10.1371/journal.pone.0049334"

    all_paths = """
    /otu_table.csv
    /otu_table_norm.csv
    /centers.fasta
    /graphs/
    /stats/
    /comp_phyla/
    /comp_tips/
    /comp_order/
    /comp_class/
    /blast/db_hits.xml
    /blast/stdout.txt
    /blast/stderr.txt
    /crest/assignments.tsv
    /crest/composition.tsv
    /crest/tree.txt
    /crest/Relative_Abundance.tsv
    /crest/Richness.tsv
    /crest/stdout.txt
    /crest/stderr.txt
    """

    def __nonzero__(self): return bool(self.p.assignments)

    def __init__(self, centers, result_dir):
        # Parent #
        self.centers    = centers
        self.database   = database
        self.result_dir = result_dir
        # Auto paths #
        self.base_dir = self.result_dir + self.short_name + '/'
        self.p = AutoPaths(self.base_dir, self.all_paths)
        # Find where the database is on the file system #
        self.database_path  = which('classify').physical_path.directory.directory
        self.database_path += 'parts/flatdb/%s/%s.fasta' % (database, database)
        # Graphs #
        self.graphs = [getattr(plots, cls_name)(self) for cls_name in plots.__all__[:-1]]
        # OTU table #
        self.otu_csv      = TSVTable(self.p.otu_csv)
        self.otu_csv_norm = TSVTable(self.p.otu_csv_norm)

    def run(self, cpus=None):
        # Number of cores #
        if cpus is None: cpus = min(multiprocessing.cpu_count(), 32)
        # Run #
        sh.blastn('-task',              'megablast',
                  '-num_threads',       cpus,
                  '-query',             self.centers,
                  '-db',                self.database_path,
                  '-out',               self.p.db_hits,
                  '-max_target_seqs',   '100',
                  '-outfmt',            '5',
                  _out=self.p.blast_stdout, _err=self.p.blast_stderr)
        # Check #
        if os.path.getsize(self.p.db_hits) == 0:
            raise Exception("Hits file empty. The MEGABLAST process was probably killed.")
        # CREST #
        self.p.crest_dir.remove()
        sh.classify('--verbose', '--rdp',
                    '-o', self.base_dir + 'crest/',
                    '-d', self.database,
                    self.p.db_hits,
                    out=self.p.crest_stdout.path, _err=self.p.crest_stderr.path)
        # Move #
        shutil.move(self.p.db_hits.prefix_path + '_Composition.tsv', self.p.crest_composition)
        shutil.move(self.p.db_hits.prefix_path + '_Tree.txt', self.p.crest_tree)
        shutil.move(self.p.db_hits.prefix_path + '_Assignments.tsv', self.p.crest_assignments)
        # Clean up #
        if os.path.exists("error.log") and os.path.getsize("error.log") == 0: os.remove("error.log")

    @property_cached
    def results(self):
        results = CrestResults(self)
        message = "You can't access results from Crest before running the algorithm."
        if not results: raise Exception(message)
        return results

###############################################################################
class CrestResults(object):

    def __nonzero__(self): return bool(self.p.assignments)
    def __len__(self):     return len(self.p.assembled)

    def __init__(self, crest):
        self.crest  = crest
        self.p      = crest.p

    @property_cached
    def assignments(self):
        result = {}
        with open(self.p.assignments, 'r') as handle:
            for line in handle:
                code, species = line.split('\t')
                result[code] = tuple(species.strip('\n').split(';'))[:8]
        return result

    @property
    def count_assigned(self):
        """How many got a position"""
        return len([s for s in self.assignments.values() if s != ('No hits',)])