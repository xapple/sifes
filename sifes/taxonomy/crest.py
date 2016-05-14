# Built-in modules #
import os, shutil, multiprocessing

# Internal modules #
from illumitag.clustering.statistics import StatsOnTaxonomy
from illumitag.clustering.taxonomy import Taxonomy, SimpleTaxonomy
from illumitag.clustering.taxonomy import plots
from illumitag.clustering.composition.custom_rank import CompositionPhyla, CompositionOrder, CompositionClass
from illumitag.clustering.composition.custom_rank import CompositionFamily, CompositionGenus
from illumitag.clustering.composition.tips import CompositionTips

# First party modules #
from fasta import FASTA
from plumbing.autopaths import AutoPaths
from plumbing.slurm import num_processors
from plumbing.common import which
from plumbing.cache import property_cached
from plumbing.csv_tables import TSVTable

# Third party modules #
import sh

# Constants #
home = os.environ['HOME'] + '/'

###############################################################################
class CrestTaxonomy(Taxonomy):
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
    /db_hits.xml
    /crest/assignments.tsv
    /crest/composition.tsv
    /crest/tree.txt
    /crest/Relative_Abundance.tsv
    /crest/Richness.tsv
    /crest_stdout.txt
    """

    short_name = 'crest'
    title = 'LCAClassifier'
    article = "http://dx.plos.org/10.1371/journal.pone.0049334"
    version = "version 2.0 - March 2014"

    def __init__(self, fasta_path, parent, database='silvamod', base_dir=None):
        # Parent #
        self.otu, self.parent = parent, parent
        # Inherited #
        self.samples = self.parent.samples
        # FASTA #
        self.fasta = FASTA(fasta_path)
        # The database to use #
        self.database = database
        # Find where it is on the file system
        self.database_path  = which('classify').physical_path.directory.directory
        self.database_path += 'parts/flatdb/%s/%s.fasta' % (database, database)
        # Dir #
        if base_dir is None: self.base_dir = self.parent.p.crest_dir
        else: self.base_dir = base_dir
        self.p = AutoPaths(self.base_dir, self.all_paths)
        # Graphs #
        self.graphs = [getattr(plots, cls_name)(self) for cls_name in plots.__all__[:-1]]
        # OTU table #
        self.otu_csv = TSVTable(self.p.otu_csv)
        self.otu_csv_norm = TSVTable(self.p.otu_csv_norm)
        # Filtered centers file #
        self.centers = FASTA(self.p.centers)
        # Composition tables #
        self.comp_phyla  = CompositionPhyla(self,  self.p.comp_phyla)
        self.comp_order  = CompositionOrder(self,  self.p.comp_order)
        self.comp_class  = CompositionClass(self,  self.p.comp_class)
        self.comp_family = CompositionFamily(self, self.p.comp_order)
        self.comp_genus  = CompositionGenus(self,  self.p.comp_class)
        # The complete one #
        self.comp_tips  = CompositionTips(self,  self.p.comp_tips)
        # Stats #
        self.stats = StatsOnTaxonomy(self)

    def assign(self, cpus=None):
        # Number of cores #
        if cpus is None: cpus = min(multiprocessing.cpu_count(), 32)
        # Run #
        sh.blastn('-task', 'megablast', '-num_threads', cpus, '-query', self.fasta, '-db', self.database_path, '-max_target_seqs', '100', '-outfmt', '5' ,'-out', self.p.db_hits)
        if os.path.getsize(self.p.db_hits) == 0: raise Exception("Hits file empty. The MEGABLAST process was probably killed.")
        # CREST #
        self.p.crest_dir.remove()
        sh.classify('--verbose', '--rdp', '-o', self.base_dir + 'crest/', '-d', self.database, self.p.db_hits, _out=self.p.stdout.path)
        shutil.move(self.p.db_hits.prefix_path + '_Composition.tsv', self.p.crest_composition)
        shutil.move(self.p.db_hits.prefix_path + '_Tree.txt', self.p.crest_tree)
        shutil.move(self.p.db_hits.prefix_path + '_Assignments.tsv', self.p.crest_assignments)
        # Clean up #
        if os.path.exists("error.log") and os.path.getsize("error.log") == 0: os.remove("error.log")
        # Return #
        return self

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

###############################################################################
class SimpleCrestTaxonomy(SimpleTaxonomy, CrestTaxonomy):
    short_name = 'crest'

    def __init__(self, fasta, base_dir, database='silvamod'):
        SimpleTaxonomy.__init__(self, fasta, base_dir)
        # The database to use #
        self.database = database
        # Find where it is on the file system
        self.database_path  = which('classify').physical_path.directory.directory
        self.database_path += 'parts/flatdb/%s/%s.fasta' % (database, database)
