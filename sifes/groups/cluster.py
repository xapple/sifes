# Built-in modules #
import shutil

# Internal modules #
import illumitag
from plumbing.tmpstuff import new_temp_path
from plumbing.autopaths import AutoPaths
from fasta import FASTA
from illumitag.running.cluster_runner import ClusterRunner
from illumitag.clustering.otu.uparse import UparseOTUs
from illumitag.clustering.otu.uclust import UclustOTUs
from illumitag.clustering.otu.cdhit import CdhitOTUs
from illumitag.clustering.reporting import ClusterReporter
from illumitag.reporting.clusters import ClusterReport

# Third party modules #
import pandas
from tqdm import tqdm
from shell_command import shell_output

###############################################################################
class Cluster(object):
    """Analyzes a group of samples."""

    all_paths = """
    /reads/all_reads.fasta
    /otus/
    /logs/
    /report/report.pdf
    /metadata.csv
    """

    def __repr__(self): return '<%s object "%s" with %i samples>' % (self.__class__.__name__, self.name, len(self.samples))
    def __iter__(self): return iter(self.samples)
    def __len__(self): return len(self.samples)
    def __getitem__(self, key):
        if isinstance(key, basestring): return [c for c in self.children if c.short_name == key.lower()][0]
        elif isinstance(key, int) and hasattr(self.first, 'num'): return [c for c in self.children if c.num == key][0]
        else: return self.children[key]

    @property
    def first(self): return self.children[0]

    @property
    def count_seq(self):
        return sum([len(sample) for sample in self])

    def __init__(self, samples, name, base_dir=None):
        # Save samples #
        self.name = name
        self.samples, self.children = samples, samples
        # Check names are unique #
        names = [s.short_name for s in samples if s.used]
        assert len(names) == len(set(names))
        # Figure out pools #
        self.pools = list(set([s.pool for s in self.samples]))
        self.pools.sort(key = lambda x: x.id_name)
        # Directory #
        if base_dir: self.base_dir = base_dir
        else: self.base_dir = illumitag.view_dir + "clusters/" + self.name + '/'
        # Loaded #
        self.loaded = False

    def load(self):
        """A second __init__ that is delayed and called only if needed"""
        # Load the pools and samples #
        for p in self.pools: p.load()
        for s in self.samples: s.load()
        # Dir #
        self.p = AutoPaths(self.base_dir, self.all_paths)
        # Figure out if it's a project #
        if set(self.samples) == set(self.first.pool.project.samples): self.project = self.first.pool.project
        else: self.project = None
        # Runner #
        self.runner = ClusterRunner(self)
        # FASTA #
        self.reads = FASTA(self.p.all_reads_fasta)
        # OTU picking #
        self.otu_uparse = UparseOTUs(self)
        self.otu_uclust = UclustOTUs(self)
        self.otu_cdhit  = CdhitOTUs(self)
        # Preferred #
        self.otus = self.otu_uparse
        # Simple reporting #
        self.reporter = ClusterReporter(self)
        # Full report #
        self.report = ClusterReport(self)
        # Loaded #
        self.loaded = True
        # Return self for convenience #
        return self

    def run(self, *args, **kwargs):
        self.runner.run(*args, **kwargs)

    def run_slurm(self, *args, **kwargs):
        self.runner.run_slurm(*args, **kwargs)

    def process_samples(self):
        for sample in tqdm(self): sample.process()

    def combine_reads(self):
        """This is the first function should call. It will combine all the
        reads of all the samples of this cluster into one big FASTA file."""
        paths = [sample.fasta.path for sample in self]
        shell_output('cat %s > %s' % (' '.join(paths), self.reads))
        return self.reads

    def set_size(self, length):
        """Trim all sequences to a specific length starting from the end."""
        self.size_trimmed = FASTA(new_temp_path())
        def trim_iterator(reads):
            for read in reads:
                if len(read) < length: continue
                yield read[-length:]
        self.size_trimmed.write(trim_iterator(self.reads))
        self.size_trimmed.close()
        # Replace it #
        self.reads.remove()
        shutil.move(self.size_trimmed, self.reads)

    def run_uparse(self): self.otu_uparse.run()

    @property
    def metadata(self):
        return pandas.DataFrame([s.info for s in self], index=[s.short_name for s in self])

    def export_metadata(self):
        self.metadata.to_csv(self.p.metadata, sep='\t', encoding='utf-8')