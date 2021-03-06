# Built-in modules #

# Internal modules #
import sifes
from sifes.groups.aggregate import Aggregate

# First party modules #
from plumbing.cache     import property_cached
from plumbing.autopaths import DirectoryPath

# Third party modules #

###############################################################################
class Bundle(Aggregate):
    """Regroup result files and reports from one or several projects for delivery."""

    all_paths = """
    /projects/
    /samples/
    /metadata/samples.xlsx
    /metadata/multiplexed.xlsx
    /demultiplexing/report.pdf
    """

    def __init__(self, name, samples, out_dir=None):
        # Directory #
        if out_dir is None: out_dir = sifes.bundles_dir
        # Super #
        super(self.__class__,self).__init__(name, samples, out_dir)
        # Figure out the projects within #
        proj_names = sorted(list(set([s.project_short_name for s in samples])))
        self.projects = [sifes.projects[p] for p in proj_names]

    def run(self):
        self.base_dir.remove()
        self.base_dir.create()
        for p in self.projects:
            # Directories #
            proj_dir    = DirectoryPath(self.p.projects_dir + p.name)
            # Reports for samples #
            reports_dir = DirectoryPath(proj_dir + 'reports')
            reports_dir.create(safe=True)
            for s in p: s.report.output_path.copy(reports_dir + s.short_name + '.pdf')
            # Early exits #
            if not all(s.filter for s in p): return
            if not p.cluster: return
            # Reports for cluster #
            p.cluster.report.output_path.copy(reports_dir + 'project_report.pdf')
            # Data files #
            data_dir = DirectoryPath(proj_dir + 'data')
            data_dir.create(safe=True)
            p.cluster.centering.results.centers.copy(data_dir + 'centers.fasta')           # centers
            p.cluster.taxonomy.results.assignments_file.copy(data_dir + 'assignments.txt') # assignments
            p.cluster.otu_table.p.flat.copy(data_dir + 'otu_table.tsv')                    # otu_table
            # Optional extras #
            #p.cluster.reads.link_to(data_dir, safe=True)             # all_reads.fasta
            #p.cluster.centering.readmap.link_to(data_dir, safe=True) # readmap.uc
            # Taxa tables #
            taxa_dir = DirectoryPath(data_dir + 'taxa_tables')
            taxa_dir.remove()
            p.cluster.taxa_table.base_dir.copy(taxa_dir)

    def add_samples(self):
        """
        In some cases, for instance when the samples are originally multiplexed in
        FASTQ files, it can be useful to include the original raw FASTQ files of
        each sample to the bundle.
        """
        for s in self:
            # Directories #
            sample_dir = DirectoryPath(self.p.samples_dir + s.short_name)
            sample_dir.create(safe=True)
            # Raw reads for each sample #
            s.pair.fwd.copy(sample_dir)
            s.pair.rev.copy(sample_dir)

    @property_cached
    def results(self):
        results = BundleResults(self)
        message = "You can't access results from a bundle before running the algorithm."
        if not results: raise Exception(message)
        return results

###############################################################################
class BundleResults(object):

    def __nonzero__(self): return bool(self.p.x)

    def __init__(self, parent):
        self.parent      = parent
        self.base_dir    = parent.base_dir
        self.p           = parent.p
