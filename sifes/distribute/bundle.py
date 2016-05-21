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
    """Regroup result files and reports from one of several projects for delivery."""

    all_paths = """
    /projects/
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
        for p in self.projects:
            # Directories #
            proj_dir    = DirectoryPath(self.p.projects_dir + p.name)
            # Reports #
            reports_dir = DirectoryPath(proj_dir + 'reports')
            reports_dir.create(safe=True)
            p.cluster.report.output_path.copy(reports_dir)
            for s in p: s.report.output_path.copy(reports_dir + s.short_name)
            # Data files #
            data_dir = DirectoryPath(proj_dir + 'data')
            data_dir.create(safe=True)
            p.cluster.centering.results.centers.copy(data_dir) # centers
            p.cluster.taxonomy.p.assignments.copy(data_dir)    # assignments
            p.cluster.otu_table.p.flat.copy(data_dir)          # otu_table
            # Taxa tables #
            taxa_dir = DirectoryPath(data_dir + 'taxa_tables')
            taxa_dir.remove()
            p.cluster.taxa_table.base_dir.copy(taxa_dir)

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
