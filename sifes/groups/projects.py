# Built-in modules #

# Internal modules #
import illumitag
from illumitag.groups.aggregate import Collection, Aggregate
from plumbing.autopaths import AutoPaths

# Third party modules #

###############################################################################
class Projects(Collection):
    """A collection of projects."""
    pass

###############################################################################
class Project(Aggregate):
    """A project containing several pools possibly spanning several runs."""

    def __repr__(self): return '<%s object "%s" with %i pools/presamples>' % \
                               (self.__class__.__name__, self.name, len(self))

    def __init__(self, name, pools, dir_path):
        """Please specify the name of the project, the pools it must contain and
        finally the directory path where the results should be generated."""
        # Attributes #
        self.name = name
        self.pools, self.children = pools, pools
        self.loaded = False
        # Dir #
        self.base_dir = dir_path + self.name + '/'
        self.p = AutoPaths(self.base_dir, self.all_paths)
        # Extra #
        self.meta_data_path = illumitag.repos_dir + 'projects/' + self.name + '.csv'

    @property
    def long_name(self): return self.first.project_long_name

    @property
    def title(self):
        pool = [p for p in self.pools if 'project_title' in p.info]
        if pool: return pool[0].info['project_title']

    @property
    def abstract(self):
        pool = [p for p in self.pools if 'project_abstract' in p.info]
        if pool: return pool[0].info['project_abstract']