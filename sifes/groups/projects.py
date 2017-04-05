# Built-in modules #

# Internal modules #
from sifes.groups.aggregate import Aggregate
from sifes.groups.lump      import Lump
import sifes

# First party modules #
from plumbing.processes import prll_map
from plumbing.timer     import Timer

# Third party modules #

###############################################################################
class Projects(Lump):
    """A collection of projects."""
    def __init__(self, name="all_projects", *args, **kwargs):
        super(self.__class__,self).__init__(name, *args, **kwargs)

###############################################################################
class Project(Aggregate):
    """A project containing several samples."""

    def __init__(self, name, samples, *args, **kwargs):
        """Please specify the name of the project and the samples it must contain."""
        # Base directory #
        out_dir = sifes.project_dir + samples[0].info.get('organization', '') + '/'
        # Super #
        super(self.__class__,self).__init__(name, samples, out_dir, *args, **kwargs)

    @property
    def long_name(self): return self.first.project_long_name

    @property
    def organization(self): return self.first.info.get('organization', '')

    def run_analysis(self):
        print "Join reads"
        with Timer(): prll_map(lambda s: s.joiner.run(cpus=1), self)
        for s in self: print s.short_name, s.joiner.results.unassembled_percent
