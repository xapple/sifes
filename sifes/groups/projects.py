# Built-in modules #

# Internal modules #
import sifes
from sifes.groups.aggregate import Aggregate
from sifes.groups.lump      import Lump

# First party modules #

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