# Built-in modules #

# Internal modules #
import sifes
from sifes.groups.aggregate import Collection, Aggregate

# First party modules #

# Third party modules #

###############################################################################
class Projects(Collection):
    """A collection of projects."""
    pass

###############################################################################
class Project(Aggregate):
    """A project containing several samples."""

    def __init__(self, name, samples):
        """Please specify the name of the project and the samples it must contain."""
        # Base directory #
        out_dir = sifes.project_dir + samples[0].info.get('organization', '') + '/'
        # Super #
        super(self.__class__,self).__init__(name, samples, out_dir)

    @property
    def long_name(self): return self.first.project_long_name