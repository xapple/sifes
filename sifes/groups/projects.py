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

    def __repr__(self): return '<%s object "%s" with %i samples>' % \
                               (self.__class__.__name__, self.name, len(self))

    def __init__(self, name, samples):
        """Please specify the name of the project and the samples it must contain."""
        out_dir = sifes.project_dir + self.name + '/'
        super(self.__class__,self).__init__(name, samples, out_dir)

    @property
    def long_name(self): return self.first.project_long_name