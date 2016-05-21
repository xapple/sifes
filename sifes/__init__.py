b"""This module needs Python 2.7.x"""

# Built-in modules #
import os, sys

# Constants #
url         = 'http://xapple.github.io/sifes/'
repo_url    = 'http://github.com/xapple/sifes/'
__version__ = '2.0.0'
home        = os.environ.get('HOME', '~') + '/'

# No need for an X display #
import matplotlib
matplotlib.use('Agg', warn=False)

# Get paths to module #
self       = sys.modules[__name__]
module_dir = os.path.dirname(self.__file__)

# The module is maybe a git repository #
from plumbing.git import GitRepo
repos_dir = os.path.abspath(module_dir + '/../') + '/'
git_repo  = GitRepo(repos_dir, empty=True)

# Output directories #
view_dir     = home     + 'SIFES/views/'
project_dir  = view_dir + 'projects/'
samples_dir  = view_dir + 'samples/'
clusters_dir = view_dir + 'clusters/'
reports_dir  = home     + 'SIFES/reports/'
bundles_dir  = home     + 'SIFES/bundles/'

# Internal modules #
from sifes.groups.projects import Projects

# The main objects, empty at first, call load() to populate them #
samples    = []
_projects  = []
projects   = Projects(_projects)

# Expose functions #
from sifes.load import load