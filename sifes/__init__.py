b"""This module needs Python 2.7.x"""

# Built-in modules #
import os, sys, glob

# Constants #
url         = 'http://github.com/xapple/sifes/'
__version__ = '2.0.0'
home        = os.environ.get('HOME', '~') + '/'

# No need for an X display #
import matplotlib
matplotlib.use('Agg', warn=False)

# Get paths to module #
self       = sys.modules[__name__]
module_dir = os.path.dirname(self.__file__)

# The module is a maybe git repository #
from plumbing.git import GitRepo
repos_dir = os.path.abspath(module_dir + '/../') + '/'
git_repo  = GitRepo(repos_dir, empty=True)

# Output directories #
view_dir    = home     + 'SIFES/views/'
project_dir = view_dir + 'projects/'
samples_dir = view_dir + 'samples/'
reports_dir = home     + 'SIFES/reports/'

# Internal modules #
from sifes.groups.projects import Projects, Project
from sifes.groups.samples import Sample

###############################################################################
# The main objects, empty at first, call load() to populate them #
samples    = []
_projects  = []
projects   = Projects(_projects)

###############################################################################
def load(json_dir_path):
    """Will load all the JSON files found in the given json_dir_path.
    If all the files are from the same project, return that project."""

    # Load all found JSON files #
    json_paths  = glob.glob(json_dir_path + '*.json')
    new_samples = [Sample(j) for j in json_paths]
    samples.append(new_samples)

    # Compose into projects #
    proj_names = sorted(list(set([s.project_short_name for s in samples])))
    new_projs  = [Project(name, [s for s in samples if s.project_short_name==name]) for name in proj_names]
    _projects.append(new_projs)

    # Link the samples to their project #
    for s in samples: s.project = projects[s.project_short_name]

    # Return project #
    proj_name = set(s.project_short_name for s in samples)
    if len(proj_name) == 1: return projects[proj_name.pop()]