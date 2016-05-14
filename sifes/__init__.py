b"""This module needs Python 2.7.x"""

# Special variables #
__version__ = '2.0.0'

# Built-in modules #
import os, sys, glob

# Get paths to module #
self = sys.modules[__name__]
module_dir = os.path.dirname(self.__file__)

# The module is a git repository #
from plumbing.git import GitRepo
repos_dir = os.path.abspath(module_dir + '/../') + '/'
git_repo  = GitRepo(repos_dir, empty=True)

# No need for an X display #
import matplotlib
matplotlib.use('Agg', warn=False)

# Internal modules #
from sifes.groups.aggregate import Aggregate
from sifes.groups.projects import Projects, Project
from sifes.groups.samples import Sample

# Constants #
url  = 'http://github.com/xapple/illumitag/'
home = os.environ.get('HOME', '~') + '/'

###############################################################################
# Output directories #
view_dir    = home + 'SIFES/views/'
project_dir = view_dir + 'projects/'
reports_dir = home + 'SIFES/reports/'

# Load all standard pools #
pools_dir = repos_dir + 'json/pools/*/'
json_paths = glob.glob(pools_dir + '*.json')
pools = [Pool(j, view_dir + 'pools/') for j in json_paths]
pools.sort(key=lambda x: str(x))

# Load all presamples #
presamples_dir = repos_dir + 'json/presamples/*/'
json_paths = glob.glob(presamples_dir + '*.json')
presamples = [Presample(j, view_dir + 'presamples/') for j in json_paths]
presamples.sort(key=lambda x: str(x))

# Load all legacy pyrosamples #
pyrosamples_dir = repos_dir + 'json/pyrosamples/'
json_paths = glob.glob(pyrosamples_dir + '*.json')
pyrosamples = [Pyrosample(j, view_dir + 'pyrosamples/') for j in json_paths]
pyrosamples.sort(key=lambda x: str(x))
demultiplexer = Demultiplexer454(pyrosamples)

# All of them together #
all_objects = pools + presamples + pyrosamples

# Compose into runs #
run_nums = sorted(list(set([p.run_num for p in all_objects]))) # [1,2,3,4,5]
runs = [Run(num, [p for p in all_objects if p.run_num==num], view_dir + 'runs/') for num in run_nums]
runs = Runs(runs)
for p in all_objects: p.run = runs[p.run_num]

# Compose into projects #
proj_names = sorted(list(set([p.project_short_name for p in all_objects])))
projects = [Project(name, [p for p in all_objects if p.project_short_name==name], project_dir) for name in proj_names]
projects = Projects(projects)
for p in all_objects: p.project = projects[p.project_short_name]

# Make an aggregate with all pools #
aggregate = Aggregate('all', pools, view_dir + 'aggregates/')

# Make our favorite clusters #
__import__("sifes.clustering.favorites")