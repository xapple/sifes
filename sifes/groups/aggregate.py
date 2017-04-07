# Futures #
from __future__ import division

# Built-in modules #
import warnings

# Internal modules #

# First party modules #
from plumbing.cache             import property_cached
from plumbing.autopaths         import AutoPaths, DirectoryPath

###############################################################################
class Aggregate(object):
    """An arbitrary aggregate of several samples."""

    all_paths = """
    /report/report.pdf
    /graphs/
    /plexfiles/
    /logs/
    /cluster/
    /results/slurm_report.csv
    """

    def __repr__(self): return '<%s object "%s" with %i samples>' % \
                               (self.__class__.__name__, self.name, len(self))

    def __str__(self):             return self.short_name
    def __iter__(self):            return iter(self.children)
    def __len__(self):             return len(self.children)
    def __contains__(self, item):  return item in self.children

    def __add__(self, other):
        name     = self.name + " and " + other.name
        children = self.children + other.children
        return self.__class__(name, children, sort=False)

    def __getitem__(self, key):
        if   isinstance(key, basestring):  return [c for c in self.children if str(c) == key][0]
        elif isinstance(key, int):
            if key < 0:                    return self.children[key]
            if hasattr(self.first, 'num'): return [c for c in self.children if int(c.num) == key][0]
            else:                          return self.children[key]
        elif isinstance(key, slice):       return self.children[key]
        elif isinstance(key, tuple):       return (self[k] for k in key)
        elif isinstance(key, list):        return [self[k] for k in key]
        else:                              raise TypeError("The key '%s' is not found" % key)

    @property
    def first(self): return self.children[0]

    def __init__(self, name, samples, out_dir, sort=True):
        # Attributes #
        self.name       = name
        self.short_name = name
        self.samples    = samples
        self.children   = samples
        # Check names are unique #
        names = [s.short_name for s in self.samples]
        assert len(names) == len(set(names))
        # Are the samples numbered #
        have_numbers = all(s.info.get('sample_num') for s in samples)
        if not have_numbers: warnings.warn("Not all samples of project '%s' were numbered." % self)
        # Sort the samples #
        if sort:
            if have_numbers: samples.sort(key=lambda s: int(s.info['sample_num']))
            else:            samples.sort(key=lambda s: s.short_name)
        # Base directory #
        self.base_dir = DirectoryPath(out_dir + self.name + '/')

    @property_cached
    def p(self): return AutoPaths(self.base_dir, self.all_paths)

    @property_cached
    def cluster(self):
        from sifes.groups.cluster import Cluster # To avoid a circular import
        return Cluster(self.name, self.samples[:], self.p.cluster_dir)