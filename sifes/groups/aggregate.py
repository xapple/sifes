# Futures #
from __future__ import division

# Built-in modules #

# Internal modules #

# First party modules #
from plumbing.cache             import property_cached
from plumbing.autopaths         import AutoPaths

###############################################################################
class Collection(object):
    """A collection of aggregates."""

    def __repr__(self): return 'Collection: %s' % self.children
    def __iter__(self): return iter(self.children)
    def __len__(self): return len(self.children)
    def __add__(self, other): return self.__class__(self.children + other.children)

    def __init__(self, children):
        self.children = children

    @property
    def first(self): return self.children[0]

    def __getitem__(self, key):
        if   isinstance(key, basestring):  return [c for c in self.children if str(c) == key][0]
        elif isinstance(key, int):
            if key < 0:                    return self.children[key]
            if hasattr(self.first, 'num'): return [c for c in self.children if int(c.num) == key][0]
            else:                          return self.children[key]
        elif isinstance(key, slice):       return self.children[key]
        else:                              raise TypeError('key')

###############################################################################
class Aggregate(object):
    """An arbitrary aggregate of several samples."""

    all_paths = """
    /graphs/
    /logs/
    /cluster/
    /results/slurm_report.csv
    """

    def __repr__(self): return '<%s object "%s" with %i samples>' % \
                               (self.__class__.__name__, self.name, len(self))

    def __str__(self):  return self.short_name
    def __iter__(self): return iter(self.children)
    def __len__(self):  return len(self.children)

    def __getitem__(self, key):
        if   isinstance(key, basestring):  return [c for c in self.children if str(c) == key][0]
        elif isinstance(key, int):
            if hasattr(self.first, 'num'): return [c for c in self.children if int(c.num) == key][0]
            else:                          return self.children[key]
        elif isinstance(key, slice):       return self.children[key]
        else:                              raise TypeError('key')

    @property
    def first(self): return self.children[0]

    def __init__(self, name, samples, out_dir):
        # Attributes #
        self.name       = name
        self.short_name = name
        self.samples    = samples
        self.children   = samples
        # Check names are unique #
        names = [s.short_name for s in self.samples]
        assert len(names) == len(set(names))
        # Sort the samples #
        if all(s.info.get('num') for s in samples): samples.sort(key=lambda s: s.info['num'])
        else:                                       samples.sort(key=lambda s: s.short_name)
        # Dir #
        self.base_dir = out_dir + self.name + '/'
        self.p = AutoPaths(self.base_dir, self.all_paths)

    @property_cached
    def cluster(self):
        from sifes.groups.cluster import Cluster # To avoid circular import
        return Cluster(self.samples, self.name, self.p.cluster_dir)