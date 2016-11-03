# Futures #

# Built-in modules #

# Internal modules #
import sifes

# First party modules #
from plumbing.autopaths import AutoPaths
from plumbing.cache     import property_cached

###############################################################################
class Lump(object):
    """Put several Aggregates objects (or Projects objects) together."""

    all_paths = """
    """

    def __repr__(self): return '<%s object "%s" with %i children>' % \
                               (self.__class__.__name__, self.name, len(self))
    def __iter__(self):            return iter(self.children)
    def __len__(self):             return len(self.children)
    def __contains__(self, item):  return item in self.children

    def __add__(self, other):
        return self.__class__(self.name + " and " + other.name, self.children + other.children)

    def __getitem__(self, key):
        if   isinstance(key, basestring):  return [c for c in self.children if str(c) == key][0]
        elif isinstance(key, int):
            if key < 0:                    return self.children[key]
            if hasattr(self.first, 'num'): return [c for c in self.children if int(c.num) == key][0]
            else:                          return self.children[key]
        elif isinstance(key, slice):       return self.children[key]
        else:                              raise TypeError('key')

    def __init__(self, name, aggregates, base_dir=sifes.lumps_dir):
        # Attributes #
        self.name = name
        # The list #
        self.aggregates, self.children = aggregates, aggregates
        # Sort them #
        if all(hasattr(c, "num") for c in self.children): self.children.sort(key = lambda x: x.num)
        # Base directory #
        self.base_dir = base_dir

    @property
    def first(self): return self.children[0]

    @property
    def samples(self):
        return [s for a in self.aggregates for s in a]

    @property_cached
    def p(self): return AutoPaths(self.base_dir, self.all_paths)