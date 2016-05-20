# Built-in modules #
from collections import defaultdict

# Internal modules #
from plumbing.cache import property_cached

# Third party modules #
import pandas

###############################################################################
class TaxonomicRank(object):

    @property_cached
    def taxa_table(self):
        """Uses the otu_table and the taxonomy.assignments to make a new table."""
        # Build a new frame #
        result = defaultdict(lambda: defaultdict(int))
        for sample_name, column in self.parent.otu_table.iterrows():
            for otu_name, count in column.iteritems():
                taxa = self.parent.assignments[otu_name]
                result[taxa[-1]][sample_name] += count
        # Fill the holes #
        result = pandas.DataFrame(result)
        result = result.fillna(0)
        result = result.astype(int)
        # Sort the table by sum #
        sums = result.sum()
        sums = sums.sort_values(ascending=False)
        result = result.reindex_axis(sums.keys(), axis=1)
        # Return #
        return result

################################################################################
class CompositionPhyla(TaxonomicRank):
    base_rank = 2

class CompositionClass(TaxonomicRank):
    base_rank = 3

class CompositionOrder(TaxonomicRank):
    base_rank = 4

class CompositionFamily(TaxonomicRank):
    base_rank = 5

class CompositionGenus(TaxonomicRank):
    base_rank = 6

class CompositionSpecies(TaxonomicRank):

    @property_cached
    def graphs(self):
        return [TaxaBarstackTips(self)]

    def count_otus(self, speices):
        """How many OTUs got this classification"""
        return len([1 for s in self.taxonomy.assignments.values() if s[-1] == speices])
