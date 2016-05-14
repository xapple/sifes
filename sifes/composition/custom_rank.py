# Built-in modules #
from collections import defaultdict

# Internal modules #
from illumitag.clustering.composition import Composition
from plumbing.cache import property_cached
from illumitag.clustering.composition.plots import TaxaBarstack

# Third party modules #
import pandas

###############################################################################
class CustomRank(Composition):
    """Depending on the variable `base_rank`, this custom taxa table contains
    different relative level of a group of organisms in the taxonomic hierarchy.
    For instance if `base_rank` is 2 it will only include the phyla level,
    and sometimes the class levels, when phyla are over a certain threshold.
    Additionally, low abundance taxa are grouped into 'Others'."""

    @property_cached
    def taxa_table(self):
        # Build a new frame #
        result = defaultdict(lambda: defaultdict(int))
        for sample_name, column in self.parent.otu_table.iterrows():
            for otu_name, count in column.iteritems():
                taxa = self.parent.assignments[otu_name]
                rank = taxa[self.base_rank] if len(taxa) > self.base_rank else taxa[-1]
                result[rank][sample_name] += count
        # Fill the holes #
        result = pandas.DataFrame(result)
        result = result.fillna(0)
        result = result.astype(int)
        # Ungroup the high abundance ranks #
        high_abundance  = result.sum() > self.high_abundance_threshold
        self.high_ranks = [r for r, true_false in high_abundance.iteritems() if true_false]
        for rank in self.high_ranks:
            # Find OTUs, but don't consider the ones that were dropped because of zero counts #
            cond = lambda taxa: len(taxa) > self.base_rank and taxa[self.base_rank] == rank
            otus = [name for name in self.parent.otu_table.columns if cond(self.parent.assignments[name])]
            # Make new columns #
            new_columns = defaultdict(lambda: defaultdict(int))
            for otu_name in otus:
                taxa = self.parent.assignments[otu_name]
                lower_rank = taxa[self.base_rank + 1] if len(taxa) > self.base_rank + 1 else taxa[-1]
                lower_rank += '*'
                for sample_name, column in self.parent.otu_table.iterrows():
                    count = column[otu_name]
                    if count == 0: continue
                    new_columns[lower_rank][sample_name] += count
            new_columns = pandas.DataFrame(new_columns)
            # Check #
            assert (new_columns.sum(axis=1) == result[rank]).all()
            # Switch them in place #
            result = result.drop(rank, axis=1)
            result = result.join(new_columns)
        # Group low abundant into 'others' #
        low_abundance = result.sum() < self.low_abundance_threshold
        self.low_ranks = [r for r, true_false in low_abundance.iteritems() if true_false]
        the_others_count = result.loc[:, low_abundance].sum(axis=1)
        result = result.loc[:, ~low_abundance]
        result['Others'] = the_others_count
        # Sort the table by sum #
        sums = result.sum()
        sums = sums.sort_values(ascending=False)
        result = result.reindex_axis(sums.keys(), axis=1)
        # Return result #
        return result

    @property_cached
    def graphs(self):
        return [TaxaBarstack(self)]

################################################################################
class CompositionPhyla(CustomRank):
    base_rank = 2
    high_abundance_threshold = 300000
    low_abundance_threshold  = 3000

class CompositionClass(CustomRank):
    base_rank = 3
    high_abundance_threshold = 100000
    low_abundance_threshold  = 1000

class CompositionOrder(CustomRank):
    base_rank = 4
    high_abundance_threshold = 30000
    low_abundance_threshold  = 300

class CompositionFamily(CustomRank):
    base_rank = 5
    high_abundance_threshold = 1000
    low_abundance_threshold  = 50

class CompositionGenus(CustomRank):
    base_rank = 6
    high_abundance_threshold = 500
    low_abundance_threshold  = 10
