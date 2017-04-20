# Futures #
from __future__ import division

# Built-in modules #

# Internal modules #
from sifes.composition.taxa_table import TaxaTable

# First party modules #
from plumbing.autopaths import AutoPaths
from plumbing.cache     import property_cached

###############################################################################
class SubTaxaTable(TaxaTable):
    """Takes only a specific phylum or clade and makes focused taxa tables and
    associated graphs, throwing away all other phyla or clades."""

    def __init__(self, taxa_table, filter_rank, filter_name, result_dir):
        # Attributes #
        self.taxa_table  = taxa_table
        self.filter_rank = filter_rank
        self.filter_name = filter_name
        self.result_dir  = result_dir
        # Short cuts #
        self.rank_names = self.taxa_table.rank_names
        # Auto paths #
        self.base_dir = self.result_dir + self.filter_name.lower().replace(' ', '_') + '/'
        self.p        = AutoPaths(self.base_dir, self.all_paths)

    @property_cached
    def otu_df(self):
        """We just need to fake a new self.otu_df property and we are golden."""
        # Main objects #
        otu_df       = self.taxa_table.otu_df.copy()
        assignments  = self.taxa_table.taxonomy.results.assignments
        # Remove unwanted #
        for otu_name in otu_df:
            species = assignments[otu_name]
            if not self.filter_name in species[self.filter_rank]:
                otu_df = otu_df.drop(otu_name, 1)
        # Return #
        return otu_df
