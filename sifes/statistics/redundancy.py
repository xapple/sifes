# Futures #
from __future__ import division

# Internal modules #

# First party modules #
from plumbing.cache     import property_cached
from plumbing.autopaths import AutoPaths

# Third party modules #
from skbio.stats.ordination import rda

# Third party modules #
import pandas

################################################################################
class RedundancyAnalysis(object):
    """See https://sites.google.com/site/mb3gustame/constrained-analyses/rda
    for more information about RDA.
    
    * X is a matrix representing the explanatory variables (temp, pH, depth).
    * Y is a matrix representing the response variables (species abundance)."""

    explanatory_keys = ('custom_grouping', 'custom_attribute', 'replicate_id',
                        'depth',           'temperature',      'salinity', 'grain_size')

    # Paths #
    all_paths = """
    /explanatory.tsv
    /response.tsv
    """

    def __init__(self, cluster, result_dir):
        # Attributes #
        self.cluster, self.parent = cluster, cluster
        self.result_dir  = result_dir
        # Auto paths #
        self.base_dir = self.result_dir
        self.p = AutoPaths(self.base_dir, self.all_paths)

    def run(self):
        # The matrices #
        X = self.explanatory_variables
        Y = self.response_variables
        # Save the X and Y matrix to files for other programs #
        X.to_csv(self.p.explanatory.path, sep='\t', float_format='%.5g', encoding='utf-8')
        Y.to_csv(self.p.response.path,    sep='\t', float_format='%.5g', encoding='utf-8')
        # Use the sci-kit bio implementation #
        ordination_result = rda(Y, X, scale_Y=True)
        # Return #
        return ordination_result

    @property_cached
    def explanatory_variables(self):
        """This is the X matrix."""
        # Pass a list of dicts to the constructor and separately the row names #
        return pandas.DataFrame(({k: s.info[k] for k in self.explanatory_keys} for s in self.cluster),
                                index=(s.short_name for s in self.cluster))

    @property_cached
    def response_variables(self):
        """This is the Y matrix."""
        # Filter the OTU table to keep only the top 500 #
        N   = 500
        df  = self.cluster.otu_table.results.otu_table
        top = df.sum().sort_values(ascending=False).index[0:N]
        return df[top]