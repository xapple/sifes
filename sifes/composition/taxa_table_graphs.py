# Futures #
from __future__ import division

# Built-in modules #

# Internal modules #
from plumbing.graphs import Graph, cool_colors

# Third party modules #
from matplotlib import pyplot

# Constants #
__all__ = ['CompositionDomain',
           'CompositionPhyla',
           'CompositionClass',
           'CompositionOrder',
           'CompositionFamily',
           'CompositionGenus',
           'CompositionSpecies']

################################################################################
class TaxaBarstack(Graph):
    """Distribution of named taxa by sample (at different ranks)."""
    base_rank  = -1

    short_name    = 'taxa_barstack'
    width         = 16.0
    height        = 10.0
    bottom        = 0.35
    top           = 0.97
    left          = 0.04
    right         = 0.98
    legend_anchor = -0.15
    x_label       = 'Sample short names'
    y_label       = 'Relative abundances in percent'

    def plot(self, **kwargs):
        # Data #
        taxa_table = self.parent.taxa_tables_by_rank[self.base_rank - 1]
        self.frame = taxa_table.apply(lambda x: 100*x/x.sum(), axis=1)
        # Special case where there is only one taxa e.g. only 'Bacteria' #
        if len(self.frame.columns) < 2 : colors = 'gray'
        else:                            colors = cool_colors
        # Plot #
        axes = self.frame.plot(kind='bar', stacked=True, color=colors)
        fig = pyplot.gcf()
        # Other #
        title = 'Taxonomic relative abundances per sample at rank %i (%s).'
        title = title % (self.base_rank, self.parent.taxonomy.results.rank_names[self.base_rank-1])
        axes.set_title(title)
        axes.set_ylim([0,100])
        # Put a legend below current axis
        axes.legend(loc='upper center', bbox_to_anchor=(0.5, self.legend_anchor), fancybox=True, shadow=True, ncol=5)
        # Change font of sample names #
        axes.set_xticklabels(axes.get_xticklabels(), fontname='Menlo', fontweight='bold', size= 'large')
        #for label in (plot.get_xticklabels() + plot.get_yticklabels()):  label.set_fontproperties(font)
        # Save it #
        self.save_plot(fig, axes, **kwargs)
        pyplot.close(fig)

################################################################################
class CompositionDomain(TaxaBarstack):
    base_rank  = 1
    label      = "Domain"
    short_name = 'taxa_barstack_domain'

class CompositionPhyla(TaxaBarstack):
    base_rank  = 2
    label      = "Phyla"
    short_name = 'taxa_barstack_phyla'

class CompositionClass(TaxaBarstack):
    base_rank  = 3
    label      = "Class"
    short_name = 'taxa_barstack_class'

class CompositionOrder(TaxaBarstack):
    base_rank  = 4
    label      = "Order"
    short_name = 'taxa_barstack_order'

class CompositionFamily(TaxaBarstack):
    base_rank  = 5
    label      = "Family"
    short_name = 'taxa_barstack_family'

class CompositionGenus(TaxaBarstack):
    base_rank  = 6
    label      = "Genus"
    short_name = 'taxa_barstack_genus'

class CompositionSpecies(TaxaBarstack):
    base_rank  = 7
    label      = "Species"
    short_name = 'taxa_barstack_species'
