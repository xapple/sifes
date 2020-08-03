# Futures #
from __future__ import division

# Built-in modules #
import inspect, shutil
from collections import OrderedDict

# Internal modules #
import sifes
from sifes.report import ReportTemplate
from sifes.metadata.correspondence import reverse_corr

# First party modules #
from plumbing.autopaths import FilePath
from plumbing.common    import split_thousands, andify
from plumbing.cache     import property_cached
from pymarktex.figures  import ScaledFigure
from pymarktex          import Document

# Third party modules #
from tabulate import tabulate

###############################################################################
class ClusterReport(Document):
    """A full report generated in PDF for every Cluster object."""

    # Parameters #
    default_taxa_graph_levels = (2, 3, 4)

    def __repr__(self): return '<%s object on %s>' % (self.__class__.__name__, self.parent)

    def __init__(self, cluster):
        # Attributes #
        self.cluster, self.parent = cluster, cluster
        # Automatic paths #
        self.base_dir    = self.cluster.p.report_dir
        self.output_path = self.cluster.p.report_pdf
        # Basic export path #
        self.copy_base = sifes.reports_dir + self.cluster.project.name + '/' + self.cluster.name + '.pdf'
        self.copy_base = FilePath(self.copy_base)

    @property_cached
    def template(self): return ClusterTemplate(self)

    def generate(self):
        # Message #
        print("Making report for cluster '%s'" % self.cluster.name)
        # Dynamic templates #
        self.markdown = unicode(self.template)
        # Render to latex #
        self.make_body()
        self.make_latex({'title': 'Cluster report'})
        self.make_pdf(safe=True)
        # Copy to reports directory #
        self.copy_base.directory.create(safe=True)
        shutil.copy(self.output_path, self.copy_base)
        # Return #
        return self.output_path

###############################################################################
class ClusterTemplate(ReportTemplate):
    """All the parameters to be rendered in the markdown template"""
    delimiters = (u'{{', u'}}')

    def __repr__(self): return '<%s object on %s>' % (self.__class__.__name__, self.parent)

    def __init__(self, report):
        # Attributes #
        self.report, self.parent = report, report
        self.cluster = self.parent.cluster
        self.project = self.cluster.project
        self.samples = self.cluster.samples
        # Convenience #
        self.centering  = self.cluster.centering
        self.taxonomy   = self.cluster.taxonomy
        self.otu_table  = self.cluster.otu_table
        self.taxa_table = self.cluster.taxa_table

    # General information #
    def cluster_name(self):  return self.cluster.name
    def count_samples(self): return len(self.cluster)
    def project_sentence(self):
        if self.cluster.project is None: return ""
        msg = "It corresponds to project code '%s' ('%s')"
        return msg % (self.cluster.project.name, self.cluster.project.long_name)

    # Dropped samples #
    def dropped_samples(self):
        if not self.cluster.bad_samples: return False
        params = ('count_dropped_samples', 'read_count_cutoff', 'read_count_cutoff_percentile')
        return {p:getattr(self, p) for p in params}
    def count_dropped_samples(self):        return len(self.cluster.bad_samples)
    def read_count_cutoff(self):            return self.cluster.read_count_cutoff
    def read_count_cutoff_factor(self):     return self.cluster.read_count_cutoff_factor

    # Samples #
    def sample_table(self):
        # The columns #
        info = OrderedDict((
            ('Name',        lambda s: "**" + s.short_name + "**"),
            ('Description', lambda s: s.long_name),
            ('Reads lost',  lambda s: "%.1f%%" % (100 - ((len(s.filter.results.clean) / len(s))*100))),
            ('Reads left',  lambda s: split_thousands(len(s.filter.results.clean))),
        ))
        # The table #
        table = [[i+1] + [f(s) for f in info.values()] for i,s in enumerate(self.samples)]
        # Make it as text #
        table = tabulate(table, headers=['#'] + info.keys(), numalign="right", tablefmt="pipe")
        # Add caption #
        return table + "\n\n   : Summary information for all samples."

    # Location #
    def location(self):
        # Check GPS data #
        if not all(s.latitude  for s in self.samples): return False
        if not all(s.longitude for s in self.samples): return False
        # There are as many maps as there are custom groupings #
        string = ""
        for loc_map in self.cluster.locations_maps:
            caption = "Map of samples in the group '%s'" % loc_map.group_name
            path    = loc_map()
            label   = "location_" + loc_map.group_name.lower()
            string += str(ScaledFigure(path, caption, label)) + '\n\n'
        # Return the long string with all the figures #
        return {'location': string}

    # Replicates #
    def replicates(self):
        """A user inputted optional figure."""
        # Check it exists #
        if not self.cluster.p.replicates: return False
        # Return #
        caption = "Sampling strategy and replicates."
        return {'replicates': str(ScaledFigure(self.cluster.p.replicates, caption, "replicates"))}

    # Input data #
    def count_sequences(self): return split_thousands(len(self.cluster.reads))
    def input_length_dist(self):
        caption = "Distribution of sequence lengths at input"
        path    = self.cluster.reads.graphs.length_dist()
        label   = "input_length_dist"
        return str(ScaledFigure(path, caption, label))

    # Clustering #
    def clustering_citation(self):    return "the %s method (%s)" % (self.centering.long_name, self.centering.version)
    def clustering_publication(self): return self.centering.article
    def clustering_threshold(self):   return "%.1f" % self.centering.threshold
    def otus_total(self):             return split_thousands(len(self.centering.results.centers))

    # Classification #
    def classify_citation(self):    return "the '%s' method" % self.taxonomy.long_name
    def classify_database(self):    return "'" + self.taxonomy.database.long_name + "'"
    def otu_classified_table(self):
        info = OrderedDict((('Rank',         lambda i: "**" + self.taxonomy.results.rank_names[i] + "**"),
                            ('Classified',   lambda i:        self.taxonomy.results.count_assigned[i]),
                            ('Unclassified', lambda i:        self.taxonomy.results.count_unassigned[i])))
        table = [[i+1] + [f(i) for f in info.values()] for i in range(len(self.taxonomy.results.rank_names))]
        table = tabulate(table, headers=['#'] + info.keys(), numalign="right", tablefmt="pipe")
        return table + "\n\n   : Summary information for all samples."

    # OTU table filtering #
    def unwanted_taxa(self):        return andify(self.otu_table.unwanted_taxa)
    def otus_filtered(self):        return split_thousands(len(self.otu_table.results.centers))
    def otu_sizes_graph(self):
        caption = "Distribution of OTU sizes"
        path    = self.otu_table.results.graphs.otu_sizes_dist()
        label   = "otu_sizes_graph"
        return str(ScaledFigure(path, caption, label))

    # OTU table graphs #
    def otu_sums_graph(self):
        caption = "Distribution of OTU presence per OTU"
        path    = self.otu_table.results.graphs.otu_sums_graph()
        label   = "otu_sums_graph"
        return str(ScaledFigure(path, caption, label))
    def sample_sums_graph(self):
        caption = "Distribution of OTU presence per sample"
        path    = self.otu_table.results.graphs.sample_sums_graph()
        label   = "sample_sums_graph"
        return str(ScaledFigure(path, caption, label))
    def cumulative_presence(self):
        caption = "Cumulative number of reads by OTU presence"
        path    = self.otu_table.results.graphs.cumulative_presence()
        label   = "cumulative_presence"
        return str(ScaledFigure(path, caption, label))

    # Composition #
    def level_one_barstack(self):
        rank    = self.parent.default_taxa_graph_levels[0]
        graph   = [g for g in self.taxa_table.results.graphs.by_rank if g.base_rank == rank][0]
        caption = "Relative abundances per sample on the '%s' level" % graph.label
        label   = "level_one_barstack"
        path    = graph()
        return str(ScaledFigure(path, caption, label))
    def level_two_barstack(self):
        rank    = self.parent.default_taxa_graph_levels[1]
        graph   = [g for g in self.taxa_table.results.graphs.by_rank if g.base_rank == rank][0]
        caption = "Relative abundances per sample on the '%s' level" % graph.label
        label   = "level_two_barstack"
        path    = graph()
        return str(ScaledFigure(path, caption, label))
    def level_three_barstack(self):
        rank    = self.parent.default_taxa_graph_levels[2]
        graph   = [g for g in self.taxa_table.results.graphs.by_rank if g.base_rank == rank][0]
        caption = "Relative abundances per sample on the '%s' level" % graph.label
        label   = "level_three_barstack"
        path    = graph()
        return str(ScaledFigure(path, caption, label))

    # Comparison #
    def comparison(self):
        if len(self.cluster) < 2: return False
        else: return {'otu_nmds': self.otu_nmds()}
    def otu_nmds(self):
        caption = "NMDS using the OTU table for %i samples" % len(self.cluster)
        path    = self.cluster.nmds_graph()
        label   = "otu_nmds"
        return str(ScaledFigure(path, caption, label))

    # Alpha diversity #
    def alpha_diversity(self):
        return {'alpha_diversity_table': self.alpha_diversity_table(),
                'down_sampled_to':       self.down_sampled_to()}
    def down_sampled_to(self): return split_thousands(self.cluster.down_to)
    def alpha_diversity_table(self):
        from skbio.diversity import alpha_diversity  as alphadiv
        from skbio.stats     import subsample_counts as subsample
        k    = self.cluster.down_to
        info = OrderedDict((
            ('Name',     lambda s: "**" + s.short_name + "**"),
            ('Chao1',    lambda s: alphadiv('chao1',   subsample(s.otu_counts, k))),
            ('Ace',      lambda s: alphadiv('ace',     subsample(s.otu_counts, k))),
            ('Shannon',  lambda s: alphadiv('shannon', subsample(s.otu_counts, k))),
            ('Simpson',  lambda s: alphadiv('simpson', subsample(s.otu_counts, k)))))
        table = [[i+1] + [f(s) for f in info.values()] for i,s in enumerate(self.samples)]
        table = tabulate(table, headers=['#'] + info.keys(), numalign="right", tablefmt="pipe")
        return table + "\n\n   : Summary of diversity estimates for all samples."

    # Seqenv #
    def seqenv(self):
        if not self.cluster.seqenv: return False
        params = ('seqenv_citation', 'seqenv_database', 'seqenv_publication',
                  'seqenv_heatmap')
        return {p:getattr(self, p) for p in params}

    def seqenv_citation(self):    return self.cluster.seqenv.long_name
    def seqenv_database(self):    return self.cluster.seqenv.database
    def seqenv_publication(self): return self.cluster.seqenv.article
    def seqenv_heatmap(self):
        caption = "Heatmap of most significant 'isolation source' controlled vocabulary"
        path    = self.cluster.seqenv.results.graphs.seqenv_heatmap()
        label   = "seqenv_heatmap"
        return str(ScaledFigure(path, caption, label))

    # Subtaxa #
    def sub_taxa(self):
        """Optional extra figures."""
        # Check it exists #
        if not self.cluster.sub_taxa: return False
        # There can be many figures #
        string = ""
        for sub in self.cluster.sub_taxa_tables:
            rank    = sub.filter_rank + 2
            graph   = [g for g in sub.results.graphs.by_rank if g.base_rank == rank][0]
            caption = "Subtaxonomic filter on '%s' (rank %i) displayed on the '%s' level (rank %i)"
            caption = caption % (sub.filter_name, sub.filter_rank, graph.label, rank)
            path    = graph()
            label   = "subtaxa_" + sub.filter_name.lower().replace(' ', '_')
            string += str(ScaledFigure(path, caption, label)) + '\n\n'
        # Return the long string with all the figures #
        return {'sub_taxa': string}

    # Beta-dispersion #
    pass

    # Diversity regression #
    def diversity_reg(self):
        if not self.cluster.otu_table: return False
        params = ('diversity_reg_metric', 'diversity_reg_down_to',
                  'diversity_reg_chao1', 'diversity_reg_ace',
                  'diversity_reg_shannon', 'diversity_reg_simpson',)
        return {p:getattr(self, p) for p in params}

    def diversity_reg_metric(self):
        metric = self.cluster.graphs.diversity_reg_chao1.custom_metadata
        return reverse_corr.get(metric, metric)
    def diversity_reg_down_to(self): return split_thousands(self.cluster.down_to)

    def diversity_reg_chao1(self):
        caption = "Chao1 diversity sample comparison"
        path = self.cluster.graphs.diversity_reg_chao1()
        return str(ScaledFigure(path, caption, inspect.stack()[0][3]))
    def diversity_reg_ace(self):
        caption = "Ace diversity sample comparison"
        path = self.cluster.graphs.diversity_reg_ace()
        return str(ScaledFigure(path, caption, inspect.stack()[0][3]))
    def diversity_reg_shannon(self):
        caption = "Shannon diversity sample comparison"
        path = self.cluster.graphs.diversity_reg_shannon()
        return str(ScaledFigure(path, caption, inspect.stack()[0][3]))
    def diversity_reg_simpson(self):
        caption = "Simpson diversity sample comparison"
        path = self.cluster.graphs.diversity_reg_simpson()
        return str(ScaledFigure(path, caption, inspect.stack()[0][3]))

    # Redundancy Analysis #
    def redundancy_analysis(self):
        params = ('rda_citation'
                  'rda_heatmap')
        return {p:getattr(self, p) for p in params}

    def rda_citation(self):    return self.cluster.seqenv.long_name
    def rda_heatmap(self):
        caption = "Heatmap of most significant 'isolation source' controlled vocabulary"
        path    = self.cluster.seqenv.results.graphs.seqenv_heatmap()
        label   = "seqenv_heatmap"
        return str(ScaledFigure(path, caption, label))
