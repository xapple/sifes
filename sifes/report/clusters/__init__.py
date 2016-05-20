# Futures #
from __future__ import division

# Built-in modules #
import shutil
from collections import OrderedDict

# Internal modules #
import sifes
from sifes.report import ReportTemplate

# First party modules #
from pymarktex         import Document
from plumbing.common   import split_thousands, andify
from pymarktex.figures import ScaledFigure

# Third party modules #
from tabulate import tabulate

###############################################################################
class ClusterReport(Document):
    """A full report generated in PDF for every Cluster object."""

    def __repr__(self): return '<%s object on %s>' % (self.__class__.__name__, self.parent)

    def __init__(self, cluster):
        # Attributes #
        self.cluster, self.parent = cluster, cluster
        # Automatic paths #
        self.base_dir    = self.cluster.p.report_dir
        self.output_path = self.cluster.p.report_pdf

    def generate(self):
        # Dynamic templates #
        self.markdown = unicode(ClusterTemplate(self))
        # Render to latex #
        self.make_body()
        self.make_latex()
        self.make_pdf(safe=True)
        # Copy to reports directory #
        shutil.copy(self.output_path, self.copy_base)
        # Return #
        return self.output_path

    copy_base = property(lambda self: sifes.reports_dir + self.cluster.project.name + '/' + self.cluster.name + '.pdf')

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
        if not self.cluster.count_dropped_samples: return False
        params = ('count_dropped_samples', 'read_count_cutoff', 'read_count_cutoff_percentile')
        return {p:getattr(self, p) for p in params}
    def count_dropped_samples(self):        return self.cluster.count_dropped_samples
    def read_count_cutoff(self):            return self.cluster.read_count_cutoff
    def read_count_cutoff_percentile(self): return self.cluster.read_count_cutoff_percentile

    # Samples #
    def sample_table(self):
        # The columns #
        info = OrderedDict((
            ('Name',        lambda s: "**" + s.short_name + "**"),
            ('Reference',   lambda s: "`" + s.name + "`"),
            ('Description', lambda s: s.long_name),
            ('Reads lost',  lambda s: "%.1f%%" % (100 - ((len(s.filter.results.clean) / len(s))*100))),
            ('Reads left',  lambda s: split_thousands(len(s.filter.results.clean))),
        ))
        # The table #
        table = [[i+1] + [f(self.cluster.samples[i]) for f in info.values()] for i in range(len(self.cluster))]
        # Make it as text #
        table = tabulate(table, headers=['#'] + info.keys(), numalign="right", tablefmt="pipe")
        # Add caption #
        return table + "\n\n   : Summary information for all samples."

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
    def clustering_threshold(self):   return "%.1f%%" % self.centering.threshold
    def otus_total(self):             return split_thousands(len(self.centering.results.centers))

    # Classification #
    def classify_citation(self):    return "the %s method (%s)" % (self.taxonomy.long_name, self.taxonomy.version)
    def classify_database(self):    return self.taxonomy.database
    def otu_classified_table(self):
        info = OrderedDict((('Rank',         lambda i: "**" + self.taxonomy.results.rank_names[i] + "**"),
                            ('Classified',   lambda i:        self.taxonomy.results.count_assigned[i]),
                            ('Unclassified', lambda i:        self.taxonomy.results.count_unassigned[i])))
        table = [[i+1] + [f(i) for f in info.values()] for i in range(len(self.taxonomy.results.rank_names))]
        table = tabulate(table, headers=['#'] + info.keys(), numalign="right", tablefmt="pipe")
        return table + "\n\n   : Summary information for all samples."

    # OTU table filtering #
    def unwanted_phyla(self):       return andify(self.otu_table.unwanted_phyla)
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
    def phyla_barstack(self):
        caption = "Relative abundances per sample on the phyla level"
        path    = self.taxa_table.results.graphs.taxa_barstack_phyla
        label   = "phyla_barstack"
        return str(ScaledFigure(path, caption, label))
    def class_barstack(self):
        caption = "Relative abundances per sample on the class level"
        path    = self.taxa_table.results.graphs.taxa_barstack_class
        label   = "class_barstack"
        return str(ScaledFigure(path, caption, label))
    def order_barstack(self):
        caption = "Relative abundances per sample on the order level"
        path    = self.taxa_table.results.graphs.taxa_barstack_order
        label   = "order_barstack"
        return str(ScaledFigure(path, caption, label))

    # Comparison #
    def otu_nmds(self):
        caption = "NMDS using the OTU table for %i samples" % len(self.cluster)
        path    = self.cluster.nmds_graph()
        label   = "otu_nmds"
        return str(ScaledFigure(path, caption, label))

    # Diversity #
    pass

    # Beta-dispersion #
    pass
