# Futures #
from __future__ import division

# Built-in modules #
import shutil, inspect
from collections import OrderedDict

# Internal modules #
import sifes
from sifes.report import ReportTemplate

# First party modules #
from pymarktex          import Document
from pymarktex.figures  import ScaledFigure
from plumbing.autopaths import FilePath
from plumbing.common    import split_thousands
from plumbing.cache     import property_cached

# Third party modules #
from tabulate import tabulate

###############################################################################
class MultiplexReport(Document):
    """A full report generated in PDF for every Demultiplexer object.
    Replaces the report for that project."""

    def __repr__(self): return '<%s object on %s>' % (self.__class__.__name__, self.parent)

    def __init__(self, plexer):
        # Attributes #
        self.plexer, self.parent = plexer, plexer
        # The project #
        self.plexproj = plexer.plexproj
        # Automatic paths #
        self.base_dir    = self.plexproj.p.report_dir
        self.output_path = self.plexproj.p.report_pdf
        # Basic export path #
        self.copy_base = sifes.reports_dir + self.plexproj.short_name + '/' + self.parent.short_name + '.pdf'
        self.copy_base = FilePath(self.copy_base)

    @property_cached
    def template(self): return MultiplexTemplate(self)

    def generate(self):
        # Dynamic templates #
        self.markdown = unicode(self.template)
        # Render to latex #
        self.make_body()
        self.make_latex(params={'title': 'Demultiplexing report'})
        self.make_pdf(safe=True, include_src=True)
        # Copy to reports directory #
        self.copy_base.directory.create(safe=True)
        shutil.copy(self.output_path, self.copy_base)
        # Return #
        return self.output_path

###############################################################################
class MultiplexTemplate(ReportTemplate):
    """All the parameters to be rendered in the markdown template"""
    delimiters = (u'{{', u'}}')

    def __repr__(self): return '<%s object on %s>' % (self.__class__.__name__, self.parent)

    def get_associated_template(self):
        return False

    def __init__(self, report):
        # Attributes #
        self.report, self.parent = report, report
        self.plexer    = self.parent.plexer
        self.plexproj  = self.parent.plexproj
        self.plexfiles = self.parent.plexer.plexfiles
        self.samples   = self.parent.plexer.samples

    ############## General information ##############
    def plex_short_name(self): return self.plexproj.short_name
    def plex_long_name(self):  return self.plexproj.long_name

    ############## Input ##############
    def count_input_files(self):  return len(self.plexproj)
    def count_pools(self):        return len(self.plexfiles)

    def input_table(self):
        # The columns #
        info = OrderedDict((
            ('Name',    lambda p: "**" + p.name + "**"),
            ('Files',   lambda p: "`" + str(len(p.inputs)) + "`"),
            ('Samples', lambda p: "`" + str(len(p.samples)) + "`"),
            ('Reads',   lambda p: split_thousands(p.pair.count)),
            ('Loss',    lambda p: "%.1f%%" % (100.0 * (1 - sum(map(len, p.samples)) / p.pair.count))),
        ))
        # The table #
        table = [[i+1] + [f(self.plexfiles[i]) for f in info.values()] for i in range(len(self.plexfiles))]
        # Make it as text #
        table = tabulate(table, headers=['#'] + info.keys(), numalign="right", tablefmt="pipe")
        # Add caption #
        return table + "\n\n   : Summary information for all multiplexing pools."

    ############## Output ##############
    def count_real_samples(self): return len(self.samples)
    def count_loss(self):
        reads_lost = sum(map(len, self.plexproj.samples)) - sum(map(len, self.samples))
        percent = 100.0 * reads_lost / sum(map(len, self.plexproj.samples))
        return "%s (%.1f%%)" % (split_thousands(reads_lost), percent)

    def output_table(self):
        # The columns #
        info = OrderedDict((
            ('Name',  lambda s: "**" + s.short_name + "**"),
            ('Group', lambda s: "`" + s.info['multiplexed_in'] + "`"),
            ('Reads', lambda s: split_thousands(s.pair.count)),
        ))
        # The table #
        table = [[i+1] + [f(self.samples[i]) for f in info.values()] for i in range(len(self.samples))]
        # Make it as text #
        table = tabulate(table, headers=['#'] + info.keys(), numalign="right", tablefmt="pipe")
        # Add caption #
        return table + "\n\n   : Summary information for all final samples."

    ############## Misstags ##############
    def mistags(self):
        if not self.plexer.results: return False
        params = ('not_both_primers', 'unknown_fwd_barcode', 'unknown_rev_barcode',
                  'primer_mismatches', 'barcode_mismatches', 'mistag_heatmap')
        return {p:getattr(self, p) for p in params}

    def not_both_primers(self):     return split_thousands(self.plexer.results.not_both_primers)
    def unknown_fwd_barcode(self):  return split_thousands(self.plexer.results.unknown_fwd_barcode)
    def unknown_rev_barcode(self):  return split_thousands(self.plexer.results.unknown_rev_barcode)

    def primer_mismatches(self):   return self.plexer.first.primer_mismatches
    def barcode_mismatches(self):  return self.plexer.first.barcode_mismatches

    def mistag_heatmap(self):
        caption = "Heatmap representing mistaged read pairs. Existing samples are outlined in black."
        graph = self.plexer.results.graphs.mistag_heatmap(rerun=True)
        return str(ScaledFigure(graph.path, caption, inspect.stack()[0][3]))

    ############## Predictions ##############
    def predictions(self):
        return False
        output = ""
        for i, p in enumerate(self.plexfiles):
            barcodes = p.guess_barcodes().most_common(len(p.samples)+2)
            output += "\n\n### %s:\n\n" % p.name
            for k,v in barcodes: output += "* %s: %s\n" % (k, v)
        return output
