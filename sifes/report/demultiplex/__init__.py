# Futures #
from __future__ import division

# Built-in modules #
import shutil
from collections import OrderedDict

# Internal modules #
import sifes
from sifes.report import ReportTemplate

# First party modules #
from pymarktex          import Document
from plumbing.autopaths import FilePath
from plumbing.common    import split_thousands

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
        self.plexed = plexer.plexed
        # Automatic paths #
        self.base_dir    = self.plexed.p.report_dir
        self.output_path = self.plexed.p.report_pdf
        # Basic export path #
        self.copy_base = sifes.reports_dir + self.plexed.short_name + '/' + self.parent.short_name + '.pdf'
        self.copy_base = FilePath(self.copy_base)

    def generate(self):
        # Dynamic templates #
        self.main = MultiplexTemplate(self)
        self.markdown = unicode(self.main)
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
        self.plexer  = self.parent.plexer
        self.plexed  = self.parent.plexed
        self.pools   = self.parent.plexer.pools
        self.samples = self.parent.plexer.samples

    ############## General information ##############
    def plex_short_name(self): return self.plexed.short_name
    def plex_long_name(self):  return self.plexed.long_name

    ############## Input ##############
    def count_input_files(self):  return len(self.plexed)
    def count_pools(self):        return len(self.pools)

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
        table = [[i+1] + [f(self.pools[i]) for f in info.values()] for i in range(len(self.pools))]
        # Make it as text #
        table = tabulate(table, headers=['#'] + info.keys(), numalign="right", tablefmt="pipe")
        # Add caption #
        return table + "\n\n   : Summary information for all multiplexing pools."

    ############## Output ##############
    def count_real_samples(self): return len(self.samples)
    def count_loss(self):
        reads_lost = sum(map(len, self.plexed.samples)) - sum(map(len, self.samples))
        percent = 100.0 * reads_lost / sum(map(len, self.plexed.samples))
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

    ############## Predictions ##############
    def predictions_table(self):
        output = ""
        for i, p in enumerate(self.pools):
            barcodes = p.guess_barcodes().most_common(len(p.samples)+2)
            output += "\n\n### %s:\n\n" % p.name
            for k,v in barcodes: output += "* %s: %s\n" % (k, v)
        return output
