# Futures #
from __future__ import division

# Built-in modules #
import shutil
from collections import OrderedDict

# Internal modules #
import sifes
from sifes.report import ReportTemplate

# First party modules #
from plumbing.common import split_thousands
from pymarktex import Document

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

    def generate(self):
        # Dynamic templates #
        self.main = MultiplexTemplate(self)
        self.markdown = unicode(self.main)
        # Render to latex #
        self.make_body()
        self.make_latex()
        self.make_pdf(safe=True)
        # Copy to reports directory #
        shutil.copy(self.output_path, self.copy_base)
        # Return #
        return self.output_path

    copy_base = property(lambda self: sifes.reports_dir + self.plexed.short_name + '/' + self.parent.short_name + '.pdf')

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
        self.pools   = self.parent.parent.pools
        self.samples = self.parent.parent.samples

    ############## General information ##############
    def plex_short_name(self): return self.plexed.short_name
    def plex_long_name(self):  return self.plexed.long_name

    ############## Input ##############
    def count_input_files(self):  return len(self.plexed)
    def count_pools(self):        return len(self.pools)

    def input_table(self):
        # The columns #
        info = OrderedDict((
            ('Name',  lambda p: "**" + p.name + "**"),
            ('Files', lambda p: "`" + str(len(p.inputs)) + "`"),
            ('Reads', lambda p: split_thousands(p.pair.count)),
            ('Loss',  lambda p: "%.1f%%" % (100 * sum(map(len, p.inputs)) / p.pair.count)),
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
        percent = 100 * reads_lost / sum(map(len, self.plexed.samples))
        return "%s (%.1f%%)" % (reads_lost, percent)

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