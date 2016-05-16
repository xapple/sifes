# Futures #
from __future__ import division

# Built-in modules #
import shutil

# Internal modules #
import sifes
from sifes.report import ReportTemplate

# First party modules #
from pymarktex import Document, HeaderTemplate, FooterTemplate

# Third party modules #

###############################################################################
class MultiplexReport(Document):
    """A full report generated in PDF for every Sample object."""

    def __repr__(self): return '<%s object on %s>' % (self.__class__.__name__, self.parent)

    def __init__(self, demultiplexer):
        # Attributes #
        self.demultiplexer, self.parent = demultiplexer, demultiplexer
        # Automatic paths #
        self.base_dir    = self.demultiplexer.p.report_dir
        self.output_path = self.demultiplexer.p.report_pdf

    def generate(self):
        # Dynamic templates #
        self.markdown = unicode(MultiplexTemplate(self))
        self.header = HeaderTemplate()
        self.footer = FooterTemplate()
        # Render to latex #
        self.make_body()
        self.make_latex()
        self.make_pdf(safe=True)
        # Copy to reports directory #
        shutil.copy(self.output_path, self.copy_base)
        # Return #
        return self.output_path

    copy_base = property(lambda self: sifes.reports_dir + self.parent.project.name + '/' + self.parent.short_name + '.pdf')

###############################################################################
class MultiplexTemplate(ReportTemplate):
    """All the parameters to be rendered in the markdown template"""
    delimiters = (u'{{', u'}}')

    def __repr__(self): return '<%s object on %s>' % (self.__class__.__name__, self.parent)

    def __init__(self, report):
        # Attributes #
        self.report  = report
        self.sample  = self.parent.sample
        self.project = self.sample.project

    ############## General information ##############
    def sample_short_name(self):     return self.sample.short_name
    def sample_long_name(self):      return self.sample.short_name
    def project_short_name(self):    return self.sample.project_short_name
    def project_long_name(self):     return self.sample.project_long_name
    def project_other_samples(self): return len(self.project) - 1