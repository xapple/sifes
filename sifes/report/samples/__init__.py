# Futures #
from __future__ import division

# Built-in modules #
import re, shutil
from collections import OrderedDict

# Internal modules #
import sifes
from sifes.report import ReportTemplate

# First party modules #
from pymarktex         import Document
from plumbing.common   import split_thousands
from pymarktex.figures import ScaledFigure, DualFigure

# Third party modules #
import pandas
from tabulate import tabulate

###############################################################################
class SampleReport(Document):
    """A full report generated in PDF for every Sample object."""

    def __repr__(self): return '<%s object on %s>' % (self.__class__.__name__, self.parent)

    def __init__(self, sample):
        # Attributes #
        self.sample, self.parent = sample, sample
        # Automatic paths #
        self.base_dir    = self.sample.p.report_dir
        self.output_path = self.sample.p.report_pdf

    def generate(self):
        # Dynamic templates #
        self.markdown = unicode(SampleTemplate(self))
        # Render to latex #
        self.make_body()
        self.make_latex()
        self.make_pdf(safe=True)
        # Copy to reports directory #
        shutil.copy(self.output_path, self.copy_base)
        # Return #
        return self.output_path

    copy_base = property(lambda self: sifes.reports_dir + self.sample.project.name + '/' + self.sample.short_name + '.pdf')

###############################################################################
class SampleTemplate(ReportTemplate):
    """All the parameters to be rendered in the markdown template"""
    delimiters = (u'{{', u'}}')

    def __repr__(self): return '<%s object on %s>' % (self.__class__.__name__, self.parent)

    def __init__(self, report):
        # Attributes #
        self.report  = report
        self.sample  = self.parent.sample
        self.project = self.sample.project

    def values_with_percent(self, val):
        percentage = lambda x,y: (len(x)/len(y))*100 if len(y) != 0 else 0
        percent    = percentage(val, self.sample)
        return "%s (%.1f%%)" % (split_thousands(len(val)), percent)

    ############## General information ##############
    def sample_short_name(self):     return self.sample.short_name
    def sample_long_name(self):      return self.sample.long_name
    def project_short_name(self):    return self.sample.project_short_name
    def project_long_name(self):     return self.sample.project_long_name
    def project_other_samples(self): return len(self.project) - 1

    ############## JSON ##############
    def json_content(self):
        content = self.sample.json_path.read('utf-8')
        #TODO: remove blank lines
        pass
        # Remove the contacts #
        content = re.sub('\A(.+?)^    },$', '', content, flags=re.M|re.DOTALL)
        # Remove the last brace #
        return content.strip('\n }')

    ############## Raw data ##############
    def fwd_size(self):  return str(self.sample.pair.fwd.size)
    def fwd_count(self): return split_thousands(self.sample.pair.fwd.count)
    def fwd_qual(self):  return "%.2f" % self.sample.pair.fwd.avg_quality
    def rev_size(self):  return str(self.sample.pair.rev.size)
    def rev_count(self): return split_thousands(self.sample.pair.rev.count)
    def rev_qual(self):  return "%.2f" % self.sample.pair.rev.avg_quality

    def per_base_qual(self):
        params = [self.sample.pair.fwd.fastqc.results.per_base_qual,
                  self.sample.pair.rev.fastqc.results.per_base_qual]
        params += ["Forward", "Reverse"]
        params += ["fwd_per_base_qual", "rev_fwd_per_base_qual"]
        params += ["Per base quality", "per_base_qual"]
        return str(DualFigure(*params))

    def per_seq_qual(self):
        params = [self.sample.pair.fwd.fastqc.results.per_seq_qual,
                  self.sample.pair.rev.fastqc.results.per_seq_qual]
        params += ["Forward", "Reverse"]
        params += ["fwd_per_seq_qual", "rev_per_seq_qual"]
        params += ["Per sequence quality", "per_seq_qual"]
        return str(DualFigure(*params))

    ############## Joining ##############
    def joiner_version(self): self.sample.joiner.long_name

    def assembled_count(self):
        return self.values_with_percent(self.sample.joiner.results.assembled)
    def unassembled_count(self):
        return self.values_with_percent(self.sample.joiner.results.unassembled)

    def joined_len_dist(self):
        caption = "Distribution of sequence lengths after joining"
        path    = self.sample.joiner.results.assembled.graphs.length_dist()
        label   = "joined_len_dist"
        return str(ScaledFigure(path, caption, label))

    ############## Filtering ##############
    def primer_max_dist(self):    return self.sample.filter.primer_max_dist
    def mismatches_allowed(self): return self.sample.filter.primer_mismatches
    def primer_discard(self):
        before = self.sample.joining.results.assembled
        after  = self.sample.filter.results.primers_fasta
        return split_thousands(len(before) - len(after))
    def primer_left(self):
        return split_thousands(len(self.sample.filter.results.primers_fasta))

    def n_base_discard(self):
        before = self.sample.filter.results.primers_fasta
        after  = self.sample.filter.results.n_base_fasta
        return split_thousands(len(before) - len(after))
    def n_base_left(self):
        return split_thousands(len(self.sample.filter.results.n_base_fasta))

    def min_read_length(self): return self.sample.filter.min_read_length
    def max_read_length(self): return self.sample.filter.max_read_length
    def length_discard(self):
        before = self.sample.filter.results.n_base_fasta
        after  = self.sample.filter.results.length_fasta
        return split_thousands(len(before) - len(after))
    def length_left(self):
        return split_thousands(len(self.sample.filter.results.length_fasta))

    def percent_remaining(self):
        return "%.1f%%" % ((len(self.sample.filter.results.clean)/len(self.sample))*100)

    ############## Taxonomy ##############
    def abundant_table(self):
        if not self.project.otus: return False
        # The data #
        row = self.sample.counts
        frame = pandas.DataFrame(index=range(len(row)))
        frame['Rank']  = range(1, len(row)+1)
        frame['Clade'] = row.index
        frame['Reads'] = [split_thousands(r) for r in row.values]
        frame['OTUs']  = [self.sample.project.cluster.otus.taxonomy.comp_tips.count_otus(s) for s in row.index]
        frame = frame[0:20]
        # Make it as text #
        table = tabulate(OrderedDict(frame), headers="keys", numalign="right", tablefmt="pipe")
        # Add caption #
        return table + "\n\n   : The 20 most abundant predicted species in this sample."

    ############## Diversity ##############
    def diversity(self):
        if not self.project.otus: return False
        params = ('total_otu_sum', 'total_otu_count', 'chao1_curve',
                  'ace_curve', 'shannon_curve', 'simpson_curve')
        return {p:getattr(self, p) for p in params}

    def total_otu_sum(self):   return split_thousands(sum(self.sample.counts))
    def total_otu_count(self): return split_thousands(len(self.sample.counts))
    def chao1_curve(self):
        caption = "Chao1 rarefaction curve"
        path = self.sample.diversity.chao1.path
        label = "chao1_curve"
        return str(ScaledFigure(path, caption, label))
    def ace_curve(self):
        caption = "Ace rarefaction curve"
        path = self.sample.diversity.ace.path
        label = "ace_curve"
        return str(ScaledFigure(path, caption, label))
    def shannon_curve(self):
        caption = "Shannon rarefaction curve"
        path = self.sample.diversity.shannon.path
        label = "shannon_curve"
        return str(ScaledFigure(path, caption, label))
    def simpson_curve(self):
        caption = "Simpson rarefaction curve"
        path = self.sample.diversity.simpson.path
        label = "simpson_curve"
        return str(ScaledFigure(path, caption, label))
