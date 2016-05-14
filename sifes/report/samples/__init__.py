# Futures #
from __future__ import division

# Built-in modules #
import os, re, shutil, socket
from collections import OrderedDict

# Internal modules #
import sifes
from sifes.reporting import ReportTemplate

# First party modules #
from plumbing.common import split_thousands
from pymarktex import Document, HeaderTemplate, FooterTemplate
from pymarktex.figures import ScaledFigure, DualFigure

# Third party modules #
import pandas
from tabulate import tabulate

# Constants #
ssh_header = "ssh://" + os.environ.get("FILESYSTEM_HOSTNAME", socket.getfqdn())

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

    copy_base = property(lambda self: sifes.reports_dir + self.sample.project.name + '/' + self.sample.short_name + '.pdf')

###############################################################################
class SampleTemplate(ReportTemplate):
    """All the parameters to be rendered in the markdown template"""
    delimiters = (u'{{', u'}}')

    def __repr__(self): return '<%s object on %s>' % (self.__class__.__name__, self.parent)

    def __init__(self, report):
        # Attributes #
        self.report, self.parent = report, report
        self.sample     = self.parent.sample
        self.project       = self.sample.project

    # General information #
    def sample_short_name(self):     return self.sample.short_name
    def sample_long_name(self):      return self.sample.short_name
    def project_short_name(self):    return self.sample.project_short_name
    def project_long_name(self):     return self.sample.project_long_name
    def project_other_samples(self): return len(self.project) - 1

    # JSON #
    def json_url(self):
        url = "https://github.com/limno/sifes/tree/master/json/samples"
        url += "/run%03d/run%03d-sample%03d.json"
        return url % (self.sample.run_num, self.sample.run_num, self.sample.num)
    def json_content(self):
        content = self.sample.json_path.read('utf-8')
        content = re.sub('\A(.+?)^    },$', '', content, flags=re.M|re.DOTALL)
        return content.strip('\n }')

    # Process info #
    def results_directory(self): return ssh_header + self.sample.base_dir

    # Raw data #
    def fwd_size(self):  return str(self.sample.fwd.size)
    def fwd_count(self): return split_thousands(self.sample.fwd.count)
    def fwd_qual(self):  return "%.2f" % self.sample.fwd.avg_quality
    def rev_size(self):  return str(self.sample.rev.size)
    def rev_count(self): return split_thousands(self.sample.rev.count)
    def rev_qual(self):  return "%.2f" % self.sample.rev.avg_quality
    def illumina_report(self): return self.sample.run.html_report_path
    def per_base_qual(self):
        params = [self.sample.fastq.fwd.fastqc.results.per_base_qual,
                  self.sample.fastq.rev.fastqc.results.per_base_qual]
        params += ["Forward", "Reverse"]
        params += ["fwd_per_base_qual", "rev_fwd_per_base_qual"]
        params += ["Per base quality", "per_base_qual"]
        return str(DualFigure(*params))
    def per_seq_qual(self):
        params = [self.sample.fastq.fwd.fastqc.results.per_seq_qual,
                  self.sample.fastq.rev.fastqc.results.per_seq_qual]
        params += ["Forward", "Reverse"]
        params += ["fwd_per_seq_qual", "rev_per_seq_qual"]
        params += ["Per sequence quality", "per_seq_qual"]
        return str(DualFigure(*params))

    # Joining #
    def joiner_version(self, val): self.sample.joiner.long_name
    def values_with_percent(self, val):
        percentage = lambda x,y: (len(x)/len(y))*100 if len(y) != 0 else 0
        percent = percentage(val, self.sample)
        return "%s (%.1f%%)" % (split_thousands(len(val)), percent)
    def assembled_count(self): return self.values_with_percent(self.sample.assembled)
    def unassembled_count(self): return self.values_with_percent(self.sample.unassembled)
    def low_qual_count(self):
        count = self.sample.assembled.stats['lowqual']
        return "%s (%.1f%%)" % (split_thousands(count), (count/len(self.sample))*100)
    def assembly_len_dist(self):
        caption = "Distribution of sequence lengths after joining"
        path = self.sample.assembled.length_dist.path
        label = "assembly_len_dist"
        return str(ScaledFigure(path, caption, label))
    def joined_quality(self):
        params = [self.sample.assembled.fastqc.results.per_base_qual,
                  self.sample.assembled.fastqc.results.per_seq_qual]
        params += ["Per base", "Per sequence"]
        params += ["joined_per_base_qual", "joined_per_seq_qual"]
        params += ["Joined sequence quality", "joined_quality"]
        return str(DualFigure(*params))

    # Filtering #
    def mismatches_allowed(self): return self.sample.assembled.primer_mismatches
    def primer_discard(self):
        before = self.sample.assembled
        after  = self.sample.assembled.good_primers.orig_reads
        return split_thousands(len(before) - len(after))
    def primer_left(self):
        return split_thousands(len(self.sample.assembled.good_primers.orig_reads))

    def n_base_discard(self):
        good = self.sample.assembled.good_primers
        return split_thousands(len(good.orig_reads) - len(good.n_filtered))
    def n_base_left(self):
        return split_thousands(len(self.sample.assembled.good_primers.n_filtered))

    def window_size(self): return self.sample.assembled.good_primers.qual_windowsize
    def window_threshold(self): return self.sample.assembled.good_primers.qual_threshold
    def window_discard(self):
        good = self.sample.assembled.good_primers
        return split_thousands(len(good.n_filtered) - len(good.qual_filtered))
    def window_left(self):
        return split_thousands(len(self.sample.assembled.good_primers.qual_filtered))

    def length_threshold(self): return self.sample.assembled.good_primers.min_length
    def length_discard(self):
        good = self.sample.assembled.good_primers
        return split_thousands(len(good.qual_filtered) - len(good.len_filtered))
    def length_left(self):
        return split_thousands(len(self.sample.assembled.good_primers.len_filtered))

    def percent_remaining(self):
        good = self.sample.assembled.good_primers
        return "%.1f%%" % ((len(good.len_filtered)/len(self.sample))*100)

    # Taxonomy #
    def abundant_table(self):
        # The data #
        row = self.sample.counts
        frame = pandas.DataFrame(index=range(len(row)))
        frame['Rank']  = range(1, len(row)+1)
        frame['Clade'] = row.index
        frame['Reads'] = [split_thousands(r) for r in row.values]
        frame['OTUs'] = [self.sample.project.cluster.otus.taxonomy.comp_tips.count_otus(s) for s in row.index]
        frame = frame[0:20]
        # Make it as text #
        table = tabulate(OrderedDict(frame), headers="keys", numalign="right", tablefmt="pipe")
        # Add caption #
        return table + "\n\n   : The 20 most abundant species in this sample."

    # Diversity #
    def total_otu_sum(self): return split_thousands(sum(self.sample.counts))
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
