#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""
A script to run all the procedure on Mican's first data job.
The code is: micans_v6_exp1

To first generate the JSON files:

    $ ~/repos/sifes/metadata/excel_to_json.py ~/repos/sifes/metadata/excel/projects/micans/micans_v6_exp1_plexed/metadata_plexed.xlsx
    $ ~/repos/sifes/metadata/excel_to_json.py ~/repos/sifes/metadata/excel/projects/micans/micans_v6_exp1/metadata.xlsx

"""

# Built-in modules #
import shutil

# Internal modules #
import sifes.filtering.seq_filter
from sifes.demultiplex.demultiplexer import Demultiplexer

# First party modules #
from plumbing.processes import prll_map
from plumbing.timer     import Timer
from plumbing.autopaths import FilePath

# Third party modules #
from tqdm import tqdm

###############################################################################
# Load multiplexed project #
plexed = sifes.load("~/deploy/sifes/metadata/json/projects/micans/micans_v6_exp1_plexed/")

# Load all real project #
projects = ['minican_4_5', 'posmic_olk', 'pseud_fluo', 'febex_dp', 'cosc_1', 'mysterious']
projects = [sifes.load("~/deploy/sifes/metadata/json/projects/micans/" + p, False) for p in projects]
samples  = [s for p in projects for s in p]
samples.sort(key=lambda s: s.num)

# Demultiplex - 1h00 #
demultiplexer = Demultiplexer(plexed, samples)
with Timer(): demultiplexer.run()

# Demultiplex Report #
with Timer(): demultiplexer.report.generate()

###############################################################################
# Get information for excel file #
for s in samples: print s.pair.fwd.count
for s in samples: print s.pair.rev.count
for s in samples: print s.pair.fwd.md5
for s in samples: print s.pair.rev.md5
for s in samples: print len(s.pair.fwd.first)

# Uncompress reads - 0h15 #
with Timer(): prll_map(lambda s: s.uncompressed_pair, samples)
for s in samples:
    assert s.uncompressed_pair.fwd.count == s.uncompressed_pair.rev.count
    assert s.pair.count                  == s.uncompressed_pair.rev.count

# Join reads - 1h30 #
with Timer(): prll_map(lambda s: s.joiner.run(cpus=1), samples)
for s in samples: print s.short_name, s.joiner.results.unassembled_percent

# Filter - 0h20 #
sifes.filtering.seq_filter.SeqFilter.primer_mismatches = 0
sifes.filtering.seq_filter.SeqFilter.primer_max_dist   = 25
sifes.filtering.seq_filter.SeqFilter.min_read_length   = 55
sifes.filtering.seq_filter.SeqFilter.max_read_length   = 140
with Timer(): prll_map(lambda s: s.filter.run(), samples)
for s in samples: print s.short_name, s.filter.primers_fasta.count
for s in samples: print s.short_name, s.filter.n_base_fasta.count
for s in samples: print s.short_name, s.filter.length_fasta.count
for s in samples: print s.short_name, s.filter.renamed_fasta.count

# Cluster combine reads - 0h01 #
for p in projects: print p.short_name, p.cluster.num_dropped_samples
for p in projects: print p.short_name, [s.clean.count for s in p.cluster.good_samples]
for p in projects: print p.short_name, [s.clean.count for s in p.cluster.bad_samples]
for p in projects: print p.short_name, p.cluster.read_count_cutoff
with Timer(): prll_map(lambda p: p.cluster.combine_reads(), projects)

# Make centers  - 0h12 #
with Timer(): [p.cluster.centering.run(cpus=32) for p in projects]
for p in projects: print p.short_name, p.cluster.centering.results.centers.count

# Taxonomy assignment - 0h03 #
with Timer(): prll_map(lambda p: p.cluster.taxonomy.run(cpus=1), projects)
for p in projects: print p.short_name, len(p.cluster.taxonomy.results.assignments)

# Make the OTU table - 0h03 #
with Timer(): [p.cluster.otu_table.run() for p in projects]

# Make the taxa tables - 0h01 #
with Timer(): [p.cluster.taxa_table.run() for p in projects]

# Make sample graphs - 0hxx #
def basic_graphs(s):
    s.joiner.results.assembled.graphs.length_dist(rerun=True)
    s.pair.fwd.fastqc.run()
    s.pair.rev.fastqc.run()
with Timer(): prll_map(basic_graphs, samples)

# Make diversity sample graphs - 0h22 #
def diversity_plot(s):
    s.graphs.chao1(rerun=True)
    s.graphs.ace(rerun=True)
    s.graphs.shannon(rerun=True)
    s.graphs.simpson(rerun=True)
with Timer(): prll_map(diversity_plot, samples)

# Make project graphs - 0h04 #
def otu_plot(p):
    p.cluster.otu_table.results.graphs.otu_sizes_dist(rerun=True)
    p.cluster.otu_table.results.graphs.otu_sums_graph(rerun=True)
    p.cluster.otu_table.results.graphs.sample_sums_graph(rerun=True)
    p.cluster.otu_table.results.graphs.cumulative_presence(rerun=True)
    p.cluster.reads.graphs.length_dist(rerun=True)
    for g in p.cluster.taxa_table.results.graphs.__dict__.values(): g(rerun=True)
    if len (p.cluster) < 2: return
    p.cluster.nmds_graph(rerun=True)
with Timer(): prll_map(otu_plot, projects)

# Optionally clean cache #
for s in samples: s.report.cache_dir.remove()
for s in samples: s.report.cache_dir.create()
for s in samples: print FilePath(s.report.cache_dir + 'genera_table.pickle').remove()

# The average quality is super long - 0h10 #
from sifes.report.samples import SampleTemplate
with Timer(): prll_map(lambda s: SampleTemplate(s.report).fwd_qual, samples)
with Timer(): prll_map(lambda s: SampleTemplate(s.report).rev_qual, samples)

# Make cluster reports - 0h02 #
with Timer(): prll_map(lambda p: p.cluster.report.generate(), projects)

# Make sample reports - 0h15 # slowest: irb_bc_60_16_10_2
with Timer(): prll_map(lambda s: s.report.generate(), samples)

# Bundle - 0h02 #
from sifes.distribute.bundle import Bundle
bundle = Bundle("micans_v6", samples)
with Timer(): bundle.run()

# Extra files #
path = sifes.home + "deploy/sifes/metadata/excel/projects/micans/micans_v6_exp1/metadata.xlsx"
shutil.copy(path, bundle.p.samples_xlsx)
path = sifes.home + "deploy/sifes/metadata/excel/projects/micans/micans_v6_exp1_plexed/metadata_plexed.xlsx"
shutil.copy(path, bundle.p.multiplexed)
path = sifes.reports_dir + 'micans_v6_exp1_plexed/demultiplexer.pdf'
shutil.copy(path, bundle.p.demultiplexing_report)

# Upload - 0h03 #
from sifes.distribute.upload import DropBoxUpload
dbx_upload = DropBoxUpload(bundle.base_dir, '/Micans V6 analysis delivery')
with Timer(): dbx_upload.run()