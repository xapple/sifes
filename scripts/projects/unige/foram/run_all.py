#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""
A script to run all the procedure on our foram desalt project.

To first generate the JSON files:

    $ ~/repos/sifes/metadata/excel_to_json.py ~/repos/sifes/metadata/excel/projects/unige/foram_plexed/metadata_plexed.xlsx
    $ ~/repos/sifes/metadata/excel_to_json.py ~/repos/sifes/metadata/excel/projects/unige/foram/metadata.xlsx

To clean everything up:

    $ rm -rf ~/SIFES/raw/projects/unige/foram
    $ rm -rf ~/SIFES/views/projects/unige/foram
    $ rm -rf ~/SIFES/views/projects/unige/foram_plexed
    $ rm -rf ~/SIFES/views/samples/unige/foram
    $ rm -rf ~/SIFES/views/samples/unige/foram_plexed

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
# Load multiplexed and real project #
plexed = sifes.load("~/deploy/sifes/metadata/json/projects/unige/foram_plexed/")
proj   = sifes.load("~/deploy/sifes/metadata/json/projects/unige/foram/")

# Parameters #
sifes.filtering.seq_filter.SeqFilter.primer_mismatches = 0
sifes.filtering.seq_filter.SeqFilter.primer_max_dist   = 40
sifes.filtering.seq_filter.SeqFilter.min_read_length   = 140 - 8 - 8 - 19 - 20
sifes.filtering.seq_filter.SeqFilter.max_read_length   = 250 - 8 - 8 - 19 - 20
sifes.groups.samples.Sample.default_joiner = 'pandaseq'
sifes.taxonomy.mothur_classify.MothurClassify.default_database = 'foraminifera'
sifes.report.clusters.ClusterReport.default_taxa_graph_levels  = (4, 5, 6)

print("# Demultiplex - 0h08 #")
demultiplexer = Demultiplexer(plexed, proj)
with Timer(): demultiplexer.run()

print("# Demultiplex Report #")
with Timer(): demultiplexer.report.generate()

###############################################################################
print("# Join reads - 0h01 #")
with Timer(): prll_map(lambda s: s.joiner.run(cpus=1), proj)
for s in proj: print s.short_name, s.joiner.results.unassembled_percent

print("# Make fastqc graphs - 0h01 #")
def fastqc_graphs(s):
    s.pair.fwd.fastqc.run()
    s.pair.rev.fastqc.run()
with Timer(): prll_map(fastqc_graphs, proj)

print("# Make length dist graphs - 0h01 #")
def len_dist_graphs(s): s.joiner.results.assembled.graphs.length_dist(rerun=True)
with Timer(): prll_map(len_dist_graphs, proj)

print("# Filter - 0h03 #")
with Timer(): prll_map(lambda s: s.filter.run(), proj)

print("# Cluster combine reads - 0h01 #")
print proj.short_name, proj.cluster.num_dropped_samples
print proj.short_name, [s.clean.count for s in proj.cluster.good_samples]
print proj.short_name, [s.clean.count for s in proj.cluster.bad_samples]
print proj.short_name, proj.cluster.read_count_cutoff
with Timer(): proj.cluster.combine_reads()

print("# Make centers  - 0h01 #")
with Timer(): proj.cluster.centering.run(cpus=8)
print proj.short_name, proj.cluster.centering.results.centers.count

print("# Taxonomy assignment - 0h0x #")
with Timer(): proj.cluster.taxonomy.run(cpus=8)
print proj.short_name, len(proj.cluster.taxonomy.results.assignments)

print("# Make the OTU table - 0h01 #")
with Timer(): proj.cluster.otu_table.run()

print("# Make the taxa tables - 0h01 #")
with Timer(): proj.cluster.taxa_table.run()

print("# Make sample graphs - 0h02 #")
def sample_plots(s):
    s.graphs.chao1(rerun=True)
    s.graphs.ace(rerun=True)
    s.graphs.shannon(rerun=True)
    s.graphs.simpson(rerun=True)
    s.graphs.location_map(rerun=True)
with Timer(): prll_map(sample_plots, proj)

print("# Make project graphs - 0h01 #")
def otu_plot(p):
    p.cluster.otu_table.results.graphs.otu_sizes_dist(rerun=True)
    p.cluster.otu_table.results.graphs.otu_sums_graph(rerun=True)
    p.cluster.otu_table.results.graphs.sample_sums_graph(rerun=True)
    p.cluster.otu_table.results.graphs.cumulative_presence(rerun=True)
    p.cluster.reads.graphs.length_dist(rerun=True)
    for g in p.cluster.taxa_table.results.graphs.by_rank: g(rerun=True)
    for g in p.cluster.locations_maps: g(rerun=True)
    if len (p.cluster) < 2: return
    p.cluster.nmds_graph(rerun=True)
with Timer(): otu_plot(proj)

###############################################################################
print("# Make cluster reports - 0h01 #")
with Timer(): proj.cluster.report.generate()

print("# Make sample reports  - 0h02 #")
for s in proj: s.report.purge_cache()
with Timer(): prll_map(lambda s: s.report.generate(), proj)

###############################################################################
print("# Bundle - 0h02 #")
from sifes.distribute.bundle import Bundle
bundle = Bundle("desalt_37f", proj.samples)
with Timer(): bundle.run()

print("# Extra files #")
path = sifes.home + "deploy/sifes/metadata/excel/projects/unige/foram/metadata.xlsx"
shutil.copy(path, bundle.p.samples_xlsx)
path = sifes.reports_dir + 'foram_plexed/demultiplexer.pdf'
shutil.copy(path, bundle.p.demultiplexing_report)

print("# Upload - 0h02 #")
from sifes.distribute.upload import DropBoxRclone
dbx_sync = DropBoxRclone(bundle.base_dir, '/Desalt 37F delivery')
with Timer(): dbx_sync.run()
print("Total delivery: %s" % bundle.base_dir.size)
