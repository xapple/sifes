#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""
A script to run all the procedure on the hundred lakes project.
"""

# Built-in modules #
import shutil

# Internal modules #
import sifes.filtering.seq_filter

# First party modules #
from plumbing.processes import prll_map
from plumbing.timer     import Timer
from plumbing.autopaths import FilePath

# Third party modules #
from tqdm import tqdm

###############################################################################
# Load two projects #
proj = sifes.load("~/deploy/sifes/metadata/json/projects/uppsala_universitet/hundred_lakes_euk/")

# Parameters #
sifes.filtering.seq_filter.SeqFilter.primer_mismatches = 1
sifes.filtering.seq_filter.SeqFilter.primer_max_dist   = 70
sifes.filtering.seq_filter.SeqFilter.min_read_length   = 300
sifes.filtering.seq_filter.SeqFilter.max_read_length   = 560
for s in proj: s.default_joiner = 'pandaseq'

###############################################################################
# Get information for excel file #
for s in proj: print s.pair.fwd.count
for s in proj: print s.pair.rev.count
for s in proj: print s.pair.fwd.md5
for s in proj: print s.pair.rev.md5
for s in proj: print len(s.pair.fwd.first)

# Uncompress reads - 0h01 #
with Timer(): prll_map(lambda s: s.uncompressed_pair, proj)
for s in proj:
    assert s.uncompressed_pair.fwd.count == s.uncompressed_pair.rev.count
    assert s.pair.count                  == s.uncompressed_pair.rev.count

# Join reads - 0h0x #
with Timer(): prll_map(lambda s: s.joiner.run(cpus=1), proj)
for s in proj: print s.short_name, s.joiner.results.unassembled_percent

# Filter - 0h0x #
with Timer(): prll_map(lambda s: s.filter.run(), proj)
for s in proj: print s.short_name, s.filter.primers_fasta.count
for s in proj: print s.short_name, s.filter.n_base_fasta.count
for s in proj: print s.short_name, s.filter.length_fasta.count
for s in proj: print s.short_name, s.filter.renamed_fasta.count

# Cluster combine reads - 0h0x #
print proj.short_name, proj.cluster.num_dropped_samples
print proj.short_name, [s.clean.count for s in proj.cluster.good_samples]
print proj.short_name, [s.clean.count for s in proj.cluster.bad_samples]
print proj.short_name, proj.cluster.read_count_cutoff
with Timer(): proj.cluster.combine_reads()

# Make centers  - 0h15 #
with Timer(): proj.cluster.centering.run(cpus=32)
print proj.short_name, proj.cluster.centering.results.centers.count

# Taxonomy assignment - 0h0x #
with Timer(): proj.cluster.taxonomy.run(cpus=1)
print proj.short_name, len(proj.cluster.taxonomy.results.assignments)

# Make the OTU table - 0h0x #
with Timer(): proj.cluster.otu_table.run()

# Make the taxa tables - 0h0x #
with Timer(): proj.cluster.taxa_table.run()

# Make sample graphs - 0hxx #
def basic_graphs(s):
    s.pair.fwd.fastqc.run()
    s.pair.rev.fastqc.run()
    s.joiner.results.assembled.graphs.length_dist(rerun=True)
with Timer(): prll_map(basic_graphs, proj)

# Make diversity sample graphs - 0hxx #
def diversity_plot(s):
    s.graphs.chao1(rerun=True)
    s.graphs.ace(rerun=True)
    s.graphs.shannon(rerun=True)
    s.graphs.simpson(rerun=True)
with Timer(): prll_map(diversity_plot, proj)

# Make project graphs - 0h0x #
def otu_plot(p):
    p.cluster.otu_table.results.graphs.otu_sizes_dist(rerun=True)
    p.cluster.otu_table.results.graphs.otu_sums_graph(rerun=True)
    p.cluster.otu_table.results.graphs.sample_sums_graph(rerun=True)
    p.cluster.otu_table.results.graphs.cumulative_presence(rerun=True)
    p.cluster.reads.graphs.length_dist(rerun=True)
    for g in p.cluster.taxa_table.results.graphs.__dict__.values(): g(rerun=True)
    if len (p.cluster) < 2: return
    p.cluster.nmds_graph(rerun=True)
with Timer(): otu_plot(proj)

# Optionally clean cache #
for s in proj: s.report.cache_dir.remove()
for s in proj: s.report.cache_dir.create()
for s in proj: print FilePath(s.report.cache_dir + 'genera_table.pickle').remove()

# Make cluster reports - 0h0x #
with Timer(): proj.cluster.report.generate()

# Make sample reports - 0h0x #
with Timer(): prll_map(lambda s: s.report.generate(), proj)

# Bundle - 0h02 #
from sifes.distribute.bundle import Bundle
bundle = Bundle("under_ice", proj.samples)
with Timer(): bundle.run()

# Upload to CKAN #
from sifes.distribute.ckan_samples import CkanSamples
ckan = CkanSamples(proj.samples)
ckan.run()