#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""
A script to run all the procedure on the hundred lakes project.
"""

import os

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
proj = sifes.load("~/deploy/sifes/metadata/json/projects/uppsala_universitet/merce_temp_change/")

# Parameters #
sifes.filtering.seq_filter.SeqFilter.primer_mismatches = 2
sifes.filtering.seq_filter.SeqFilter.primer_max_dist   = 35
sifes.filtering.seq_filter.SeqFilter.min_read_length   = 400
sifes.filtering.seq_filter.SeqFilter.max_read_length   = 530
for s in proj: s.default_joiner = 'pandaseq'

###############################################################################
print("# Get information for excel file #")
for s in proj: print s.short_name
for s in proj: print s.pair.fwd.count
for s in proj: print s.pair.rev.count
for s in proj: print s.pair.fwd.md5
for s in proj: print s.pair.rev.md5
for s in proj: print len(s.pair.fwd.first)

print("# Uncompress reads - 0h02 #")
with Timer(): prll_map(lambda s: s.uncompressed_pair, proj)
for s in proj:
    assert s.uncompressed_pair.fwd.count == s.uncompressed_pair.rev.count
    assert s.pair.count                  == s.uncompressed_pair.rev.count

print("# Join reads - 0h02 #")
with Timer(): prll_map(lambda s: s.joiner.run(cpus=1), proj)
for s in proj: print s.joiner.results.unassembled_percent

print("# Filter - 0h11 #")
with Timer(): prll_map(lambda s: s.filter.run(), proj)
for s in proj: print s.short_name, s.filter.primers_fasta.count
for s in proj: print s.short_name, s.filter.n_base_fasta.count
for s in proj: print s.short_name, s.filter.length_fasta.count
for s in proj: print s.short_name, s.filter.renamed_fasta.count

print("# Cluster combine reads - 0h01 #")
print proj.short_name, proj.cluster.num_dropped_samples
print proj.short_name, [s.clean.count for s in proj.cluster.good_samples]
print proj.short_name, [s.clean.count for s in proj.cluster.bad_samples]
print proj.short_name, proj.cluster.read_count_cutoff
with Timer(): proj.cluster.combine_reads()

print("# Make centers  - 0h16 #")
with Timer(): proj.cluster.centering.run(cpus=32)
print proj.short_name, proj.cluster.centering.results.centers.count

print("# Taxonomy assignment - 0h03 #")
with Timer(): proj.cluster.taxonomy.run(cpus=1)
print proj.short_name, len(proj.cluster.taxonomy.results.assignments)

print("# Make the OTU table - 0h01 #")
with Timer(): proj.cluster.otu_table.run()

print("# Make the taxa tables - 0h01 #")
with Timer(): proj.cluster.taxa_table.run()

print("# Make fastqc graphs - 0h50 #")
for s in tqdm(proj):
    print "fastqc101 on sample '%s'" % s
    s.pair.fwd.fastqc.run()
    s.pair.rev.fastqc.run()

print("# Make diversity sample graphs - 0h03 #")
def diversity_plot(s):
    s.graphs.chao1(rerun=True)
    s.graphs.ace(rerun=True)
    s.graphs.shannon(rerun=True)
    s.graphs.simpson(rerun=True)
with Timer(): prll_map(diversity_plot, proj)

print("# Make project graphs - 0h12 #")
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

print("# Make length dist graphs - 0h23 #")
for s in tqdm(proj): print s.joiner.results.assembled.graphs.length_dist(rerun=True)

#print("# Optionally clean cache #")
#for s in proj: s.report.cache_dir.remove()
#for s in proj: s.report.cache_dir.create()
#for s in proj: print FilePath(s.report.cache_dir + 'genera_table.pickle').remove()

print("# Make cluster reports - 0h05 #")
with Timer(): proj.cluster.report.generate()

print("# Attribute deletion because of odd pickling parallelization problem #")
for s in proj: del s.joiner.results.assembled.graphs

print("# Make sample reports - 0h05 #")
with Timer(): prll_map(lambda s: s.report.generate(), proj)

print("# Bundle - 0h02 #")
from sifes.distribute.bundle import Bundle
bundle = Bundle("merce_temp_change", proj.samples)
with Timer(): bundle.run()

print("# Extra files #")
path = sifes.home + "deploy/sifes/metadata/excel/projects/uppsala_universitet/merce_temp_change/metadata.xlsx"
shutil.copy(path, bundle.p.samples_xlsx)

print("# Upload - 0h02 #")
from sifes.distribute.upload import DropBoxSync
dbx_sync = DropBoxSync(bundle.base_dir, '/Mercè temperature change')
with Timer(): dbx_sync.run()
print("Total delivery: %s" % bundle.base_dir.size)
