#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""
A script to run all the procedure on our V1V2 desalt project.

To first generate the JSON files:

    $ ~/repos/sifes/metadata/excel_to_json.py ~/repos/sifes/metadata/excel/projects/unige/desalt_plexed/metadata_plexed.xlsx
    $ ~/repos/sifes/metadata/excel_to_json.py ~/repos/sifes/metadata/excel/projects/unige/desalt/metadata.xlsx

To clean everything up:

    $ rm -rf ~/SIFES/raw/projects/unige/desalt
    $ rm -rf ~/SIFES/views/projects/unige/desalt
    $ rm -rf ~/SIFES/views/projects/unige/desalt_plexed
    $ rm -rf ~/SIFES/views/samples/unige/desalt
    $ rm -rf ~/SIFES/views/samples/unige/desalt_plexed

To run everything:

    time ~/deploy/sifes/scripts/other/running_runner.py ~/deploy/sifes/scripts/projects/unige/desalt/run_all.py --start_at=0
    
To get report:

    rsync -avz --update edna:/home/sinclair/SIFES/views/projects/unige/desalt_plexed/report/report.pdf ~/Desktop/current_report.pdf; open ~/Desktop/current_report.pdf
    rsync -avz --update edna:/home/sinclair/SIFES/views/samples/unige/desalt/as1a/report/report.pdf ~/Desktop/current_report.pdf; open ~/Desktop/current_report.pdf
    rsync -avz --update edna:/home/sinclair/SIFES/views/projects/unige/desalt/cluster/desalt/report/report.pdf ~/Desktop/current_report.pdf; open ~/Desktop/current_report.pdf
    
"""

import os
home = os.environ.get('HOME', '~') + '/'
execfile(home + "deploy/sifes/scripts/projects/unige/desalt/load.py")

###############################################################################
print("# Demultiplex - 0h55 #")
with Timer(): demultiplexer.run()

print("# Demultiplex Report #")
with Timer(): demultiplexer.report.generate()

###############################################################################
print("# Join reads - 0h03 #")
with Timer(): prll_map(lambda s: s.joiner.run(cpus=1), proj)
for s in proj: print s.short_name, s.joiner.results.unassembled_percent

print("# Make fastqc graphs - 0h02 #")
def fastqc_graphs(s):
    s.pair.fwd.fastqc.run()
    s.pair.rev.fastqc.run()
with Timer(): prll_map(fastqc_graphs, proj)

print("# Make length dist graphs - 0h08 #")
def len_dist_graphs(s): s.joiner.results.assembled.graphs.length_dist(rerun=True)
with Timer(): prll_map(len_dist_graphs, proj)

print("# Filter - 0h11 #")
with Timer(): prll_map(lambda s: s.filter.run(), proj)

print("# Cluster combine reads - 0h01 #")
print proj.short_name, proj.cluster.num_dropped_samples
print proj.short_name, [s.clean.count for s in proj.cluster.good_samples]
print proj.short_name, [s.clean.count for s in proj.cluster.bad_samples]
print proj.short_name, proj.cluster.read_count_cutoff
with Timer(): proj.cluster.combine_reads()

print("# Make centers  - 0h05 #")
with Timer(): proj.cluster.centering.run(cpus=8)
print proj.short_name, proj.cluster.centering.results.centers.count

print("# Taxonomy assignment - 0h0x #")
with Timer(): proj.cluster.taxonomy.run(cpus=8)
print proj.short_name, len(proj.cluster.taxonomy.results.assignments)

print("# Make the OTU table - 0h01 #")
with Timer(): proj.cluster.otu_table.run()

print("# Make the taxa tables - 0h01 #")
with Timer(): proj.cluster.taxa_table.run()

print("# Make sample graphs - 0h05 #")
def sample_plots(s):
    s.graphs.chao1(rerun=True)
    s.graphs.ace(rerun=True)
    s.graphs.shannon(rerun=True)
    s.graphs.simpson(rerun=True)
    #s.graphs.location_map(rerun=True)
with Timer(): prll_map(sample_plots, proj)

print("# Make project graphs - 0h04 #")
def otu_plot(p):
    p.cluster.otu_table.results.graphs.otu_sizes_dist(rerun=True)
    p.cluster.otu_table.results.graphs.otu_sums_graph(rerun=True)
    p.cluster.otu_table.results.graphs.sample_sums_graph(rerun=True)
    p.cluster.otu_table.results.graphs.cumulative_presence(rerun=True)
    p.cluster.reads.graphs.length_dist(rerun=True)
    for g in p.cluster.taxa_table.results.graphs.by_rank: g(rerun=True)
    #for g in p.cluster.locations_maps: g(rerun=True)
    if len (p.cluster) < 2: return
    p.cluster.nmds_graph(rerun=True)
with Timer(): otu_plot(proj)

print("# Make sub taxa graphs - 0h0x #")
with Timer():
    map(lambda sub:sub.run(), proj.cluster.sub_taxa_tables)
    graphs = (g for sub in proj.cluster.sub_taxa_tables for g in sub.results.graphs.by_rank)
    map(lambda g:g(rerun=True), graphs)

###############################################################################
print("# Make cluster reports - 0h01 #")
with Timer(): proj.cluster.report.generate()

print("# Make sample reports - 0h0x #")
for s in proj: s.report.purge_cache()
with Timer(): prll_map(lambda s: s.report.generate(), proj)

###############################################################################
print("# Bundle and upload - 0h0x #")
bundle   = Bundle("desalt_v1v2", proj.samples)
dbx_sync = DropBoxRclone(bundle.base_dir, '/Desalt V1V2 delivery')

with Timer():
    bundle.run()
    path = sifes.home + "deploy/sifes/metadata/excel/projects/unige/desalt/metadata.xlsx"
    shutil.copy(path, bundle.p.samples_xlsx)
    path = sifes.reports_dir + 'desalt_plexed/demultiplexer.pdf'
    shutil.copy(path, bundle.p.demultiplexing_report)
    dbx_sync.run()

print("Total delivery: %s" % bundle.base_dir.size)
