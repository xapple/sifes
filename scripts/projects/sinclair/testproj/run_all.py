#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""
A script to run all the procedure on our testproj project.

To first generate the JSON files:

    $ ~/repos/sifes/metadata/excel_to_json.py ~/repos/sifes/metadata/excel/projects/sinclair/testproj_plexed/metadata_plexed.xlsx
    $ ~/repos/sifes/metadata/excel_to_json.py ~/repos/sifes/metadata/excel/projects/sinclair/testproj/metadata.xlsx

To clean everything up:

    $ rm -rf ~/SIFES/raw/projects/sinclair/testproj
    $ rm -rf ~/SIFES/views/projects/sinclair/testproj_plexed
    $ rm -rf ~/SIFES/views/projects/sinclair/testproj
    $ rm -rf ~/SIFES/views/samples/sinclair/testproj
    $ rm -rf ~/SIFES/views/samples/sinclair/testproj_plexed
    
# Get report #
rsync -avz --update edna:/home/sinclair/SIFES/views/projects/sinclair/testproj_plexed/report/report.pdf ~/Desktop/current_report.pdf; open ~/Desktop/current_report.pdf
rsync -avz --update edna:/home/sinclair/SIFES/views/samples/sinclair/testproj/as1a/report/report.pdf ~/Desktop/current_report.pdf; open ~/Desktop/current_report.pdf
rsync -avz --update edna:/home/sinclair/SIFES/views/projects/sinclair/testproj/cluster/testproj/report/report.pdf ~/Desktop/current_report.pdf; open ~/Desktop/current_report.pdf
"""

execfile("~/deploy/sifes/scripts/projects/sinclair/testproj/load.py")

###############################################################################
print("# Demultiplex #")
with Timer(): demultiplexer.run()

print("# Demultiplex Report #")
with Timer(): demultiplexer.report.generate()

###############################################################################
print("# Join reads #")
with Timer(): prll_map(lambda s: s.joiner.run(cpus=1), proj)
for s in proj: print s.short_name, s.joiner.results.unassembled_percent

print("# Make fastqc graphs #")
def fastqc_graphs(s):
    s.pair.fwd.fastqc.run()
    s.pair.rev.fastqc.run()
with Timer(): prll_map(fastqc_graphs, proj)

print("# Make length dist graphs #")
def len_dist_graphs(s): s.joiner.results.assembled.graphs.length_dist(rerun=True)
with Timer(): prll_map(len_dist_graphs, proj)

print("# Filter #")
with Timer(): prll_map(lambda s: s.filter.run(), proj)

print("# Cluster combine reads #")
print proj.short_name, proj.cluster.num_dropped_samples
print proj.short_name, [s.clean.count for s in proj.cluster.good_samples]
print proj.short_name, [s.clean.count for s in proj.cluster.bad_samples]
print proj.short_name, proj.cluster.read_count_cutoff
with Timer(): proj.cluster.combine_reads()

print("# Make centers #")
with Timer(): proj.cluster.centering.run(cpus=8)
print proj.short_name, proj.cluster.centering.results.centers.count

print("# Taxonomy assignment #")
with Timer(): proj.cluster.taxonomy.run(cpus=8)
print proj.short_name, len(proj.cluster.taxonomy.results.assignments)

print("# Make the OTU table #")
with Timer(): proj.cluster.otu_table.run()

print("# Make the taxa tables #")
with Timer(): proj.cluster.taxa_table.run()

print("# Make sample graphs #")
def sample_plots(s):
    s.graphs.chao1(rerun=True)
    s.graphs.ace(rerun=True)
    s.graphs.shannon(rerun=True)
    s.graphs.simpson(rerun=True)
    #s.graphs.location_map(rerun=True)
with Timer(): prll_map(sample_plots, proj)

print("# Make project graphs #")
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

###############################################################################
print("# Make cluster reports #")
with Timer(): proj.cluster.report.generate()

print("# Make sample reports #")
for s in proj: s.report.purge_cache()
with Timer(): prll_map(lambda s: s.report.generate(), proj)

###############################################################################
print("# Bundle and upload #")
from sifes.distribute.bundle import Bundle
bundle = Bundle("testproj", proj.samples)
with Timer(): bundle.run()

path = sifes.home + "deploy/sifes/metadata/excel/projects/sinclair/testproj/metadata.xlsx"
shutil.copy(path, bundle.p.samples_xlsx)
path = sifes.reports_dir + 'testproj_plexed/demultiplexer.pdf'
shutil.copy(path, bundle.p.demultiplexing_report)

from sifes.distribute.upload import DropBoxRclone
dbx_sync = DropBoxRclone(bundle.base_dir, '/Testproj delivery')
with Timer(): dbx_sync.run()
print("Total delivery: %s" % bundle.base_dir.size)