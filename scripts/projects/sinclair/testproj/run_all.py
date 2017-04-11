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
    $ rm -rf ~/SIFES/views/samples/sinclair/testproj
    $ rm -rf ~/SIFES/views/samples/sinclair/testproj_plexed
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
plexed = sifes.load("~/deploy/sifes/metadata/json/projects/sinclair/testproj_plexed/", raw_files_must_exist=False)
proj   = sifes.load("~/deploy/sifes/metadata/json/projects/sinclair/testproj/",        raw_files_must_exist=False)

# Parameters #
sifes.filtering.seq_filter.SeqFilter.primer_mismatches = 0
sifes.filtering.seq_filter.SeqFilter.primer_max_dist   = 50
sifes.filtering.seq_filter.SeqFilter.min_read_length   = 370
sifes.filtering.seq_filter.SeqFilter.max_read_length   = 450
for s in proj: s.default_joiner = 'pandaseq'

# Demultiplex - 0h55 #
demultiplexer = Demultiplexer(plexed, proj)
with Timer(): demultiplexer.run()

# Demultiplex Report #
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
for s in proj: print s.short_name, s.filter.primers_fasta.count
for s in proj: print s.short_name, s.filter.n_base_fasta.count
for s in proj: print s.short_name, s.filter.length_fasta.count
for s in proj: print s.short_name, s.filter.renamed_fasta.count

###############################################################################
print("# Make sample reports #")
with Timer(): prll_map(lambda s: s.report.generate(), proj)

###############################################################################
print("# Bundle #")
from sifes.distribute.bundle import Bundle
bundle = Bundle("testproj", proj.samples)
with Timer(): bundle.run()

print("# Extra files #")
path = sifes.home + "deploy/sifes/metadata/excel/projects/sinclair/testproj/metadata.xlsx"
shutil.copy(path, bundle.p.samples_xlsx)

print("# Upload #")
from sifes.distribute.upload import DropBoxRclone
dbx_sync = DropBoxRclone(bundle.base_dir, '/Testproj delivery')
with Timer(): dbx_sync.run()
print("Total delivery: %s" % bundle.base_dir.size)
