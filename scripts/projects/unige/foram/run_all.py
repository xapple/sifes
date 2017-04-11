#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""
A script to run all the procedure on our foram desalt project.

To first generate the JSON files:

    $ ~/repos/sifes/metadata/excel_to_json.py ~/repos/sifes/metadata/excel/projects/unige/foram_plexed/metadata_plexed.xlsx
    $ ~/repos/sifes/metadata/excel_to_json.py ~/repos/sifes/metadata/excel/projects/unige/foram/metadata.xlsx

To clean everything up:

    $ rm -rf ~/SIFES/raw/projects/unige/foram
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
sifes.filtering.seq_filter.SeqFilter.min_read_length   = 140
sifes.filtering.seq_filter.SeqFilter.max_read_length   = 250
for s in proj: s.default_joiner = 'pandaseq'

# Demultiplex - 0h08 #
demultiplexer = Demultiplexer(plexed, proj)
with Timer(): demultiplexer.run()

# Demultiplex Report #
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

print("# Filter - 0hxx #")
with Timer(): prll_map(lambda s: s.filter.run(), proj)
for s in proj: print s.short_name, s.filter.primers_fasta.count
for s in proj: print s.short_name, s.filter.n_base_fasta.count
for s in proj: print s.short_name, s.filter.length_fasta.count
for s in proj: print s.short_name, s.filter.renamed_fasta.count

###############################################################################
print("# Make sample reports - 0h0x #")
with Timer(): prll_map(lambda s: s.report.generate(), proj)

###############################################################################
print("# Bundle - 0h02 #")
from sifes.distribute.bundle import Bundle
bundle = Bundle("desalt_37f", proj.samples)
with Timer(): bundle.run()

print("# Extra files #")
path = sifes.home + "deploy/sifes/metadata/excel/projects/unige/foram/metadata.xlsx"
shutil.copy(path, bundle.p.samples_xlsx)

print("# Upload - 0h02 #")
from sifes.distribute.upload import DropBoxRclone
dbx_sync = DropBoxRclone(bundle.base_dir, '/Desalt 37F delivery')
with Timer(): dbx_sync.run()
print("Total delivery: %s" % bundle.base_dir.size)
