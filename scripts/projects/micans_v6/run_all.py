#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""
A script to run all the procedure on Mican's first data job.
The code is: micans_v6_exp1
"""

# Built-in modules #

# Internal modules #
import sifes.filtering.seq_filter
from sifes.demultiplex.demultiplexer import Demultiplexer

# Third party modules #
from tqdm import tqdm

###############################################################################
# Load multiplexed project #
plexed = sifes.load("~/repos/sifes/metadata/json/projects/micans_v6_exp1_plexed/")

# Load real project #
proj = sifes.load("~/repos/sifes/metadata/json/projects/micans_v6_exp1/")

# Demultiplex #
demultipler = Demultiplexer(plexed, proj)
demultipler.run()

###############################################################################
# Get information for excel file #
for s in proj: print s.pair.fwd.count
for s in proj: print s.pair.rev.count
for s in proj: print s.pair.fwd.md5
for s in proj: print s.pair.rev.md5
for s in proj: print s.seq_len

# Join reads #
for s in tqdm(proj): s.joiner.run(cpus=1)

# Filter #
sifes.filtering.seq_filter.SeqFilter.primer_mismatches = 0
sifes.filtering.seq_filter.SeqFilter.primer_max_dist   = 25
sifes.filtering.seq_filter.SeqFilter.min_read_length   = 100
sifes.filtering.seq_filter.SeqFilter.max_read_length   = 160
for s in tqdm(proj): s.filter.run()

# Cluster #
proj.cluster.combine_reads()
proj.cluster.otus.run()

# Make report #
for s in tqdm(proj): s.report.generate()
