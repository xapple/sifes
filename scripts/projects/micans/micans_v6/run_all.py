#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""
A script to run all the procedure on Mican's first data job.
The code is: micans_v6_exp1
"""

# Built-in modules #

# Internal modules #
import sifes
import sifes.filtering.seq_filter
from sifes.demultiplex.demultiplexer import Demultiplexer

# Third party modules #
from tqdm import tqdm

###############################################################################
# Load multiplexed project #
plexed = sifes.load("~/deploy/sifes/metadata/json/projects/micans/micans_v6_exp1_plexed/")

# Load all real project #
projects = ['minican_4_5', 'posmic_olk', 'pseud_fluo', 'febex_dp', 'cosc_1']
projects = [sifes.load("~/deploy/sifes/metadata/json/projects/micans/" + p, False) for p in projects]
samples  = [s for p in projects for s in p]

# Demultiplex #
demultiplexer = Demultiplexer(plexed, samples)
demultiplexer.run()

# Guess barcodes #
barcodes_guessed = demultiplexer.pools[1].guess_barcodes(stop_at=40000)
for k,v in barcodes_guessed.most_common(20): print k, ':', v

# Demultiplex Report #
demultiplexer.report.generate()

###############################################################################
# Get information for excel file #
for s in samples: print s.pair.fwd.count
for s in samples: print s.pair.rev.count
for s in samples: print s.pair.fwd.md5
for s in samples: print s.pair.rev.md5
for s in samples: print len(s.pair.fwd.first_read)

# Join reads #
for s in tqdm(samples): s.joiner.run(cpus=1)

# Filter #
sifes.filtering.seq_filter.SeqFilter.primer_mismatches = 0
sifes.filtering.seq_filter.SeqFilter.primer_max_dist   = 25
sifes.filtering.seq_filter.SeqFilter.min_read_length   = 100
sifes.filtering.seq_filter.SeqFilter.max_read_length   = 160
for s in tqdm(samples): s.filter.run()

# Cluster #
for proj in projects:
    proj.cluster.combine_reads()
    proj.cluster.otus.run()

# Make report #
for s in tqdm(proj): s.report.generate()
