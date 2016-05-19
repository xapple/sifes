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

# Internal modules #
import sifes.filtering.seq_filter
from sifes.demultiplex.demultiplexer import Demultiplexer

# First party modules #
from plumbing.processes import prll_map

# Third party modules #
from tqdm import tqdm

###############################################################################
# Load multiplexed project #
plexed = sifes.load("~/deploy/sifes/metadata/json/projects/micans/micans_v6_exp1_plexed/")

# Load all real project #
projects = ['minican_4_5', 'posmic_olk', 'pseud_fluo', 'febex_dp', 'cosc_1']
projects = [sifes.load("~/deploy/sifes/metadata/json/projects/micans/" + p, False) for p in projects]
samples  = [s for p in projects for s in p]
samples.sort(key=lambda s: s.num)

# Demultiplex #
demultiplexer = Demultiplexer(plexed, samples)
demultiplexer.run()

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
prll_map(lambda s: s.joiner.run(cpus=1), samples)

# Filter #
sifes.filtering.seq_filter.SeqFilter.primer_mismatches = 0
sifes.filtering.seq_filter.SeqFilter.primer_max_dist   = 25
sifes.filtering.seq_filter.SeqFilter.min_read_length   = 100
sifes.filtering.seq_filter.SeqFilter.max_read_length   = 160
for s in tqdm(samples): s.filter.run()
prll_map(lambda s: s.filter.run(), samples)

# Cluster #
for proj in projects:
    proj.cluster.combine_reads()
    proj.cluster.otus.run()

# Make report #
for s in tqdm(proj): s.report.generate()
