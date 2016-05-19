#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""
A script to run all the procedure on the test project.
The code is: test_proj

To first generate the JSON files:

    $ ~/repos/sifes/metadata/excel_to_json.py ~/repos/sifes/metadata/excel/projects/envonautics/test_proj/metadata.xlsx

"""

# Built-in modules #

# Internal modules #
import sifes.filtering.seq_filter

# First party modules #
from plumbing.processes import prll_map

# Third party modules #
from tqdm import tqdm

###############################################################################
# Load multiplexed project #
test_proj = sifes.load("~/deploy/sifes/metadata/json/projects/envonautics/test/")

# Join reads #
prll_map(lambda s: s.joiner.run(cpus=1), test_proj)
for s in tqdm(test_proj): s.joiner.run(cpus=1)

# Filter #
sifes.filtering.seq_filter.SeqFilter.primer_mismatches = 1
sifes.filtering.seq_filter.SeqFilter.primer_max_dist   = 70
sifes.filtering.seq_filter.SeqFilter.min_read_length   = 400
sifes.filtering.seq_filter.SeqFilter.max_read_length   = 560
prll_map(lambda s: s.filter.run(), test_proj)
for s in tqdm(test_proj): s.filter.run()

# Cluster #
test_proj.cluster.combine_reads()
test_proj.cluster.otus.run()
test_proj.cluster.taxonomy.run(cpus=8)

# Make report #
for s in tqdm(test_proj): s.report.generate()
test_proj.report.generate()