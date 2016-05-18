#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""
A script to run small snippets of code.
"""

# Built-in modules #

# Internal modules #
import sifes.filtering.seq_filter

# Third party modules #
from tqdm import tqdm

###############################################################################
test_proj = sifes.load("~/deploy/sifes/metadata/json/projects/envonautics/test/")
sifes.filtering.seq_filter.SeqFilter.primer_mismatches = 1
sifes.filtering.seq_filter.SeqFilter.primer_max_dist   = 70
sifes.filtering.seq_filter.SeqFilter.min_read_length   = 400
sifes.filtering.seq_filter.SeqFilter.max_read_length   = 560
for s in tqdm(test_proj): s.filter.run()
