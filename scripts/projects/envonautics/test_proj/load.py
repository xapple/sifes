#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""
A script to load the basic objects.
"""

# Built-in modules #

# First party modules #
from plumbing.processes import prll_map
from plumbing.timer import Timer

# Internal modules #
import sifes

###############################################################################
# Load multiplexed project #
test_proj = sifes.load("~/deploy/sifes/metadata/json/projects/envonautics/test/")

# Parameters #
sifes.filtering.seq_filter.SeqFilter.primer_mismatches = 1
sifes.filtering.seq_filter.SeqFilter.primer_max_dist   = 70
sifes.filtering.seq_filter.SeqFilter.min_read_length   = 400
sifes.filtering.seq_filter.SeqFilter.max_read_length   = 560
