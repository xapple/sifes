#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""
A script to load the basic objects for the under ice project.
"""

# Built-in modules #
import shutil

# First party modules #
from plumbing.processes import prll_map
from plumbing.timer import Timer

# Internal modules #
import sifes

###############################################################################
# Load multiplexed project #
proj = sifes.load("~/deploy/sifes/metadata/json/projects/envonautics/under_ice/")

# Parameters #
sifes.filtering.seq_filter.SeqFilter.primer_mismatches = 1
sifes.filtering.seq_filter.SeqFilter.primer_max_dist   = 70
sifes.filtering.seq_filter.SeqFilter.min_read_length   = 400
sifes.filtering.seq_filter.SeqFilter.max_read_length   = 560