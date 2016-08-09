#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""
A script to load the basic objects for the hundred lakes project.
"""

# Built-in modules #
import shutil

# First party modules #
from plumbing.processes import prll_map
from plumbing.timer import Timer

# Internal modules #
import sifes

###############################################################################
# Load two projects #
proj = sifes.load("~/deploy/sifes/metadata/json/projects/uppsala_universitet/hundred_lakes_euk/")

# Parameters #
sifes.filtering.seq_filter.SeqFilter.primer_mismatches = 1
sifes.filtering.seq_filter.SeqFilter.primer_max_dist   = 70
sifes.filtering.seq_filter.SeqFilter.min_read_length   = 300
sifes.filtering.seq_filter.SeqFilter.max_read_length   = 560
for s in proj: s.default_joiner = 'pandaseq'
