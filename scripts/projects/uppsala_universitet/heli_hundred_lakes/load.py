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
sifes.filtering.seq_filter.SeqFilter.primer_mismatches = 2
sifes.filtering.seq_filter.SeqFilter.primer_max_dist   = 35
sifes.filtering.seq_filter.SeqFilter.min_read_length   = 400
sifes.filtering.seq_filter.SeqFilter.max_read_length   = 530
for s in proj: s.default_joiner = 'pandaseq'
