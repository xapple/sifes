#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""
A script to load the basic objects for Mercè's temperature change project.
"""

# Built-in modules #
import shutil

# First party modules #
from plumbing.processes import prll_map
from plumbing.timer import Timer

# Internal modules #
import sifes

###############################################################################
# Load one project #
proj = sifes.load("~/deploy/sifes/metadata/json/projects/uppsala_universitet/merce_temp_change/")

# Parameters #
sifes.filtering.seq_filter.SeqFilter.primer_mismatches = 2
sifes.filtering.seq_filter.SeqFilter.primer_max_dist   = 35
sifes.filtering.seq_filter.SeqFilter.min_read_length   = 400
sifes.filtering.seq_filter.SeqFilter.max_read_length   = 530
for s in proj: s.default_joiner = 'pandaseq'

print "Hi"
proj.cluster.unifrac_matrix.run()