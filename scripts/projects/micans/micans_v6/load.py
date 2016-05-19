#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""
A script to load the basic objects.
"""

# Built-in modules #

# First party modules #
from plumbing.processes import prll_map

# Internal modules #
import sifes
from sifes.demultiplex.demultiplexer import Demultiplexer

###############################################################################
# Load multiplexed project #
plexed = sifes.load("~/deploy/sifes/metadata/json/projects/micans/micans_v6_exp1_plexed/")

# Load all real project #
projects = ['minican_4_5', 'posmic_olk', 'pseud_fluo', 'febex_dp', 'cosc_1', 'mysterious']
projects = [sifes.load("~/deploy/sifes/metadata/json/projects/micans/" + p, False) for p in projects]
samples  = [s for p in projects for s in p]
samples.sort(key=lambda s: s.num)

# Demultiplex #
demultiplexer = Demultiplexer(plexed, samples)

# Parameters #
sifes.filtering.seq_filter.SeqFilter.primer_mismatches = 0
sifes.filtering.seq_filter.SeqFilter.primer_max_dist   = 25
sifes.filtering.seq_filter.SeqFilter.min_read_length   = 100
sifes.filtering.seq_filter.SeqFilter.max_read_length   = 160

