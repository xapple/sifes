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

###############################################################################
# Load multiplexed project #
test_proj = sifes.load("~/deploy/sifes/metadata/json/projects/envonautics/test/")
