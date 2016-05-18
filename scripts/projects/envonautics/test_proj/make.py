#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""
A script to run create a fake test project with a few sequences
"""

# Built-in modules #

# Internal modules #
import sifes

# Third party modules #
from tqdm import tqdm

###############################################################################
# Load the project #
test_proj = sifes.load("~/deploy/sifes/metadata/json/projects/envonautics/test_proj/", False)

# Load the sequences to copy #
import illumitag
illumi_proj = illumitag.projects['ice']
illumi_proj.load()
illumi_proj.cluster.load()

# Copy #
for s in test_proj:
    pass