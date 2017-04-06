#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""
A script to run create a fake test project with a few sequences
"""

# Built-in modules #
from itertools import izip, islice

# Internal modules #
import sifes

# Third party modules #
from tqdm import tqdm

###############################################################################
# Load the project #
test_proj = sifes.load("~/deploy/sifes/metadata/json/projects/envonautics/test/", False)

# Load the sequences to copy #
import illumitag
illumi_proj = illumitag.projects['ice']
illumi_proj.load()
illumi_proj.cluster.load()

# Samples #
illumi_samples = illumi_proj[8:16]

# Copy #
for i_s, s in tqdm(izip(illumi_samples, test_proj)):
    s.pair.fwd.write(islice(i_s.fwd, 10000))
    s.pair.rev.write(islice(i_s.rev, 10000))
