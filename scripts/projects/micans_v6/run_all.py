#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""
A script to run all the procedure on Mican's first data job.
The code is: micans_v6_exp1
"""

# Built-in modules #

# Internal modules #
import sifes

# Third party modules #
from tqdm import tqdm

###############################################################################
proj = sifes.load("~/repos/sifes/metadata/json/projects/micans_v6_exp1/")

# Get information for excel file #
for s in proj: print s.pair.fwd.count
for s in proj: print s.pair.rev.count
for s in proj: print s.pair.fwd.md5
for s in proj: print s.pair.rev.md5
for s in proj: print s.seq_len

# Join reads #
for s in tqdm(proj): s.joiner.run()

# Filter #
for s in tqdm(proj): s.filter.run()

# Cluster #
proj.cluster.combine_reads()
proj.cluster.otus.run()

# Make report #
for s in tqdm(proj): s.report.generate()
