#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""
A script to run small snippets of code.

# Get report #
rsync -avz --update edna:/home/sinclair/SIFES/views/projects/unige/foram_plexed/report/report.pdf ~/Desktop/current_report.pdf; open ~/Desktop/current_report.pdf
"""

# Built-in modules #

# Internal modules #
import sifes
from sifes.demultiplex.demultiplexer import Demultiplexer

# First party modules #
from plumbing.processes import prll_map
from plumbing.timer     import Timer

# Third party modules #
from tqdm import tqdm

###############################################################################
# Load multiplexed and real project #
plexed = sifes.load("~/deploy/sifes/metadata/json/projects/unige/foram_plexed/", raw_files_must_exist=False)
proj   = sifes.load("~/deploy/sifes/metadata/json/projects/unige/foram/",        raw_files_must_exist=False)

###############################################################################
