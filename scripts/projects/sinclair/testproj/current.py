#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""
A script to run small snippets of code.

# Get report #
rsync -avz --update edna:/home/sinclair/SIFES/views/projects/sinclair/testproj_plexed/report/report.pdf ~/Desktop/current_report.pdf; open ~/Desktop/current_report.pdf

To clean everything up:

    $ rm -rf ~/SIFES/raw/projects/sinclair/testproj
    $ rm -rf ~/SIFES/views/projects/sinclair/testproj_plexed
    $ rm -rf ~/SIFES/views/samples/sinclair/testproj
    $ rm -rf ~/SIFES/views/samples/sinclair/testproj_plexed
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
plexed = sifes.load("~/deploy/sifes/metadata/json/projects/sinclair/testproj_plexed/", raw_files_must_exist=False)
proj   = sifes.load("~/deploy/sifes/metadata/json/projects/sinclair/testproj/",        raw_files_must_exist=False)

###############################################################################
demultiplexer = Demultiplexer(plexed, proj)
demultiplexer.run()
demultiplexer.report.generate()
proj.first.joiner.run()
print proj.first.joiner.results.unassembled_percent