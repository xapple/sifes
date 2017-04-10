#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""
A script to run all the procedure on our foram desalt project.

To first generate the JSON files:

    $ ~/repos/sifes/metadata/excel_to_json.py ~/repos/sifes/metadata/excel/projects/unige/foram_plexed/metadata_plexed.xlsx
    $ ~/repos/sifes/metadata/excel_to_json.py ~/repos/sifes/metadata/excel/projects/unige/foram/metadata.xlsx

To clean everything up:

    $ rm -rf ~/SIFES/raw/projects/unige/desalt
    $ rm -rf ~/SIFES/raw/projects/unige/foram
    $ rm -rf ~/SIFES/views/projects/unige/desalt_plexed
    $ rm -rf ~/SIFES/views/projects/unige/foram_plexed
    $ rm -rf ~/SIFES/views/samples/unige/desalt
    $ rm -rf ~/SIFES/views/samples/unige/foram
    $ rm -rf ~/SIFES/views/samples/unige/desalt_plexed
    $ rm -rf ~/SIFES/views/samples/unige/foram_plexed

"""

# Built-in modules #
import shutil

# Internal modules #
import sifes.filtering.seq_filter
from sifes.demultiplex.demultiplexer import Demultiplexer

# First party modules #
from plumbing.processes import prll_map
from plumbing.timer     import Timer
from plumbing.autopaths import FilePath

# Third party modules #
from tqdm import tqdm

###############################################################################
# Load multiplexed and real project #
plexed = sifes.load("~/deploy/sifes/metadata/json/projects/unige/foram_plexed/")
proj   = sifes.load("~/deploy/sifes/metadata/json/projects/unige/foram/")

# Demultiplex - 0h08 #
demultiplexer = Demultiplexer(plexed, proj)
with Timer(): demultiplexer.run()

# Demultiplex Report #
with Timer(): demultiplexer.report.generate()

###############################################################################
