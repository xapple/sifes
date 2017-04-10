#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""
A script to run all the procedure on our V1V2 desalt project.

To first generate the JSON files:

    $ ~/repos/sifes/metadata/excel_to_json.py ~/repos/sifes/metadata/excel/projects/unige/desalt_plexed/metadata_plexed.xlsx
    $ ~/repos/sifes/metadata/excel_to_json.py ~/repos/sifes/metadata/excel/projects/unige/desalt/metadata.xlsx

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
plexed = sifes.load("~/deploy/sifes/metadata/json/projects/unige/desalt_plexed/", raw_files_must_exist=False)
proj   = sifes.load("~/deploy/sifes/metadata/json/projects/unige/desalt/",        raw_files_must_exist=False)

# Demultiplex - xh00 #
demultiplexer = Demultiplexer(plexed, proj)
with Timer(): demultiplexer.run()

# Demultiplex Report #
with Timer(): demultiplexer.report.generate()

