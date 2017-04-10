#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""
A script to run small snippets of code.
"""

# Built-in modules #

# Internal modules #
import sifes
from sifes.demultiplex.demultiplexer import Demultiplexer

# Third party modules #
from tqdm import tqdm

###############################################################################
plexed        = sifes.load("~/deploy/sifes/metadata/json/projects/unige/desalt_plexed/", raw_files_must_exist=False)
proj          = sifes.load("~/deploy/sifes/metadata/json/projects/unige/desalt/", raw_files_must_exist=False)
demultiplexer = Demultiplexer(plexed, proj)

demultiplexer.run()
demultiplexer.report.generate()
