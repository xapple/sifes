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
plexed = sifes.load("~/deploy/sifes/metadata/json/unige/desalt_plexed/")
proj = sifes.load("~/deploy/sifes/metadata/json/unige/desalt/")
demultiplexer = Demultiplexer(plexed, proj.samples)
