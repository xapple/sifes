#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""
A script to run small snippets of code.
"""

# Built-in modules #

# Internal modules #
import sifes

# Third party modules #
from tqdm import tqdm

###############################################################################
plexed = sifes.load("~/deploy/sifes/metadata/json/projects/micans/micans_v6_exp1_plexed/")
projects = ['minican_4_5', 'posmic_olk', 'pseud_fluo', 'febex_dp', 'cosc_1']
projects = [sifes.load("~/deploy/sifes/metadata/json/projects/micans/" + p + '/', False) for p in projects]
samples  = [s for p in projects for s in p]
from sifes.demultiplex.demultiplexer import Demultiplexer
demultiplexer = Demultiplexer(plexed, samples)
demultiplexer.run()
