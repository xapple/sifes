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
execfile(sifes.home + "deploy/sifes/scripts/projects/micans/micans_v6/load.py")
projects[0].first.report.generate()