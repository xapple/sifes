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
execfile(sifes.home + "deploy/sifes/scripts/projects/envonautics/under_ice/load.py")
from sifes.distribute.ckan_samples import CkanSamples
samples = proj[['bt1', 'rl1', 'lb1']]
ckan = CkanSamples(samples, groups={'under-ice': samples})
ckan.run()
