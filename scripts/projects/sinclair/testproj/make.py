#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""
    $ ~/repos/sifes/metadata/excel_to_json.py ~/repos/sifes/metadata/excel/projects/sinclair/testproj_plexed/metadata_plexed.xlsx
    $ ~/repos/sifes/metadata/excel_to_json.py ~/repos/sifes/metadata/excel/projects/sinclair/testproj/metadata.xlsx
"""

# Built-in modules #
from itertools import izip, islice

# Internal modules #
import sifes

# Third party modules #
from tqdm import tqdm

###############################################################################
# Load the project #
test_proj = sifes.load("~/deploy/sifes/metadata/json/projects/sinclair/testproj_plexed/", False)
test_pair = test_proj[0].pair

# Load the sequences to copy #
real_proj = sifes.load("~/deploy/sifes/metadata/json/projects/unige/desalt_plexed/")
real_pair = real_proj[0].pair

# Copy #
test_pair.fwd.write(islice(real_pair.fwd, 10000))
test_pair.rev.write(islice(real_pair.rev, 10000))
