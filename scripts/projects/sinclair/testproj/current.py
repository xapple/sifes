#!/usr/bin/env python2
# -*- coding: utf-8 -*-

# Built-in modules #

# Internal modules #
import sifes
from sifes.demultiplex.demultiplexer import Demultiplexer

# Third party modules #
from tqdm import tqdm

###############################################################################
plexed        = sifes.load("~/deploy/sifes/metadata/json/projects/sinclair/testproj_plexed/", raw_files_must_exist=False)
proj          = sifes.load("~/deploy/sifes/metadata/json/projects/sinclair/testproj/", raw_files_must_exist=False)
demultiplexer = Demultiplexer(plexed, proj.samples)

demultiplexer.run()