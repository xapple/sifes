#!/usr/bin/env python2
# -*- coding: utf-8 -*-

# Built-in modules #

# Internal modules #
import sifes
from sifes.demultiplex.demultiplexer import Demultiplexer

# Third party modules #
from tqdm import tqdm

###############################################################################
plexed = sifes.load("~/deploy/sifes/metadata/json/projects/sinclair/testproj_plexed/", raw_files_must_exist=False)
proj   = sifes.load("~/deploy/sifes/metadata/json/projects/sinclair/testproj/", raw_files_must_exist=False)

demultiplexer = Demultiplexer(plexed, proj)

#demultiplexer.plexfiles[0].guess_barcodes()
#print demultiplexer.plexfiles[0].primers.fwd_seq
#print demultiplexer.plexfiles[0].primers.rev_seq
#
#demultiplexer.plexfiles[0].merge_lanes()
#demultiplexer.plexfiles[0].primer_statistics()

demultiplexer.run()
demultiplexer.report.generate()