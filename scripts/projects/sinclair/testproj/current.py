#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""
A script to run small snippets of code.

# Get report #
rsync -avz --update edna:/home/sinclair/SIFES/views/projects/sinclair/testproj_plexed/report/report.pdf ~/Desktop/current_report.pdf; open ~/Desktop/current_report.pdf
rsync -avz --update edna:/home/sinclair/SIFES/views/samples/sinclair/testproj/as1a/report/report.pdf ~/Desktop/current_report.pdf; open ~/Desktop/current_report.pdf
rsync -avz --update edna:/home/sinclair/SIFES/views/projects/sinclair/testproj/cluster/testproj/report/report.pdf ~/Desktop/current_report.pdf; open ~/Desktop/current_report.pdf

"""

# Built-in modules #
import shutil

# Internal modules #
import sifes
from sifes.taxonomy import mothur_classify
from sifes.report   import clusters
from sifes.groups   import cluster
from sifes.demultiplex.demultiplexer import Demultiplexer

# First party modules #
from plumbing.processes import prll_map
from plumbing.timer     import Timer

# Third party modules #
from tqdm import tqdm

###############################################################################
# Load multiplexed and real project #
plexed = sifes.load("~/deploy/sifes/metadata/json/projects/sinclair/testproj_plexed/", raw_files_must_exist=False)
proj   = sifes.load("~/deploy/sifes/metadata/json/projects/sinclair/testproj/",        raw_files_must_exist=False)
demultiplexer = Demultiplexer(plexed, proj)

###############################################################################
# Parameters #
sifes.filtering.seq_filter.SeqFilter.primer_mismatches = 0
sifes.filtering.seq_filter.SeqFilter.primer_max_dist   = 50
sifes.filtering.seq_filter.SeqFilter.min_read_length   = 370 - 8 - 8 - 21 - 18
sifes.filtering.seq_filter.SeqFilter.max_read_length   = 450 - 8 - 8 - 21 - 18

sifes.groups.samples.Sample.default_joiner = 'pandaseq'
sifes.report.clusters.ClusterReport.default_taxa_graph_levels  = (3, 4, 5)

sifes.groups.cluster.Cluster.default_taxonomy = 'qiime'
sifes.taxonomy.qiime_classify.QiimeClassify.default_database = 'pr_two'
sifes.taxonomy.mothur_classify.MothurClassify.default_database = 'pr_two'

###############################################################################
#demultiplexer = Demultiplexer(plexed, proj)
#demultiplexer.run()
#demultiplexer.report.generate()

###############################################################################
#proj.first.filter.run()
#proj.first.report.generate()

###############################################################################
#proj.cluster.report.generate()