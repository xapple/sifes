#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""
A script to load settings and import modules.
"""

# Built-in modules #
import shutil

# Internal modules #
import sifes
from sifes.taxonomy                  import mothur_classify
from sifes.report                    import clusters
from sifes.groups                    import cluster
from sifes.filtering                 import seq_filter
from sifes.demultiplex.demultiplexer import Demultiplexer
from sifes.distribute.bundle         import Bundle
from sifes.distribute.upload         import DropBoxRclone
from sifes.location.seqenv_graphs    import SeqenvHeatmap
from sifes.groups.cluster_graphs     import DiversityRegression

# First party modules #
from plumbing.timer     import Timer

# Third party modules #
from tqdm import tqdm

# Parallelization strategy #
from plumbing.processes import prll_map

###############################################################################
# Load multiplexed and real project #
plexed = sifes.load("~/deploy/sifes/metadata/json/projects/unige/foram_plexed/", raw_files_must_exist=True)
proj   = sifes.load("~/deploy/sifes/metadata/json/projects/unige/foram/",        raw_files_must_exist=False)
demultiplexer = Demultiplexer(plexed, proj)

###############################################################################
# Parameters #
sifes.filtering.seq_filter.SeqFilter.primer_mismatches = 0
sifes.filtering.seq_filter.SeqFilter.primer_max_dist   = 40
sifes.filtering.seq_filter.SeqFilter.min_read_length   = 140 - 8 - 8 - 19 - 20
sifes.filtering.seq_filter.SeqFilter.max_read_length   = 250 - 8 - 8 - 19 - 20
sifes.filtering.seq_filter.SeqFilter.search_for_region = "GACAG"

sifes.groups.samples.Sample.default_joiner = 'pandaseq'
sifes.report.clusters.ClusterReport.default_taxa_graph_levels  = (4, 5, 6)

sifes.groups.cluster.Cluster.default_centering                 = 'swarm'
sifes.groups.cluster.Cluster.default_taxonomy                  = 'qiime'
sifes.taxonomy.qiime_classify.QiimeClassify.default_database   = 'foraminifera'
sifes.taxonomy.mothur_classify.MothurClassify.default_database = 'foraminifera'

sifes.location.seqenv_graphs.SeqenvHeatmap.custom_metadata      = 'depth'
sifes.groups.cluster_graphs.DiversityRegression.custom_metadata = 'dist_exhaust'