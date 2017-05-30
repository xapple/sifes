#!/usr/bin/env python2
# -*- coding: utf-8 -*-

# Futures #
from __future__ import division

# Load #
import os
home = os.environ.get('HOME', '~') + '/'
execfile(home + "deploy/sifes/scripts/projects/unige/desalt/load.py")

# Percent unclassified at Phylum level for V1V2 #
table = proj.cluster.taxa_table.results.taxa_tables_by_rank[2].sum()
print table['Unassigned'] / table.sum()
