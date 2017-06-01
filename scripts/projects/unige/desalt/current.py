#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""
A script to run small snippets of code.
"""

import os
home = os.environ.get('HOME', '~') + '/'
execfile(home + "deploy/sifes/scripts/projects/unige/desalt/load.py")

###############################################################################
if False:
    proj.cluster.graphs.diversity_reg_chao1(rerun=True)
    proj.cluster.graphs.diversity_reg_ace(rerun=True)
    proj.cluster.graphs.diversity_reg_shannon(rerun=True)
    proj.cluster.graphs.diversity_reg_simpson(rerun=True)

    with Timer(): proj.cluster.report.generate()

###############################################################################
if False:
    print("# Bundle and upload - 0h0x #")
    bundle   = Bundle("desalt_v1v2", proj.samples)
    dbx_sync = DropBoxRclone(bundle.base_dir, '/Desalt V1V2 delivery')

    with Timer():
        bundle.run()
        path = sifes.home + "deploy/sifes/metadata/excel/projects/unige/desalt/metadata.xlsx"
        shutil.copy(path, bundle.p.samples_xlsx)
        path = sifes.reports_dir + 'desalt_plexed/demultiplexer.pdf'
        shutil.copy(path, bundle.p.demultiplexing_report)
        dbx_sync.run()

    print("Total delivery: %s" % bundle.base_dir.size)

###############################################################################
if True:
    print("# Run redundancy - 0h00 #")
    proj.cluster.redundancy.run()
