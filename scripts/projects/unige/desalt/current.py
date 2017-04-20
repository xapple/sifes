#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""
A script to run small snippets of code.
"""

execfile("~/deploy/sifes/scripts/projects/unige/desalt/load.py")

###############################################################################

for g in p.cluster.taxa_table.results.graphs.by_rank: g(rerun=True)
with Timer(): proj.cluster.report.generate()

###############################################################################
print("# Bundle and upload - 0h0x #")
bundle = Bundle("desalt_v1v2", proj.samples)
with Timer(): bundle.run()

path = sifes.home + "deploy/sifes/metadata/excel/projects/unige/desalt/metadata.xlsx"
shutil.copy(path, bundle.p.samples_xlsx)
path = sifes.reports_dir + 'desalt_plexed/demultiplexer.pdf'
shutil.copy(path, bundle.p.demultiplexing_report)

dbx_sync = DropBoxRclone(bundle.base_dir, '/Desalt V1V2 delivery')
with Timer(): dbx_sync.run()
print("Total delivery: %s" % bundle.base_dir.size)
