#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""
A script to run small snippets of code.
"""

import os
home = os.environ.get('HOME', '~') + '/'
execfile(home + "deploy/sifes/scripts/projects/sinclair/testproj/load.py")

###############################################################################
#demultiplexer = Demultiplexer(plexed, proj)
#demultiplexer.run()
#demultiplexer.report.generate()

###############################################################################
#proj.first.filter.run()
#proj.first.report.generate()

###############################################################################
#proj.cluster.seqenv.results.graphs.seqenv_heatmap(rerun=True)
#print proj.cluster.report.generate()
#proj.cluster.graphs.diversity_reg_ace(rerun=True)

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
    bundle = Bundle("testproj", proj.samples)
    dbx_sync = DropBoxRclone(bundle.base_dir, '/Testproj delivery')

    with Timer():
        bundle.run()
        path = sifes.home + "deploy/sifes/metadata/excel/projects/sinclair/testproj/metadata.xlsx"
        shutil.copy(path, bundle.p.samples_xlsx)
        path = sifes.reports_dir + 'testproj_plexed/demultiplexer.pdf'
        shutil.copy(path, bundle.p.demultiplexing_report)
        dbx_sync.run()

    print("Total delivery: %s" % bundle.base_dir.size)

###############################################################################
if True:
    import sifes
    from sifes.groups.projects import Projects
    desalt = sifes.load("~/deploy/sifes/metadata/json/projects/unige/desalt/")
    foram  = sifes.load("~/deploy/sifes/metadata/json/projects/unige/foram/")
    lump = Projects("brine_disposal", [desalt, foram])
    lump.sra.write_bio_tsv()
    lump.sra.write_sra_tsv()

