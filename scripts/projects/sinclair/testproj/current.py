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
#proj.cluster.report.generate()