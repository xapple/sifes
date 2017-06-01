#!/usr/bin/env python2
# -*- coding: utf-8 -*-


# Load #
import os
home = os.environ.get('HOME', '~') + '/'
execfile(home + "deploy/sifes/scripts/projects/unige/desalt/load.py")


proj.cluster.redundancy.run()