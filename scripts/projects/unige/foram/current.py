#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""
A script to run small snippets of code.
"""

import os
home = os.environ.get('HOME', '~') + '/'
execfile(home + "deploy/sifes/scripts/projects/unige/foram/load.py")

###############################################################################