#!/usr/bin/env python2
# -*- coding: utf-8 -*-

# Modules #
from geopy.distance import vincenty

# Load #
import os
home = os.environ.get('HOME', '~') + '/'
execfile(home + "deploy/sifes/scripts/projects/unige/desalt/load.py")

# Constant #
sl = """as1a
as1b
as1c
as2a
as2b
as2c
as5a
as5b
as5c
vmso1a
vmso1b
vmso1c
vmso2a
vmso2b
vmso2c
vmso3a
vmso3b
vmso3c
had1a
had1b
had1c
had2a
had2b
had2c
had3a
had3b
had3c"""

# Calculate distance to exhaust #
exhausts = {'Ashqelon': (31.63212,  34.5167),
            'Soreq':    (31.941850, 34.687533),
            'Hadera':   (32.46502,  34.882583)}

# Loop #
names = set(s.short_name for s in proj)
for n in sl.split('\n'):
    if n not in names:
        print ''
        continue
    s = proj[n]
    exhaust  = exhausts[s.grouping]
    location = (float(s.info['latitude'][0]), float(s.info['longitude'][0]))
    result = vincenty(exhaust, location)
    print "%.2f" % result.meters

