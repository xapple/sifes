#!/usr/bin/env python2
# -*- coding: utf-8 -*-

# Internal modules #
import sifes

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

# Load all real project #
proj = sifes.load("~/deploy/sifes/metadata/json/projects/unige/desalt/")
names = set(s.short_name for s in proj)

# Loop #
for n in sl.split('\n'):
    if n not in names:
        print ''
        continue
    s = proj[n]
    print len(s.pair.fwd.first)

