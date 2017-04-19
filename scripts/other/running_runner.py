#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""
A script that tears down the ipython environment between each step of the "run_all.py" files.
"""

# Modules #
import os, sys, re

# Internal modules #

# First party modules #
from plumbing.autopaths import FilePath
from plumbing.tmpstuff  import new_temp_dir
from plumbing.color     import Color

# Third party modules #
import sh

###############################################################################
# Get the shell arguments #
if len(sys.argv) < 2: sys.exit(sys.modules[__name__].__doc__)
script_path = sys.argv[1]

# Check that the path is valid #
if not os.path.exists(script_path): raise Exception("No file at %s." % script_path)

# Load data #
script = FilePath(script_path)

# Cut into sections #
sections = []
current  = []
for line in script:
   if re.findall('\Aprint\("# (.+?) #"\)', line):
       sections.append(current)
       current = []
   current.append(line)

# Get the header #
header = sections.pop(0)

# Make a temporary directory #
directory = new_temp_dir()

# Restore them every time #
out, err = sys.stdout, sys.stderr

# Run every section #
for i, section in enumerate(sections):
    path = directory + "section_%i.py" % i
    path.writelines(header + ['\n\n\n\n\n\n\n'] + section)
    print Color.f_ylw + "-- Section %i --" % i + Color.end
    sh.ipython(path, _out=sys.stdout, _err=sys.stderr)
    sys.stdout.flush()
    sys.stderr.flush()

# Clean up #
print "killall firefox geckodriver Xvfb"