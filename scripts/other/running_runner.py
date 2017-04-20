#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""
A script that tears down the ipython environment between each step of the "run_all.py" files.
"""

# Modules #
import os, sys, re, argparse
from argparse import RawTextHelpFormatter

# Internal modules #

# First party modules #
from plumbing.autopaths import FilePath
from plumbing.tmpstuff  import new_temp_dir
from plumbing.color     import Color

# Third party modules #
import sh

# Make a shell arguments parser #
parser = argparse.ArgumentParser(description="Running runner", formatter_class=RawTextHelpFormatter)

# All the required arguments #
parser.add_argument("script_path", help="The path to file to process", type=str)

# All the optional arguments #
parameters = {
    "start_at"      : "A single integer.",
    "include"       : "A list of comma separated integers.",
}

# Parse it #
for param, hlp in parameters.items(): parser.add_argument("--" + param, help=hlp)
args        = parser.parse_args()
script_path = args.script_path
start_at    = int(args.start_at)                if args.start_at else None
include     = map(int, args.include.split(',')) if args.include  else None

###############################################################################
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
else: sections.append(current)

# Get the header #
header = sections.pop(0)

# Make a temporary directory #
directory = new_temp_dir()

# Restore them every time #
out, err = sys.stdout, sys.stderr

# Run every section #
for i, section in enumerate(sections):
    # Optional filters #
    if start_at and i < start_at:     continue
    if include  and i not in include: continue
    # Do it #
    path = directory + "section_%i.py" % i
    path.writelines(header + ['\n\n\n\n\n\n\n'] + section)
    print Color.f_ylw + "-- Section %i/%i --" % (i, len(sections)-1) + Color.end
    sh.ipython(path, _out=sys.stdout, _err=sys.stderr)
    sys.stdout.flush()
    sys.stderr.flush()

# Clean up #
print "killall firefox geckodriver Xvfb"