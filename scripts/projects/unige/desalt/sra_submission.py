#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import sifes
from sifes.groups.projects import Projects

desalt = sifes.load("~/deploy/sifes/metadata/json/projects/unige/desalt/")
foram  = sifes.load("~/deploy/sifes/metadata/json/projects/unige/foram/")

lump = Projects("brine_disposal", [desalt, foram])
lump.sra.write_bio_tsv()
lump.sra.write_sra_tsv()