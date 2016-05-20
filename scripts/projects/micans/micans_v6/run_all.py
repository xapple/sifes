#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""
A script to run all the procedure on Mican's first data job.
The code is: micans_v6_exp1

To first generate the JSON files:

    $ ~/repos/sifes/metadata/excel_to_json.py ~/repos/sifes/metadata/excel/projects/micans/micans_v6_exp1_plexed/metadata_plexed.xlsx
    $ ~/repos/sifes/metadata/excel_to_json.py ~/repos/sifes/metadata/excel/projects/micans/micans_v6_exp1/metadata.xlsx

"""

# Built-in modules #

# Internal modules #
import sifes.filtering.seq_filter
from sifes.demultiplex.demultiplexer import Demultiplexer

# First party modules #
from plumbing.processes import prll_map

# Third party modules #
from tqdm import tqdm

###############################################################################
# Load multiplexed project #
plexed = sifes.load("~/deploy/sifes/metadata/json/projects/micans/micans_v6_exp1_plexed/")

# Load all real project #
projects = ['minican_4_5', 'posmic_olk', 'pseud_fluo', 'febex_dp', 'cosc_1']
projects = [sifes.load("~/deploy/sifes/metadata/json/projects/micans/" + p, False) for p in projects]
samples  = [s for p in projects for s in p]
samples.sort(key=lambda s: s.num)

# Demultiplex #
demultiplexer = Demultiplexer(plexed, samples)
demultiplexer.run()

# Demultiplex Report #
demultiplexer.report.generate()

###############################################################################
# Get information for excel file #
for s in samples: print s.pair.fwd.count
for s in samples: print s.pair.rev.count
for s in samples: print s.pair.fwd.md5
for s in samples: print s.pair.rev.md5
for s in samples: print len(s.pair.fwd.first_read)

# Join reads #
prll_map(lambda s: s.joiner.run(cpus=1), samples)
for s in samples: print s.short_name, s.joiner.results.assembled.count

# Filter #
sifes.filtering.seq_filter.SeqFilter.primer_mismatches = 0
sifes.filtering.seq_filter.SeqFilter.primer_max_dist   = 25
sifes.filtering.seq_filter.SeqFilter.min_read_length   = 60
sifes.filtering.seq_filter.SeqFilter.max_read_length   = 140
prll_map(lambda s: s.filter.run(), samples)
for s in samples: print s.short_name, s.filter.primers_fasta.count
for s in samples: print s.short_name, s.filter.n_base_fasta.count
for s in samples: print s.short_name, s.filter.length_fasta.count
for s in samples: print s.short_name, s.filter.renamed_fasta.count

# Cluster #
for p in projects: print p.short_name, p.cluster.num_dropped_samples
for p in projects: print p.short_name, [s.clean.count for s in p.cluster.good_samples]
for p in projects: print p.short_name, [s.clean.count for s in p.cluster.bad_samples]
for p in projects: print p.short_name, p.cluster.reads.count
for proj in projects: proj.cluster.combine_reads()

# Make centers #
for proj in projects: p.cluster.centering.run(cpus=32)
for p in projects: print p.short_name, p.cluster.centering.results.centers.count

# Taxonomy assignment #
for p in projects: p.cluster.taxonomy.run(cpus=32)
for p in projects: print p.short_name, len(p.cluster.taxonomy.results.assignments)

# Make the OTU table #
for p in projects: p.cluster.otu_table.run()

# Make the taxa tables #
for p in projects: p.cluster.taxa_table.run()

# Make graphs #
for s in samples:
    s.pair.fwd.fastqc.run()
    s.pair.rev.fastqc.run()
    s.graphs.chao1()
    s.graphs.ace()
    s.graphs.shannon()
    s.graphs.simpson()
for p in projects:
    p.cluster.otu_table.results.graphs.otu_sizes_dist()
    p.cluster.otu_table.results.graphs.otu_sums_graph()
    p.cluster.otu_table.results.graphs.sample_sums_graph()
    p.cluster.otu_table.results.graphs.cumulative_presence()
    for g in p.cluster.taxa_table.results.graphs.__dict__.values(): g()
    p.cluster.nmds_graph()

# Make reports #
for s in samples:  s.report.generate()
for p in projects: p.cluster.report.generate()
