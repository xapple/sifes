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
import shutil

# Internal modules #
import sifes.filtering.seq_filter
from sifes.demultiplex.demultiplexer import Demultiplexer

# First party modules #
from plumbing.processes import prll_map
from plumbing.autopaths import DirectoryPath
from plumbing.timer import Timer

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

# Uncompress reads #
for s in samples:
    assert s.uncompressed_pair.fwd.count == s.uncompressed_pair.rev.count
    assert s.pair.count                  == s.uncompressed_pair.rev.count

# Join reads #
with Timer(): prll_map(lambda s: s.joiner.run(cpus=1), samples)
for s in samples: print s.short_name, s.joiner.results.unassembled_percent

# Filter #
sifes.filtering.seq_filter.SeqFilter.primer_mismatches = 0
sifes.filtering.seq_filter.SeqFilter.primer_max_dist   = 25
sifes.filtering.seq_filter.SeqFilter.min_read_length   = 60
sifes.filtering.seq_filter.SeqFilter.max_read_length   = 140
with Timer(): prll_map(lambda s: s.filter.run(), samples)
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

# Make sample graphs #
for s in samples:
    s.pair.fwd.fastqc.run()
    s.pair.rev.fastqc.run()

# Make good sample graphs #
for s in p.cluster:
    s.graphs.chao1()
    s.graphs.ace()
    s.graphs.shannon()
    s.graphs.simpson()

# Make project graphs #
for p in projects:
    p.cluster.otu_table.results.graphs.otu_sizes_dist()
    p.cluster.otu_table.results.graphs.otu_sums_graph()
    p.cluster.otu_table.results.graphs.sample_sums_graph()
    p.cluster.otu_table.results.graphs.cumulative_presence()
    for g in p.cluster.taxa_table.results.graphs.__dict__.values(): g()
    if len (p.cluster) < 2: continue
    p.cluster.nmds_graph()

# Clean cache #
for s in samples: s.report.cache_dir.remove()
for s in samples: s.report.cache_dir.create()

# The average quality is super long #
from sifes.report.samples import SampleTemplate
prll_map(lambda s: SampleTemplate(s.report).fwd_qual, samples)
prll_map(lambda s: SampleTemplate(s.report).rev_qual, samples)

# Make reports #
prll_map(lambda s: s.report.generate(), samples)
for tqdm(s) in samples: s.report.generate()
for p in projects: p.cluster.report.generate()

# Bundle #
from sifes.distribute.bundle import Bundle
bundle = Bundle("micans_v6", samples)
bundle.run()

# Extra files #
path = sifes.home + "deploy/sifes/metadata/excel/projects/micans/micans_v6_exp1/metadata.xlsx"
shutil.copy(path, bundle.p.samples_xlsx)
path = sifes.home + "deploy/sifes/metadata/excel/projects/micans/micans_v6_exp1_plexed/metadata_plexed.xlsx"
shutil.copy(path, bundle.p.multiplexed)
path = sifes.reports_dir + 'micans_v6_exp1_plexed/demultiplexer.pdf'
shutil.copy(path, bundle.p.demultiplexing_report)

# Upload #
from sifes.distribute.upload import DropBoxUpload
dbx_upload = DropBoxUpload(bundle.base_dir)
dbx_upload.run()