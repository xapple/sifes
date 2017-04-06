"""
Special module to demultiplex samples.
"""

# Futures #
from __future__ import division

# Built-in modules #
import sys
from collections import Counter
import multiprocessing

# Internal modules #
from sifes.report.demultiplex import MultiplexReport

# First party modules #
from plumbing.autopaths import DirectoryPath
from plumbing.cache     import property_cached
from plumbing.processes import prll_map
from fasta              import FASTQ, PairedFASTQ
from fasta.primers      import iupac_pattern

# Third party modules #
import regex
from shell_command import shell_output
from tqdm import tqdm
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq

# Constants #

###############################################################################
class Demultiplexer(object):
    # Attributes #
    short_name = "demultiplexer"

    # Parameters #
    primer_mismatches  = 2
    barcode_mismatches = 0

    def __repr__(self): return '<%s object with %i samples>' % (self.__class__.__name__, len(self.samples))

    def __init__(self, plexed, samples):
        # Attributes #
        self.plexed  = plexed
        self.samples = samples
        # Check #
        assert all([s.info.get('multiplex_group') for s in self.plexed])
        assert all([s.info.get('multiplexed_in')  for s in self.samples])
        # Functions #
        get_inputs  = lambda name: [i for i in self.plexed if i.info['multiplex_group'] == name]
        get_samples = lambda name: [s for s in self.samples if s.info['multiplexed_in'] == name]
        # Plexfiles #
        self.plexfiles = sorted(list(set([s.info['multiplex_group'] for s in self.plexed])))
        self.plexfiles = [PlexFile(name, get_inputs(name), get_samples(name)) for name in self.plexfiles]
        # Report #
        self.report = MultiplexReport(self)

    def run(self, cpus=None):
        # Number of cores #
        if cpus is None: cpus = min(multiprocessing.cpu_count(), 32)
        # Merge lanes together #
        map(lambda p: p.merge_lanes(), self.plexfiles) #, cpus)
        # Search for barcodes #
        map(self.extract_reads, self.plexfiles) #, cpus)

    def extract_reads(self, plexfile):
        """Extract all reads from a plex file."""
        # Message #
        print "Demultiplexing plexfile '%s'" % plexfile.fwd_path
        # Main loop #
        for f,r in tqdm(plexfile.pair.parse_primers(plexfile.primers, self.primer_mismatches)):
            print 1/0
            if f.fwd_start_pos: print "G"

    @property_cached
    def results(self):
        results = DemultiplexerResults(self)
        message = "You can't access results from demultiplexing before running the algorithm."
        if not results: raise Exception(message)
        return results

###############################################################################
class DemultiplexerResults(object):

    def __nonzero__(self): return bool(False)

    def __init__(self, parent):
        self.parent = parent

###############################################################################
class PlexFile(object):
    """A pair of FASTQs containing multiple samples distinguished by barcodes."""

    def __repr__(self): return '<%s object with %i samples>' % (self.__class__.__name__, len(self.samples))
    def __nonzero__(self): return bool(self.pair)

    def __init__(self, name, inputs, samples):
        # Check #
        assert all([s.info.get('multiplex_group') == name for s in inputs])
        # Attributes #
        self.name      = name
        self.inputs    = inputs
        self.samples   = samples
        self.primers   = inputs[0].primers
        self.project   = inputs[0].project
        self.base_dir  = DirectoryPath(self.project.p.lane_cat_dir + self.name + '/')
        self.fwd_path  = FASTQ(self.base_dir + 'fwd.fastq.gz')
        self.rev_path  = FASTQ(self.base_dir + 'rev.fastq.gz')
        self.fwd_files = [s.pair.fwd for s in self.inputs]
        self.rev_files = [s.pair.rev for s in self.inputs]
        self.pair      = PairedFASTQ(self.fwd_path, self.rev_path)
        # Link the input samples back #
        for s in inputs: s.plexfile = self

    def merge_lanes(self):
        """We got a few runs that had several lanes in the same sample directory.
        We want to `cat` all these to files called fwd.fastq.gz and rev.fastq.gz"""
        # Make output directory #
        self.pair.fwd.directory.create_if_not_exists()
        # Case when you don't need to merge lanes #
        if len(self.inputs) == 1:
            self.pair.fwd.link_from(self.inputs[0].pair.fwd)
            self.pair.rev.link_from(self.inputs[0].pair.rev)
            self.pair.fwd.count = self.inputs[0].pair.fwd.count
            self.pair.rev.count = self.inputs[0].pair.rev.count
            return
        # Do it #
        shell_output("zcat %s |gzip > %s" % (' '.join(self.fwd_files), self.pair.fwd))
        shell_output("zcat %s |gzip > %s" % (' '.join(self.rev_files), self.pair.rev))
        # Check #
        assert self.pair.fwd.first.id == self.pair.rev.first.id

    #-------------------------------------------------------------------------#
    def primer_statistics(self):
        """Print all primer statistics."""
        inputs = []
        for name in ("fwd", "rev"):
            fastq = getattr(self.pair, name)
            for primer in ("fwd", "rev"):
                seq = getattr(self.primers, primer + "_seq")
                normal = seq
                revsre = seq.reverse_complement().complement()
                compl  = seq.complement()
                revcpl = seq.reverse_complement()
                inputs += [(name, fastq, primer, "normal    ", normal)]
                inputs += [(name, fastq, primer, "reverse   ", revsre)]
                inputs += [(name, fastq, primer, "complement", compl)]
                inputs += [(name, fastq, primer, "rev compl.", revcpl)]
        results = prll_map(self.get_primer_count, inputs, verbose=False)
        print '\n'.join(results)

    def get_primer_count(self, args):
        name, fastq, primer, sense, seq = args
        mismatches = 2
        pattern = regex.compile("(%s){s<=%i}" % (iupac_pattern(seq), mismatches))
        count   = 0
        for r in tqdm(fastq):
            if pattern.search(str(r.seq)): count += 1
        message  = name + " file, " + primer +  " primer, " + sense + ": "
        message += "%.1f%%" % (100.0 * (count / len(fastq)))
        message += " (%s)" % seq
        return message
