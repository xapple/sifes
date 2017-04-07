"""
Special module to demultiplex samples.
"""

# Futures #
from __future__ import division

# Built-in modules #
import sys
from collections import Counter, defaultdict, OrderedDict
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
    """A demultiplexer for when you have one barcode on each side"""

    # Attributes #
    short_name = "demultiplexer"

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
        map(lambda p: p.extract_reads(), self.plexfiles) #, cpus)

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
    """A pair of FASTQs (after optional merging of several pairs) containing
    multiple samples distinguished by custom barcodes."""

    # Parameters #
    primer_mismatches  = 0
    barcode_mismatches = 0

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

    def extract_reads(self):
        """Extract all reads from a plex file."""
        # Message #
        print "Demultiplexing plexfile '%s'" % self.fwd_path
        # Always the same barcode length #
        barlen = len(self.samples[0].info['forward_mid'])
        # Barcodes #
        fwd_barcodes = OrderedDict((s.info['forward_mid'], s.info['forward_num']) for s in self.samples)
        rev_barcodes = OrderedDict((s.info['reverse_mid'], s.info['reverse_num']) for s in self.samples)
        # Array of possibilities #
        array = defaultdict(lambda: defaultdict(int))
        for s in self.samples: array[s.info['forward_mid']][s.info['reverse_mid']] = s
        # Create samples #
        for s in self.samples: s.pair.create()
        # Dropped sequences #
        not_both_primers     = 0
        unknown_fwd_barcode  = 0
        unknown_rev_barcode  = 0
        # Main loop #
        for f,r in tqdm(self.pair.parse_primers(self.primers, self.primer_mismatches)):
            # Case primers not found #
            if not (f.fwd_match and r.rev_match) and not (r.fwd_match and f.rev_match):
                not_both_primers += 1
                continue
            # Regular case #
            if f.fwd_match and r.rev_match:
                fwd_barcode = str(f.read[f.fwd_start_pos-barlen:f.fwd_start_pos].seq)
                rev_barcode = str(r.read[r.rev_start_pos-barlen:r.rev_start_pos].seq)
                case = 'regular'
            # Goofy case #
            elif r.fwd_match and f.rev_match:
                fwd_barcode = str(r.read[r.fwd_start_pos-barlen:r.fwd_start_pos].seq)
                rev_barcode = str(f.read[f.rev_start_pos-barlen:f.rev_start_pos].seq)
                case = 'goofy'
            # Throw away unknown barcodes #
            if fwd_barcode not in fwd_barcodes: unknown_fwd_barcode += 1
            if rev_barcode not in rev_barcodes: unknown_rev_barcode += 1
            if fwd_barcode not in fwd_barcodes or rev_barcode not in rev_barcodes: continue
            # Get sample #
            s = array[fwd_barcode][rev_barcode]
            # Case it's a bad combination #
            if isinstance(s, int):
                array[fwd_barcode][rev_barcode] += 1
                continue
            # Case it's a good sample #
            if case == 'regular': s.pair.add(f.read, r.read.reverse_complement())
            if case == 'goofy':   s.pair.add(r.read, f.read.reverse_complement())
        # Close samples #
        for s in self.samples: s.pair.close()
        1/0

    #-------------------------------------------------------------------------#
    def primer_statistics(self):
        """Print all primer statistics."""
        # Check goofy or regular #
        regular = 0
        goofy   = 0
        for f,r in tqdm(self.pair.parse_primers(self.primers, self.primer_mismatches)):
            if f.fwd_match :
                if r.rev_match:
                    regular += 1
            if r.fwd_match :
                if f.rev_match:
                    goofy += 1
        print "Regular:", regular
        print "Goofy:", goofy
        # Try every primer in every direction on all reads #
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
        pattern = regex.compile("(%s){s<=%i}" % (iupac_pattern(seq), self.primer_mismatches))
        count   = 0
        for r in tqdm(fastq):
            if pattern.search(str(r.seq)): count += 1
        message  = name + " file, " + primer +  " primer, " + sense + ": "
        message += "%.1f%%" % (100.0 * (count / len(fastq)))
        message += " (%s)" % seq
        return message

    #-------------------------------------------------------------------------#
    def guess_barcodes(self, stop_at=30000):
        """Useful when barcodes seem wrong within a multiplexed file."""
        fwd_barcodes = Counter()
        rev_barcodes = Counter()
        barlen       = len(self.samples[0].info['forward_mid'])
        for i, pair in tqdm(enumerate(self.pair.parse_primers(self.primers, self.primer_mismatches))):
            f,r = pair
            if i == stop_at: break
            if not (f.fwd_match and r.rev_match) and not (r.fwd_match and f.rev_match): continue
            if f.fwd_match and r.rev_match:
                fwd_barcodes[str(f.read[f.fwd_start_pos-barlen:f.fwd_start_pos].seq)] += 1
                rev_barcodes[str(r.read[r.rev_start_pos-barlen:r.rev_start_pos].seq)] += 1
            elif r.fwd_match and f.rev_match:
                fwd_barcodes[str(r.read[r.fwd_start_pos-barlen:r.fwd_start_pos].seq)] += 1
                rev_barcodes[str(f.read[f.rev_start_pos-barlen:f.rev_start_pos].seq)] += 1
        print "Forward:\n", '\n'.join(x + ': ' + str(y) for x,y in fwd_barcodes.most_common(30))
        print "Reverse:\n", '\n'.join(x + ': ' + str(y) for x,y in rev_barcodes.most_common(30))
        return fwd_barcodes, rev_barcodes
