# Built-in modules #
from collections import Counter

# Internal modules #

# First party modules #
from fasta import FASTA
from plumbing.autopaths import AutoPaths
from plumbing.cache     import property_cached

###############################################################################
class SeqFilter(object):
    """- Filter primers:
          * Check that the primers are found where they should be found.
          * Check that the primers have the sequence they should have.
       - Filter based on N bases.
       - Filter based on lengths.

    You can adjust parameters like this:
        sifes.filtering.seq_filter.SeqFilter.primer_mismatches = 0
        sifes.filtering.seq_filter.SeqFilter.primer_max_dist   = 25
        sifes.filtering.seq_filter.SeqFilter.min_read_length   = 100
        sifes.filtering.seq_filter.SeqFilter.max_read_length   = 160
    """

    # Attributes #
    short_name = 'seq_filter'

    # Parameters #
    primer_mismatches = 2
    primer_max_dist   = 70
    min_read_length   = -1
    max_read_length   = -1

    all_paths = """
    /primers.fasta
    /n_base.fasta
    /length.fasta
    """

    def __repr__(self): return '<%s object on %s>' % (self.__class__.__name__, self.input)
    def __nonzero__(self): return bool(self.length_fasta)

    def __init__(self, input, result_dir, primers):
        # Save attributes #
        self.input      = input
        self.result_dir = result_dir
        self.primers    = primers
        # Auto paths #
        self.base_dir = self.result_dir + self.short_name + '/'
        self.p = AutoPaths(self.base_dir, self.all_paths)
        # Parsing the primers, returns a GenWithLength #
        self.parse_primers = lambda: self.input.parse_primers(self.primers, self.primer_mismatches)
        # The different files #
        self.primers_fasta = FASTA(self.p.primers)
        self.n_base_fasta  = FASTA(self.p.n_base)
        self.length_fasta  = FASTA(self.p.length)

    def run(self):
        # Primers #
        self.primer_filter()
        # N bases #
        self.n_base_filter()
        # Length #
        self.len_filter()

    # Primers #
    def primer_filter(self):
        """Will take only reads that have both primers and will trim the primers"""
        def good_primer_iterator(reads):
            for r in reads:
                if r.fwd_start_pos is None or r.rev_start_pos is None: continue
                if r.fwd_start_pos > self.primer_max_dist:             continue
                if r.rev_start_pos < -self.primer_max_dist:            continue
                yield r[r.fwd_end_pos:r.rev_end_pos]
        self.primers_fasta.write(good_primer_iterator(self.parse_primers()))

    # N base #
    def n_base_filter(self):
        def good_n_base_iterator(reads):
            for r in reads:
                if 'N' in r: continue
                yield r
        self.n_base_fasta.write(good_n_base_iterator(self.primers_fasta))

    # Length #
    def len_filter(self):
        def good_len_iterator(reads):
            for r in reads:
                if self.min_length > 0:
                    if len(r) < self.min_length: continue
                if self.max_length > 0:
                    if len(r) > self.max_length: continue
                yield r
        self.length_fasta.write(good_len_iterator(self.n_base_fasta))

    #-------------------------------------------------------------------------#
    @property_cached
    def primer_positions(self):
        # Count positions #
        all_fwd_pos, all_rev_pos = Counter(), Counter()
        for r in self.parse_primers():
            if r.fwd_start_pos is not None: all_fwd_pos.update((r.fwd_start_pos,))
            if r.rev_start_pos is not None: all_rev_pos.update((r.rev_start_pos,))
        # Return results #
        return all_fwd_pos, all_rev_pos

    @property_cached
    def results(self):
        results = SeqFilterResults(self)
        message = "You can't access results from filtering before running the algorithm."
        if not results: raise Exception(message)
        return results

###############################################################################
class SeqFilterResults(object):

    def __nonzero__(self): return bool(self.length_fasta)

    def __init__(self, parent):
        self.parent        = parent
        self.primers_fasta = parent.primers_fasta
        self.n_base_fasta  = parent.n_base_fasta
        self.length_fasta  = parent.length_fasta
        # The final result #
        self.clean         = self.length_fasta

###############################################################################
