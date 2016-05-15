# Built-in modules #
from collections import Counter

# Internal modules #

# First party modules #
from plumbing.autopaths import AutoPaths
from plumbing.cache     import property_cached

###############################################################################
class SeqFilter(object):
    """- Filter primers:
          * Check that the primers are found where they should be found.
          * Check that the primers have the sequence they should have.
       - Filter based on PHRED scores.
       - Filter based on lengths."""

    # Attributes #
    short_name = 'seq_filter'

    # Parameters #
    primer_mismatches = 2

    all_paths = """
    /primers.fastq
    /phred.fastq
    /length.fastq
    """

    def __repr__(self): return '<%s object on %s>' % (self.__class__.__name__, self.input)
    def __nonzero__(self): return bool(self.p.x)

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

    def run(self):
        # Primers #
        pass
        # Quality #

        # Length #


    # Primers #
    def primer_filter(self):
    #         for r in self.flipped_reads.parse_primers(mismatches=self.primer_mismatches):
    #     if r.fwd_start_pos is not None and r.rev_start_pos is not None:
    #         if r.fwd_start_pos < 70 and r.rev_start_pos > -70: self.good_primers.add_seq(r.read)
        def good_primer_iterator(reads):
            for read in reads:
                averaged = moving_average(read.letter_annotations["phred_quality"], self.qual_windowsize)
                if any([value < self.qual_threshold for value in averaged]): continue
                yield read
        self.qual_filtered.write(good_qual_iterator(self.n_filtered))

    # Quality #
    def phred_filter(self):
        def good_phred_iterator(reads):
            for read in reads:
                averaged = moving_average(read.letter_annotations["phred_quality"], self.qual_windowsize)
                if any([value < self.qual_threshold for value in averaged]): continue
                yield read
        self.qual_filtered.write(good_qual_iterator(self.n_filtered))

    # Length #
    def len_filter(self):
        def good_len_iterator(reads):
            for read in reads:
                if len(read) < self.min_length: continue
                yield read
        self.len_filtered.write(good_len_iterator(self.qual_filtered))

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

    def __nonzero__(self): return bool(self.x)

    def __init__(self, x):
        self.x    = x

###############################################################################
