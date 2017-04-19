# Built-in modules #
from __future__ import division

# Built-in modules #
import re, multiprocessing

# Internal modules #

# First party modules #
from fasta              import FASTA, FASTQ
from plumbing.autopaths import AutoPaths, DirectoryPath
from plumbing.cache     import property_cached
from plumbing.common    import tail, flatter

# Third party modules #
import sh
from shell_command import shell_call

###############################################################################
class Pandaseq(object):
    """Use PANDAseq to join read pairs together.

    You can adjust parameters like this:
        sifes.joining.pandaseq.Pandaseq.minimum_overlap = 40
        sifes.joining.pandaseq.Pandaseq.kmer_table_size = 4
    """

    # Attributes #
    short_name = 'pandaseq'
    long_name  = 'PANDAseq Version 2.11'
    executable = 'pandaseq'
    url        = 'https://github.com/neufeld/pandaseq/'
    doc        = 'http://neufeldserver.uwaterloo.ca/%7Eapmasell/pandaseq_man1.html'
    license    = 'GPLv3'
    dependencies = []

    # Parameters #
    minimum_overlap = 2
    kmer_table_size = 3

    all_paths = """
    /assembled.fasta
    /unassembled.fastq
    /stderr.out
    """

    def __repr__(self): return '<%s object on %s>' % (self.__class__.__name__, self.pair)
    def __nonzero__(self): return bool(self.p.assembled)

    def __init__(self, pair, result_dir, sample_name):
        # Save attributes #
        self.pair        = pair
        self.result_dir  = result_dir
        self.sample_name = sample_name
        # Auto paths #
        self.base_dir = DirectoryPath(self.result_dir + self.short_name)
        self.p = AutoPaths(self.base_dir, self.all_paths)

    def run(self, cpus=1, verbose=False):
        """-N: Eliminate uncalled bases
           -F: Keep qual scores though they are conceptually different.
           -T: The number of threads to spawn.
           -u: Write sequences for which the optimal alignment cannot be computed
               to a file as concatenated pairs.
           -o: Miniumum overlap.
           -k: K-mer table size."""
        # Message #
        if verbose: print "Joining sample '%s'" % self.sample_name
        # Number of cores #
        if cpus is None: cpus = min(multiprocessing.cpu_count(), 32)
        # Check version #
        assert "pandaseq 2.11" in sh.pandaseq('-v', _ok_code=[0,1]).stderr
        # Command #
        self.command = 'pandaseq -T 1 -f %s -r %s -u %s -o %s -k %s -F 1> %s 2> %s'
        self.command = self.command % (self.pair.fwd, self.pair.rev,
                                       self.p.unassembled,
                                       self.minimum_overlap,
                                       self.kmer_table_size,
                                       self.p.assembled,
                                       self.p.stderr)
        # Call it #
        shell_call(self.command)
        # Check #
        assert self.p.assembled
        # Check counts #
        count_lowqual     = self.results.stats['lowqual']
        count_noalign     = self.results.stats['noalign']
        count_assembled   = len(self.results.assembled)
        count_unassembled = len(self.results.unassembled)
        count_input       = len(self.pair)
        # Assertions #
        assert count_assembled + count_unassembled + count_lowqual == count_input
        assert count_noalign == count_unassembled

    def clean(self):
        self.p.unassembled.remove()
        self.p.stderr.remove()
        self.p.assembled.remove()

    @property_cached
    def results(self):
        results = PandaseqResults(self)
        message = "You can't access results from PANDAseq before running the algorithm."
        if not results: raise Exception(message)
        return results

###############################################################################
class PandaseqResults(object):

    def __nonzero__(self): return bool(self.p.assembled)
    def __len__(self):     return len(self.p.assembled)

    def __init__(self, pandaseq):
        # Attributes #
        self.parent      = pandaseq
        self.pandaseq    = pandaseq
        self.p           = pandaseq.p
        self.pair        = pandaseq.pair
        self.assembled   = FASTQ(self.p.assembled)
        self.unassembled = FASTA(self.p.unassembled)

    @property_cached
    def stats(self):
        # The raw text #
        result = {'raw': tail(self.p.out)}
        # Check errors #
        if "pandaseq: error" in result['raw']: raise Exception("Pandaseq did not run properly.")
        if result['raw'].startswith("ERR\t"):  raise Exception("Pandaseq did not run properly.")
        # Parse it #
        result['noalign'] = int(re.findall('STAT\tNOALGN\t(.+)$',     result['raw'], re.M)[0])
        result['lowqual'] = int(re.findall('STAT\tLOWQ\t(.+)$',       result['raw'], re.M)[0])
        result['distrib'] =     re.findall('STAT\tOVERLAPS\t(.+)$',   result['raw'], re.M)
        result['distrib'] = map(int, result['distrib'][0].split())
        result['lengths'] = flatter([[i+1]*v for i,v in enumerate(result['distrib'])])
        result['loss']    = 100 * sum(result['distrib'][100:]) / sum(result['distrib'])
        # Return #
        return result

    @property
    def unassembled_count(self):
        return self.stats['lowqual'] + self.stats['noalign']

    @property
    def unassembled_percent(self):
        percent = (self.unassembled_count / len(self.parent.pair)) * 100.0
        return "%.1f%%" % percent