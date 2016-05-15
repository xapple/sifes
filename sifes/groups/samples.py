# Built-in modules #

# Internal modules #
import sifes
from sifes.diversity.alpha      import AlphaDiversity
from sifes.report.samples       import SampleReport
from sifes.joining.pandaseq     import Pandaseq
from sifes.joining.qiime_join   import QiimeJoin
from sifes.joining.mothur_join  import MothurJoin
from sifes.filtering.seq_filter import SeqFilter

# First party modules #
from plumbing.autopaths import FilePath, DirectoryPath, AutoPaths
from plumbing.common    import load_json_path
from plumbing.cache     import property_cached
from fasta              import PairedFASTA, PairedFASTQ
from fasta.fastqc       import FastQC
from fasta.primers      import TwoPrimers

# Third party modules #

###############################################################################
class Sample(object):
    """Consists of either two FASTA files or two FASTQ files.
    It's a bunch of paired sequences all coming from the same particular
    IRL lab sample. Might or might not correspond to an Illumina MID."""

    raw_files_must_exist = True

    all_paths = """
    /logs/
    /info.json
    /fastqc/fwd/
    /fastqc/rev/
    /uncompressed/fwd.fastq
    /uncompressed/rev.fastq
    /joined/
    /filtered/
    /report/report.pdf
    """

    def __repr__(self):         return '<%s object "%s">' % (self.__class__.__name__, self.short_name)
    def __str__(self):          return self.short_name
    def __iter__(self):         return iter(self.children)
    def __len__(self):          return self.count
    def __getitem__(self, key): return self.children[key]

    def __init__(self, json_path):
        # Attributes #
        self.json_path = FilePath(json_path)
        # Parse #
        self.info = load_json_path(self.json_path)
        # Own attributes #
        self.num                = self.info.get('sample_num')
        self.short_name         = self.info.get('sample_short_name')
        self.long_name          = self.info.get('sample_long_name')
        self.num                = self.info.get('sample_num')
        # Project #
        self.project_short_name = self.info.get('project_short_name')
        self.project_long_name  = self.info.get('project_long_name')
        # Check the short name is only ASCII #
        assert all(ord(c) < 128 for c in self.short_name)
        assert self.short_name[0] not in "1234567890"
        # Automatic paths #
        self.base_dir = DirectoryPath(sifes.samples_dir + self.short_name + '/')
        self.p        = AutoPaths(self.base_dir, self.all_paths)
        # Make an alias to the json #
        self.p.info_json.link_from(self.json_path, safe=True)
        # Get the directory #
        prefix    = self.info.get('prefix',    '')
        directory = self.info.get('directory', '')
        suffix    = self.info.get('suffix',    '')
        # Get the file paths #
        self.fwd_path = FilePath(prefix + directory + suffix + self.info['fwd_filename'])
        self.rev_path = FilePath(prefix + directory + suffix + self.info['rev_filename'])
        # Is it a FASTA pair or a FASTQ pair ? #
        if "fastq" in self.fwd_path: self.pair = PairedFASTQ(self.fwd_path, self.rev_path)
        else:                        self.pair = PairedFASTA(self.fwd_path, self.rev_path)
        # Check that the files exist #
        if self.raw_files_must_exist:
            self.fwd_path.must_exist()
            self.rev_path.must_exist()
        # For speed let's update the sequence count cache if available #
        if self.info.get('forward_read_count') is not None:
            self.pair.fwd.count = self.info['forward_read_count']
        if self.info.get('reverse_read_count') is not None:
            self.pair.rev.count = self.info['reverse_read_count']
        # Automatic paths #
        self.p = AutoPaths(self.base_dir, self.all_paths)
        # Change location of first FastQC, we don't want to modify the INBOX #
        if self.pair.format == 'fastq':
            self.pair.fwd.fastqc = FastQC(self.pair.fwd, self.p.fastqc_fwd_dir)
            self.pair.rev.fastqc = FastQC(self.pair.rev, self.p.fastqc_rev_dir)
        # The primers #
        self.primer_fwd = self.info.get('primers', {}).get('forward', {}).get('sequence')
        self.primer_rev = self.info.get('primers', {}).get('reverse', {}).get('sequence')
        if self.primer_fwd and self.primer_rev:
            self.primers = TwoPrimers(self.primer_fwd, self.primer_rev)
        # Report in PDF #
        self.report = SampleReport(self)

    @property
    def seq_len(self):
        """The length of the first read."""
        return len(self.pair.fwd.first_read)

    @property_cached
    def uncompressed_pair(self):
        """Useful for a few stupid programs that don't take fastq.gz files such as mothur."""
        result = PairedFASTQ(self.p.uncompressed_fwd_fastq, self.p.uncompressed_rev_fastq)
        if not result.exists:
            self.pair.fwd.ungzip_to(result.fwd)
            self.pair.rev.ungzip_to(result.rev)
        return result

    @property_cached
    def joiner(self):
        """Will put the forward and reverse reads together."""
        return MothurJoin(self.uncompressed_pair, self.p.joined_dir)
        return Pandaseq(self.pair,   self.p.joined_dir)
        return QiimeJoin(self.pair,  self.p.joined_dir)

    @property_cached
    def filter(self):
        """Will filter out unwanted sequences."""
        return SeqFilter(self.joiner.results.assembled, self.p.filtered_dir, self.primers)

    @property_cached
    def diversity(self):
        """For all alpha diversity measures."""
        return AlphaDiversity(self)