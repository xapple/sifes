# Built-in modules #
import os, shutil

# Internal modules #
from illumitag.groups.outcomes   import BarcodeGroup
from illumitag.groups.assemble   import Assembled, Unassembled
from illumitag.groups.samples    import Samples
from illumitag.helper.trimconcat import TrimerAndConcactenater
from illumitag.helper.primers    import TwoPrimers
from illumitag.helper.mothur     import Mothur
from illumitag.helper.sra        import SampleSRA
from illumitag.graphs            import outcome_plots
from illumitag.running.presample_runner   import PresampleRunner
from illumitag.reporting.samples          import SampleReport
from illumitag.clustering.diversity.alpha import AlphaDiversity

# First party modules #
from fasta import FASTA, FASTQ
from fasta import PairedFASTQ
from fasta.fastqc import FastQC
from plumbing.autopaths import AutoPaths, FilePath, DirectoryPath
from plumbing.cache import property_cached
from plumbing.common import load_json_path

# Third party modules #
from shell_command import shell_call

# Constants #
home = os.environ['HOME'] + '/'

###############################################################################
class Presample(BarcodeGroup):
    """A Presample is a clumsy name for a new type of barcoded-sequence files.
    As we updated the lab protocol, sample are not multiplexed with
    our traditional 50 barcodes anymore, but with Illumina specific MIDs.
    The demultiplexing thus happens in their pipeline and we are left with one
    sample per file.
    This object is a bit like a *Pool*, a *BarcodeGroup* and a *Sample*
    all at the same time. In the end it inherits from BarcodeGroup and
    just emulates the behavior of the other objects."""

    all_paths = """
    /info.json
    /uncompressed/fwd.fastq
    /uncompressed/rev.fastq
    /logs/
    /assembled/
    /unassembled/
    /fastqc/fwd/
    /fastqc/rev/
    /graphs/
    /report/report.pdf
    /quality/trimmed.fastq
    /quality/renamed.fastq
    /quality/reads.fasta
    """

    kind = 'presample'

    def __repr__(self): return '<%s object "%s">' % (self.__class__.__name__, self.id_name)
    def __str__(self): return self.id_name
    def __iter__(self): return iter(self.children)
    def __len__(self): return self.count
    def __getitem__(self, key): return self.samples[key]

    @property
    def seq_len(self): return len(self.fwd.first_read)

    def __init__(self, json_path, out_dir):
        # Attributes #
        self.out_dir = out_dir
        self.json_path = FilePath(json_path)
        # Parse #
        self.info = load_json_path(self.json_path)
        # Basic #
        self.account = self.info['uppmax_id']
        self.run_num = self.info['run_num']
        self.run_label = self.info['run_id']
        self.project_short_name = self.info['project']
        self.project_long_name = self.info['project_name']
        self.fwd_name = self.info['forward_reads']
        self.rev_name = self.info['reverse_reads']
        # Own attributes #
        self.num = self.info['sample_num']
        self.label = self.info['sample_id']
        self.short_name = self.info['sample']
        self.long_name = self.info['sample_name']
        self.name = 'run%i_sample%i' % (self.run_num, self.num)
        self.group = self.info.get('group')
        self.id_name = "run%03d-sample%02d" % (self.run_num, self.num)
        self.fwd_mid = self.info['forward_mid']
        self.rev_mid = self.info['reverse_mid']
        self.used = True
        # Check name is ASCII #
        assert all(ord(c) < 128 for c in self.short_name)
        # Pool dummy #
        self.pool, self.parent = self, self
        # Second initialization #
        self.loaded = False

    def load(self):
        """A second __init__ that is delayed and called only if needed"""
        # Automatic paths #
        self.base_dir = DirectoryPath(self.out_dir + self.id_name + '/')
        self.p = AutoPaths(self.base_dir, self.all_paths)
        # Special #
        self.primers = TwoPrimers(self)
        # Samples dummy #
        self.info['samples'] = [{"name":self.short_name, "used":1, "group":self.group,
                                 "dummy":1, "num":self.num, "fwd":"", "rev":""}]
        self.samples = Samples(self)
        self.samples.load()
        self.loaded = True
        # Special submission attributes #
        self.sra = SampleSRA(self)
        # Check we can load it #
        if "hostname" == "uppmax":
            if not os.access('/proj/%s' % self.account, os.R_OK):
                raise Exception("You don't have access to the project %s" % self.account)
        # Get the file paths #
        self.fwd_path = home + "proj/%s/INBOX/%s/%s/%s" % (self.account, self.run_label, self.label, self.fwd_name)
        self.rev_path = home + "proj/%s/INBOX/%s/%s/%s" % (self.account, self.run_label, self.label, self.rev_name)
        self.gzipped = True if self.fwd_path.endswith('gz') else False
        self.fastq = PairedFASTQ(self.fwd_path, self.rev_path, self)
        self.fwd = self.fastq.fwd
        self.rev = self.fastq.rev
        # Check the files exist #
        if not self.fwd.exists: raise Exception("No file at '%s'" % self.fwd)
        if not self.rev.exists: raise Exception("No file at '%s'" % self.rev)
        # Change location of first FastQC, we don't want to modify the INBOX #
        self.fastq.fwd.fastqc = FastQC(self.fastq.fwd, self.p.fastqc_fwd_dir)
        self.fastq.rev.fastqc = FastQC(self.fastq.rev, self.p.fastqc_rev_dir)
        # Barcode length #
        self.bar_len = 0
        # Make an alias to the json #
        self.p.info_json.link_from(self.json_path, safe=True)
        # Assembly files as children #
        self.assembled = Assembled('', self)
        self.unassembled = Unassembled('', self)
        self.children = (self.assembled, self.unassembled)
        self.first = self.assembled
        # Special case, for when the two reads don't join #
        self.trim_and_concat = TrimerAndConcactenater(self)
        # Final #
        self.trimmed = FASTQ(self.p.trimmed)
        self.renamed = FASTQ(self.p.renamed)
        self.fasta = FASTA(self.p.reads_fasta)
        # Graphs #
        self.graphs = [getattr(outcome_plots, cls_name)(self) for cls_name in outcome_plots.__all__]
        # Runner #
        self.runner = PresampleRunner(self)
        # Diversity #
        self.diversity = AlphaDiversity(self)
        # Report #
        self.report = SampleReport(self)
        # Loaded #
        self.loaded = True
        # Return self for convenience #
        return self

    @property_cached
    def counts(self):
        """The OTU counts"""
        taxa_table = self.project.cluster.otus.taxonomy.comp_tips.taxa_table
        row = taxa_table.loc[self.short_name].copy()
        row = row.sort_values(ascending=False)
        return row

    def join(self):
        """Uses pandaseq 2.8 to join the foward and reverse reads together.
        See https://github.com/neufeld/pandaseq"""
        # Special case for new primers that don't join #
        rev_primer_name = self.info['primers']['reverse']['name']
        not_joining_primers = ("1132R", "1000R")
        if rev_primer_name in not_joining_primers:
            print "No overlap special case"
            self.trim_and_concat.run()
            return
        # Special case for primers that highly overlap, pandaseq doesn't work somehow #
        high_overlap_primers = ("806R",)
        if rev_primer_name in high_overlap_primers:
            print "High overlap special case, using mothur"
            mothur = Mothur(self.uncompressed_pair)
            mothur.make_contigs(dummy_scores=True)
            shutil.move(mothur.output, self.assembled.path)
            assert len(self.assembled) > 0
            self.unassembled.create()
            self.assembled.p.out.write("# We joined the reads by mothur not pandaseq!")
            return
        # Default case #
        command = 'pandaseq28 -T 1 -f %s -r %s -u %s -F 1> %s 2> %s'
        command = command % (self.fwd, self.rev, self.unassembled.path, self.assembled.path, self.assembled.p.out)
        shell_call(command) # Because it exits with status 1 https://github.com/neufeld/pandaseq/issues/40
        # Check #
        assert self.assembled.path.exists

    def process(self):
        """Remove the primers now"""
        def no_primers_iterator(reads):
            for r in reads: # reads must come from parse_primers() !
                yield r.read[r.fwd_end_pos:r.rev_end_pos]
        reads = self.assembled.good_primers.len_filtered.parse_primers(mismatches=self.assembled.primer_mismatches)
        self.trimmed.write(no_primers_iterator(reads))
        self.trimmed.rename_with_num(self.name + '_read', self.renamed)
        self.renamed.to_fasta(self.fasta)

    def make_mothur_output(self):
        pass

    def make_qiime_output(self):
        pass

    def make_presample_plots(self):
        for graph in self.graphs: graph.plot()

    @property_cached
    def uncompressed_pair(self):
        """Usefull for a few stupid programs that don't take fastq.gz files such as mothur"""
        result = PairedFASTQ(self.p.uncompressed_fwd_fastq, self.p.uncompressed_rev_fastq)
        if not result.exists:
            self.fwd.ungzip_to(result.fwd)
            self.rev.ungzip_to(result.rev)
        return result