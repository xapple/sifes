# Internal modules #
import os, shutil, glob

# First party modules #
from fasta               import FASTA
from plumbing.autopaths  import AutoPaths, FilePath
from plumbing.cache      import property_cached
from plumbing.csv_tables import CSVTable

# Third party modules #
import sh, pandas

# Constants #
home = os.environ.get('HOME', '~') + '/'
reference = FASTA(home + 'glob/16s/silva/v111/rep_set_aligned/97_Silva_111_rep_set.fasta')

################################################################################
class UniFrac(object):
    """A class to compute the UniFrac algorithm producing a distance matrix
    from a bunch of different samples and their reads.

    Step 1. Make an alignment of all the OTU centers against a reference.
    One can use:
        * clustalo
        * PyNAST
        * mothur      <- fastest
        * SINA

    Step 2. From the alignment produced make a tree.
    One can use:
        * RAxML
        * Fasttree    <- fastest

    Step 3. Feed the tree and the OTU table into a UniFrac algorithm
    One can use:
        * Pycogent
        * scikit-bio  <- fastest

    Step 4. Return distance matrix as a TSV.

    Information:
    - http://pycogent.org/examples/unifrac.html
    - http://telliott99.blogspot.se/2010/02/unifrac-analysis-introduction.html
    """

    short_name = "uni_frac"

    all_paths = """
    /clustalo/centers.align
    /pynast/centers.align
    /pynast/log.txt
    /pynast/fail.fasta
    /mothur/centers.align
    /mothur/report.txt
    /mothur/flip.accnos
    /raxml/output.tree
    /fasttree/output.tree
    /distances.csv
    """

    def __init__(self, parent, result_dir):
        # Save parent #
        self.cluster, self.parent = parent, parent
        # Auto paths #
        self.result_dir = result_dir
        self.base_dir = self.result_dir + self.short_name + '/'
        self.p = AutoPaths(self.base_dir, self.all_paths)
        # Alignments #
        self.clustalo_aligned = FASTA(self.p.clustalo_align)
        self.pynast_aligned   = FASTA(self.p.pynast_align)
        self.mothur_aligned   = FASTA(self.p.mothur_align)
        # Other files #
        self.raxml_tree    = FilePath(self.p.raxml_tree)
        self.fasttree_tree = FilePath(self.p.fasttree_tree)
        self.distances_csv = CSVTable(self.p.distances_csv)

    def run(self):
        # Step 1 #
        self.align_mothur()
        # Step 2 #
        self.tree_fasttree()
        # Step 3 #
        self.unifrac_pycogent()
        # Step 4 #
        return self.distances

    #------------------------------ Alignments -------------------------------#
    def align_clustalo(self):
        """Step 1 with clustalo"""
        self.centers_aligned.remove()
        sh.clustalo('-i', self.tax.centers, '--profile1', reference, '-o', self.clustalo_aligned, '--threads', 16)

    def align_pynast(self):
        """Step 1 with PyNAST"""
        sh.pynast('--input_fp', self.tax.centers, '--template_fp', reference,
                  '--fasta_out_fp', self.pynast_aligned,
                  '--log_fpa', self.p.pynast_log, '--failure_fp', self.p.pynast_fail,)

    def align_mothur(self):
        """Step 1 with mothur"""
        # Run it #
        sh.mothur("#align.seqs(candidate=%s, template=%s, search=kmer, flip=false, processors=%s);" \
                  % (self.tax.centers, reference, 16))
        # Move things #
        shutil.move(self.tax.centers.prefix_path + '.align',        self.mothur_aligned)
        shutil.move(self.tax.centers.prefix_path + '.align.report', self.p.mothur_report)
        path = self.tax.centers.prefix_path + '.flip.accnos'
        if os.path.exists(path): shutil.move(path, self.p.mothur_accnos)
        # Clean up #
        for p in glob.glob('mothur.*.logfile'): os.remove(p)

    def align_sina(self):
        """Step 1 with SINA: should use the NR 99% silva database and PT-SERVER ?"""
        sh.SINA('-h')

    #-------------------------------- Trees ---------------------------------#
    def tree_raxml(self):
        """Step 2 with RAxML"""
        sh.raxml('-T', 16, '-s', self.mothur_aligned, '-n', self.raxml_tree, '-m', 'LOREM')

    def tree_fasttree(self):
        """Step 2 with FastTree"""
        sh.FastTreeMP('-fastest' ,'-out', self.fasttree_tree, '-nt', self.mothur_aligned)

    #----------------------------- Algorithms ------------------------------#
    def unifrac_pycogent(self):
        """Step 3 with Pycogent"""
        tree_newick = open(self.fasttree_tree, 'r').read()
        from cogent.parse.tree import DndParser
        from cogent.maths.unifrac.fast_tree import UniFracTreeNode
        tree = DndParser(tree_newick, UniFracTreeNode)
        from cogent.maths.unifrac.fast_unifrac import fast_unifrac
        distances = fast_unifrac(tree, self.tax.otu_table.to_dict())
        # Make a dataframe #
        names = distances['distance_matrix'][1]
        df = pandas.DataFrame(distances['distance_matrix'][0], index=names, columns=names)
        df.to_csv(self.distances_csv, sep='\t', float_format='%.5g')

    def unifrac_skbio(self):
        # The OTU table #
        data = 0
        # The name of the samples (columns) #
        ids = 0
        # The name of the OTUs (rows) #
        otu_ids = ['OTU1', 'OTU2', 'OTU3', 'OTU4', 'OTU5', 'OTU6', 'OTU7']
        # Load the tree #
        from skbio import TreeNode
        tree = TreeNode.read(StringIO())
        # Do the algorithm #
        from skbio.diversity import beta_diversity
        matrix = beta_diversity("weighted_unifrac", data, ids, tree=tree, otu_ids=otu_ids)
        # Save the matrix on disk #
        pass

    #-------------------------------- Result ---------------------------------#
    @property_cached
    def distances(self):
        return pandas.io.parsers.read_csv(self.distances_csv, sep='\t', index_col=0)
