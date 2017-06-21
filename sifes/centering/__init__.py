# Built-in modules #
import os
from collections import defaultdict

# First party modules #
from plumbing.autopaths import FilePath

# Third party modules #
import sh

###############################################################################
class UClusterFile(FilePath):
    """A special format outputted by UPARSE
    An example line:
    H       1474    422     97.6    +       0       0       422M    run4_pool1_sample1_read2        OTU-1474

    The corresponding legend:
    # 1=Type, 2=ClusterNr, 3=SeqLength or ClusterSize, 4=PctId, 5=Strand, 6=QueryStart, 7=SeedStart, 8=Alignment, 9=QueryLabel, 10=TargetLabel
    """

    def __len__(self):     return self.count
    def __repr__(self):    return '<%s object on "%s">' % (self.__class__.__name__, self.path)
    def __nonzero__(self): return os.path.getsize(self.path) != 0

    @property
    def count(self):
        return int(sh.wc('-l', self.path).stdout.split()[0])

    @property
    def count_mapped(self):
        """Count only the reads that mapped to a center"""
        return sum([1 for line in self if line.split()[-1] != '*'])

    @property
    def otu_sample_counts(self):
        # Put results in a dict of dicts #
        result = defaultdict(lambda: defaultdict(int))
        # Loop #
        for line in self:
            # Parse the line (the query is the name of the read) #
            kind, num, length, similarity, strand, start, seed, alignment, query, target = line.split()
            # Skip no hits (the target is the OTU name) #
            if target == '*': continue
            # Remove read number to get sample name #
            sample_name = query.split(':')[0]
            # Count #
            result[target][sample_name] += 1
        # Return #
        return result
