# Internal modules #

# Third party modules #
from skbio.math.stats.ordination import RDA

################################################################################
class RedundancyAnalysis(object):
    """See https://sites.google.com/site/mb3gustame/constrained-analyses/rda
    for more information about RDA."""

    def __init__(self, cluster):
        # Attributes #
        self.cluster, self.parent = cluster, cluster

    def run(self):
        ordination_result = CCA(Y, X, scale_Y=True)
