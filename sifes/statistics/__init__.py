# Built-in modules #

# Internal modules #
from plumbing.autopaths import AutoPaths
from illumitag.clustering.statistics.nmds import NMDS
from illumitag.clustering.statistics.permanova import PERMANOVA
from illumitag.clustering.statistics.betadis import BetaDispersion
from illumitag.clustering.statistics.unifrac import Unifrac

# Third party modules #

###############################################################################
class StatsOnTaxonomy(object):
    """All the statistics to perform on an OTU table that has taxonomy
    associated to it."""

    all_paths = """
    /nmds/
    /permanova/
    /betadis/
    /unifrac/
    """

    def __init__(self, parent):
        # Save parent #
        self.tax, self.parent = parent, parent
        # Paths #
        self.p = AutoPaths(self.parent.p.stats_dir, self.all_paths)
        # R stuff #
        self.nmds = NMDS(self, self.parent.otu_csv)
        self.permanova = PERMANOVA(self)
        self.betadis = BetaDispersion(self)
        # Other #
        self.unifrac = Unifrac(self)

    def run(self):
        self.nmds.run()
        self.permanova.run()
        self.betadis.run()

###############################################################################
class StatsOnComposition(object):

    all_paths = """
    /nmds/
    /permanova/
    /betadis/
    """

    def __init__(self, parent):
        # Save parent #
        self.composition, self.parent = parent, parent
        # Paths #
        self.p = AutoPaths(self.parent.p.stats_dir, self.all_paths)
        # Children #
        self.nmds = NMDS(self, self.parent.taxa_csv)

    def run(self):
        self.nmds.run()
