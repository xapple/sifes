# Built-in modules #

# Internal modules #
from sifes.location.map_figure import MapFigure

# First party modules #

# Third party modules #

# Constants #
__all__ = ['ClusterLocationMap']

###############################################################################
class ClusterLocationMap(MapFigure):
    """Map of several sample locations."""
    short_name = 'cluster_location_map'

    def __init__(self, group_name, samples, *args, **kwargs):
        self.group_name = group_name
        self.samples    = samples
        super(ClusterLocationMap, self).__init__(*args, **kwargs)

    def plot(self, **kwargs):
        # Colors based on custom attributes #
        colors = list(set(s.attribute for s in self.samples))
        # Loop #
        for s in self.samples:
            self.add_marker(s.latitude, s.longitude, color=colors.index(s.attribute))
        # Do it #
        self.save_map(**kwargs)