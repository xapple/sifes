# Built-in modules #

# Internal modules #

# First party modules #
from plumbing.graphs    import Graph
from plumbing.autopaths import FilePath

# Third party modules #
import folium

###############################################################################
class MapFigure(Graph):
    """Draw a pretty map. Uses Folium."""

    default_zoom = 13

    def __init__(self, *args, **kwargs):
        self.markers = []
        super(MapFigure, self).__init__(*args, **kwargs)

    def add_marker(self, lng=47, lat=16, icon='map-pin', color='black'):
        self.markers.append(folium.Marker([lng, lat], icon=folium.Icon(color=color, icon=icon)))

    def save_map(self, **kwargs):
        # Prepare #
        self.select_lng_lat()
        self.map = folium.Map(location=(self.lng, self.lat), zoom_start=self.default_zoom)
        for marker in self.markers: marker.add_to(self.map)
        # Code duplication from plumbing.graphs.Graph #
        self.params = {}
        for key in self.default_params:
            if key in kwargs:                          self.params[key] = kwargs[key]
            elif hasattr(self, key):                   self.params[key] = getattr(self, key)
            elif self.default_params[key] is not None: self.params[key] = self.default_params[key]
        if 'path' in self.params:   path = FilePath(self.params['path'])
        elif hasattr(self, 'path'): path = FilePath(self.path)
        else:                       path = FilePath(self.short_name + '.pdf')
        # Save the map #
        self.map.save(path)
        # Free memory #
        del self.map

    def select_lng_lat(self):
        # Case only one marker #
        if len(self.markers) == 1:
            self.lng = self.markers[0].location[0]
            self.lat = self.markers[0].location[1]
        # Case find the middle #
        else:
            pass