# Built-in modules #

# Internal modules #
import time

# First party modules #
from plumbing.graphs    import Graph
from plumbing.autopaths import FilePath

# Third party modules #
import folium
from pyvirtualdisplay import Display
from selenium import webdriver

###############################################################################
class MapFigure(Graph):
    """Draw a pretty map. Uses Folium."""

    default_zoom = 13
    default_delay = 1
    browser_singleton = []

    def __init__(self, *args, **kwargs):
        super(MapFigure, self).__init__(*args, **kwargs)
        self.markers = []
        self.path = self.path.replace_extension('png')

    def add_marker(self, lng=47, lat=16, icon='map-pin', color='black'):
        self.markers.append(folium.Marker([lng, lat], icon=folium.Icon(color=color, icon=icon)))

    def save_map(self, **kwargs):
        # Prepare #
        self.select_lng_lat()
        self.map = folium.Map(location=(self.lng, self.lat), zoom_start=self.default_zoom)
        for marker in self.markers: marker.add_to(self.map)
        # Save the map in HTML #
        self.map.save(self.path.replace_extension('html'))
        # Convert to picture #
        self.browser.get('file://' + self.path.replace_extension('html'))
        time.sleep(self.default_delay)
        self.browser.save_screenshot(self.path)

    @property
    def browser(self):
        if len(self.browser_singleton) == 0:
            self.browser_singleton.append(Display(visible=0, size=(1280, 1024)))
            self.browser_singleton[0].start()
            self.browser_singleton.append(webdriver.Firefox())
        return self.browser_singleton[1]

    def __del__(self):
        self.browser.quit()
        self.display.stop()

    def select_lng_lat(self):
        # Case only one marker #
        if len(self.markers) == 1:
            self.lng = self.markers[0].location[0]
            self.lat = self.markers[0].location[1]
        # Case find the middle #
        else:
            pass