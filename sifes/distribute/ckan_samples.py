# Built-in modules #

# Internal modules #

# First party modules #

# Third party modules #
import ckanapi

###############################################################################
class CkanSamples(object):
    """Takes care of adding samples to a distant CKAN server."""

    server_address = "http://anaerobes.science"

    def __init__(self, samples):
        # Save attributes #
        self.samples = samples

    def run(self):
        pass