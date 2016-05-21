# Built-in modules #

# Internal modules #
import os

# First party modules #

# Third party modules #
from dropbox.client import DropboxClient

###############################################################################
class DropBoxUpload(object):
    """http://dropbox-sdk-python.readthedocs.io/en/master/"""

    def __init__(self, input_dir, output_dir='/Micans V6 analysis delivery'):
        self.dbx = DropboxClient('_WV-ZSiONaAAAAAAAAAAC6xOZMQp1hDcU9sHMOreqe-GwHxuFvcKJ7JXadxD08Ch')
        self.input_dir  = input_dir
        self.output_dir = output_dir

    def run(self):
        for root, dirs, files in os.walk(self.input_dir):
            for filename in files:
                # construct the full local path
                local_path    = os.path.join(root, filename)
                # construct the full Dropbox path
                relative_path = os.path.relpath(local_path, self.input_dir)
                dropbox_path  = os.path.join(self.output_dir, relative_path)
                # upload the file
                with open(local_path, 'rb') as f:
                    self.dbx.put_file(dropbox_path, f)
