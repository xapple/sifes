# Built-in modules #
import os, subprocess

# Internal modules #
import sifes

# First party modules #

# Third party modules #
import sh

# Constants #
home = os.environ.get('HOME', '~') + '/'

###############################################################################
class DropBoxRclone(object):
    """Expects that your dropbox is configured in rclone with the name "prod"
    Check with `$ rclone listremotes -l`."""

    def __init__(self, input_dir, output_dir):
        self.input_dir  = input_dir.rstrip('/')
        self.output_dir = output_dir.rstrip('/')

    @property
    def command(self):
        return 'sync', "'" + self.input_dir + "'", "'" + 'prod:' + self.output_dir + "'"

    def run(self):
        """Just rclone it"""
        self.stdout = subprocess.check_call('rclone --copy-links ' + ' '.join(self.command), shell=True)

###############################################################################
class DropBoxSync(object):
    """Expects that your dropbox is mounted at ~/Dropbox and then uses rsync locally.
    Don't forget to start the daemon:
        $ dropbox.py start
    """

    dbx_mount_point = sifes.home + "Dropbox/"

    def __init__(self, input_dir, output_dir):
        self.input_dir  = input_dir
        self.output_dir = output_dir

    def run(self):
        """Just rsync it and let the other script do the work."""
        assert "already running" in sh.Command("dropbox.py")("status")
        print sh.rsync('-a', self.input_dir, home + 'Dropbox/' + self.output_dir)

###############################################################################
class DropBoxUpload(object):
    """Uses this: http://dropbox-sdk-python.readthedocs.io/en/master/"""

    def __init__(self, input_dir, output_dir='/Micans V6 analysis delivery'):
        from dropbox.client import DropboxClient
        self.dbx = DropboxClient(self.dropbox_token)
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
                    self.dbx.put_file(dropbox_path, f, overwrite=True)

