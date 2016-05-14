# Built-in modules #
import os, socket

# Internal modules #
import sifes

# First party modules #
from pymarktex import Template
from plumbing.common import pretty_now

# Constants #
ssh_header = "ssh://" + os.environ.get("FILESYSTEM_HOSTNAME", socket.getfqdn())

###############################################################################
class ReportTemplate(Template):
    """Things that are common to most reports in sifes."""

    # Process info #
    def project_url(self):       return sifes.url
    def project_version(self):   return sifes.__version__
    def now(self):               return pretty_now()
    def git(self):
        if not sifes.git_repo: return False
        return {'git_hash'  : sifes.git_repo.hash,
                'git_tag'   : sifes.git_repo.tag,
                'git_branch': sifes.git_repo.branch}
