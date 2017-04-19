# Different databases #
from seqsearch.databases.silva_mothur import silva_mothur
from seqsearch.databases.foraminifera import foraminifera
from seqsearch.databases.pr_two       import pr_two

# First party modules #
from plumbing.autopaths import AutoPaths

###############################################################################
class Classify(object):

    default_database = 'silva_mothur'

    def __init__(self, centers, result_dir, database=None):
        # Attributes #
        self.centers    = centers
        self.result_dir = result_dir
        # Default or user specified #
        if database is None: self.database = self.default_database
        else:                self.database = database
        # The database to use #
        if   self.database == 'silva_mothur': self.database = silva_mothur
        elif self.database == 'foraminifera': self.database = foraminifera
        elif self.database == 'pr_two':       self.database = pr_two
        else: raise NotImplemented()
        # Auto paths #
        self.base_dir = self.result_dir + self.short_name + '/' + self.database.short_name + '/'
        self.p = AutoPaths(self.base_dir, self.all_paths)

    @property
    def message(self):
        msg = "Classifying file '%s' against '%s' with '%s'"
        return msg % (self.centers, self.database.short_name, self.short_name)