# Built-in modules #

# Internal modules #

# First party modules #

# Third party modules #
from tqdm import tqdm
import ckanapi, requests, pandas
import cStringIO as StringIO

###############################################################################
class CkanSamples(object):
    """Takes care of adding samples to a distant CKAN server.
    See example at:

    https://github.com/ckan/example-add-dataset/blob/master/add_example_datasets.py

    Example code:

        server = ckanapi.RemoteCKAN(self.server_address, user_agent='sifes', apikey=self.api_key)
        groups = server.action.group_list(id='data-explorer')
        print groups

    Resources are files, packages are datasets that contain one or more files.
    Groups contain one or more datasets. Organisations own datasets.
    """

    server_address = "http://anaerobes.science"
    api_key        = "fba01b2e-5ead-4d29-a58a-baba01867092"

    def __init__(self, samples):
        # Save attributes #
        self.samples = samples

    def run(self):
        # Main object #
        self.server = ckanapi.RemoteCKAN(self.server_address, user_agent='sifes', apikey=self.api_key)
        # Do it #
        #print "Making datasets"
        #for s in tqdm(self.samples): self.make_dataset(s)
        print "Making groups"
        self.make_groups()
        #print "Uploading all data"
        #for s in tqdm(self.samples): self.upload_resources(s)

    def make_dataset(self, s):
        # Identifier #
        short_name = s.short_name
        # Tag #
        tags = []
        if int(s.short_name[-1]) == 1: tags = [{'name': 'top'}]
        if int(s.short_name[-1]) == 8: tags = [{'name': 'bottom'}]
        # Extras #
        extras = s.info.copy()
        extras.pop('contacts', None)
        extras.pop('used', None)
        extras.pop('prefix', None)
        extras.pop('directory', None)
        extras.pop('suffix', None)
        extras.pop('organization', None)
        extras.pop('project_short_name', None)
        extras.pop('project_long_name', None)
        extras.pop('sample_long_name', None)
        extras.pop('sample_long_name', None)
        extras['primers'] = extras['primers']['forward']['name'] + ' and ' + extras['primers']['reverse']['name']
        extras.pop('organism', None)
        extras.pop('env_biome', None)
        extras.pop('env_feature', None)
        extras.pop('env_material', None)
        extras.pop('location', None)
        extras = [{'key':k, 'value':v} for k,v in extras.items()]
        # Try to purge #
        try: self.server.action.dataset_purge(id=short_name)
        except ckanapi.errors.NotFound: pass
        # The package #
        self.server.action.package_create(
            name                      = short_name,
            title                     = s.long_name,
            author                    = s.info['contacts']['contact_one']['name'],
            author_email              = s.info['contacts']['contact_one']['email'],
            maintainer                = s.info['contacts']['contact_two']['name'],
            maintainer_email          = s.info['contacts']['contact_two']['email'],
            license_id                = "uk-ogl", # print server.action.license_list()[11]
            #notes                     = None,
            #url                       = None,
            #version                   = None,
            state                     = 'active',
            type                      = 'dna_sequences',
            tags                      = tags,
            extras                    = extras,
            owner_org                 = 'envonautics',
        )

    def make_groups(self):
        # Groups #
        group_names = set([s.short_name[0:2] for s in self.samples])
        # Loop #
        for name in group_names:
            samples = [s for s in self.samples if s.short_name.startswith(name)]
            self.server.action.group_create(
                name   = name,
                title  = "Reactor " + name.upper(),
                packages = [{'id': s.short_name} for s in samples]
            )
        # Return success #
        return True

    def upload_resources(self, s):
        # The fasta file #
        self.server.action.resource_create(
            url           = self.server_address + "/dataset/" + short_name,
            name          = "joined_reads.fasta",
            description   = "The joined reads before filtering as a FASTA",
            package_id    = s.short_name,
            upload        = open(s.joiner.results.assembled),
            resource_type = 'joined_reads',
            format        = 'fasta',
            mimetype      = 'text/plain',
            size          = len(s.joiner.results.assembled),
        )
        # A graph #
        self.server.action.resource_create(
            url           = self.server_address + "/dataset/" + short_name,
            name          = "joined_reads_len_dist.pdf",
            description   = "Length distribution of the joined read as a PDF",
            package_id    = s.short_name,
            upload        = open(s.joiner.results.assembled.graphs.length_dist(), 'rb'),
            resource_type = 'plot',
            format        = 'pdf',
            mimetype      = 'application/pdf',
        )
        # Genera table #
        self.server.action.resource_create(
            url           = self.server_address + "/dataset/" + short_name,
            name          = "top_20_genera.tsv",
            description   = "Top 20 genera found in this sample as a TSV table",
            package_id    = s.short_name,
            upload        = self.make_genera_table(s),
            resource_type = 'otu_counts',
            format        = 'tsv',
            mimetype      = 'text/plain',
        )
        # Chao1 rarefaction #
        self.server.action.resource_create(
            url           = self.server_address + "/dataset/" + short_name,
            name          = "chao1.pdf",
            description   = "Rarefaction curve of alpha diversity (Chao1)",
            package_id    = s.short_name,
            upload        = open(s.graphs.chao1(), 'rb'),
            resource_type = 'plot',
            format        = 'pdf',
            mimetype      = 'application/pdf',
        )
        # Return success #
        return True

    def make_genera_table(self, s):
        row             = s.project.cluster.taxa_table.results.taxa_table_genus.loc[s.short_name]
        row             = row.sort_values(ascending=False)
        frame           = pandas.DataFrame(index=range(len(row)))
        frame['#']      = range(1, len(row)+1)
        frame['Genera'] = row.index
        frame['Reads']  = [r for r in row.values]
        frame = frame[0:20]
        frame = frame.to_csv(sep='\t', encoding='utf-8', float_format='%.5g')
        frame = StringIO.StringIO(frame)
        return frame

    def via_requests(self):
        """Some example code on how to post directly via HTTP."""
        path       = "1"
        extension  = 1
        short_name = 1
        long_name  = 1
        r = requests.post(self.server_address + '/api/action/resource_create',
                          data = {'package_id': short_name,
                                  'name':       long_name,
                                  'format':     extension,
                                  'url':        'upload'},
                          headers = {'Authorization': self.api_key},
                          files   = [('upload', file(path))])
        if r.status_code != 200:
            message = 'Error while creating resource: {0}'.format(r.content)
            raise Exception(message)