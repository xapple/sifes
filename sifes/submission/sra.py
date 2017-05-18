# -*- coding: utf-8 -*-

"""
1) Start with a BioProject creation here:
https://submit.ncbi.nlm.nih.gov/subs/bioproject/
Don't add any BioSamples at this step.

2) Then submit the BioSamples (use the first TSV):
https://submit.ncbi.nlm.nih.gov/subs/biosample/

3) Next, make a new submission here (use the second TSV):
https://submit.ncbi.nlm.nih.gov/subs/sra/

4) You can't finish the submission, you will have an option
to go back to the submission list and create a preload folder.
Use the FTP credentials to upload the data.

5) Finish the submission

See your projects here:
https://submit.ncbi.nlm.nih.gov/subs/
http://trace.ncbi.nlm.nih.gov/Traces/sra_sub/sub.cgi?login=pda

Example usage:

    import sifes
    from sifes.groups.projects import Projects
    desalt = sifes.load("~/deploy/sifes/metadata/json/projects/unige/desalt/")
    foram  = sifes.load("~/deploy/sifes/metadata/json/projects/unige/foram/")
    lump = Projects("brine_disposal", [desalt, foram])
    lump.sra.write_bio_tsv()
    lump.sra.write_sra_tsv()
"""

# Built-in modules #
import os, codecs

# Internal modules #
from plumbing.common    import Password, ascii
from plumbing.autopaths import AutoPaths

# Third party modules #
from ftputil import FTPHost

# Constants #
home = os.environ.get('HOME', '~') + '/'

# Old server #
ftp_server    = "ftp-private.ncbi.nih.gov"
ftp_login     = "sra"
ftp_password  = Password("SRA FTP password for user '%s':" % ftp_login)

# Lists #
header_bio = [
    '*sample_name',
    'description',
    'bioproject_id',
    'sample_title',
    '*organism',
    'host',
    'isolation_source',
    '*collection_date',
    '*depth',
    '*elevation',
    '*env_biome',
    '*env_feature',
    '*env_material',
    '*geo_loc_name',
    '*lat_lon',
    'primers',
    'biological_replicate',
]

default_bio = {
    'bioproject_id'   : "PRJNAXXXXXX",
    'organism'        : "aquatic metagenome",
    'host'            : "not applicable",
    'isolation_source': "water",
    'depth'           : "surface",
    'env_biome'       : "river",  # See https://trace.ddbj.nig.ac.jp/faq/biome-feature-material_e.html
    'env_feature'     : "river",
    'env_material'    : "water",
}

header_sra =  [
    'bioproject_accession',
    'biosample_accession',
    'library_ID',
    'title/short description',
    'library_strategy',
    'library_source',
    'library_selection',
    'library_layout',
    'platform',
    'instrument_model',
    'design_description',
    'reference_genome_assembly',
    'alignment_software',
    'forward_read_length',
    'reverse_read_length',
    'filetype1',
    'filename1',
    'MD5_checksum1',
    'filetype2',
    'filename2',
    'MD5_checksum2',
]

default_sra = {
    'library_strategy'    : "AMPLICON",
    'library_source'      : "METAGENOMIC",
    'library_selection'   : "PCR",
    'library_layout'      : "Paired",
    'platform'            : "ILLUMINA",
    'instrument_model'    : "Illumina MiSeq",
    'forward_read_length' : "250",
    'reverse_read_length' : "250",
    'forward_filetype'    : "fastq",
    'reverse_filetype'    : "fastq",
}

###############################################################################
class SampleSRA(object):
    """Every sample has an instance of this class, which enables us to
    access SRA required parameters for every sample."""

    def __init__(self, sample):
        self.s = sample
        self.base_name = '%s_%s_{}_reads.fastq.gz'
        self.base_name = self.base_name % (self.s.project_short_name, self.s.short_name)

    def upload_to_sra(self,
                      address = "ftp-private.ncbi.nlm.nih.gov",
                      usrname = "subftp",
                      passwrd = "w4pYB9VQ",
                      dirctry = "/uploads/lucas.sinclair@me.com_Pm16RLfi",
                      verbose = True):
        """They have an FTP site where you should drop the files first."""
        # Print #
        if verbose: print self.s.short_name + ' (' + self.s.name + ')'
        # Connect #
        if verbose: print "Connecting..."
        self.ftp = FTPHost(address, usrname, passwrd)
        # Test #
        assert self.s.fastq.fwd.count_bytes > 36
        assert self.s.fastq.rev.count_bytes > 36
        # Change directory #
        if verbose: print "Changing directories..."
        self.ftp.chdir(dirctry)
        # Make directory #
        if default_bio['bioproject_id'] not in self.ftp.listdir('.'):
            if verbose: print "Making directories..."
            self.ftp.mkdir(default_bio['bioproject_id'])
        self.ftp.chdir(default_bio['bioproject_id'])
        # Upload #
        if verbose: print "Uploading forward..."
        self.ftp.upload(self.s.fastq.fwd, self.base_name.format("forward"))
        if verbose: print "Uploading reverse..."
        self.ftp.upload(self.s.fastq.rev, self.base_name.format("reverse"))
        # Return #
        self.ftp.close()

    @property
    def biosample_line(self):
        """Will generate the corresponding BioSample entry (first TSV).
        Example usage:
            make_tsv.write_bio_tsv()
        You can add `from plumbing.common import gps_deg_to_float`
        """
        # sample_name #
        line = [self.s.project_short_name + '_' + self.s.short_name]
        # description #
        line += [self.s.project_long_name + ' ' + self.s.long_name]
        # bioproject_id #
        line += [self.s.info['bioproject']] or [default_bio['bioproject_id']]
        # sample_title #
        line += ["Sample '%s' (%s)" % (self.s.short_name, self.s.long_name)]
        # organism #
        line += [self.s.info.get('organism')] # default_bio["organism"]
        # host #
        line += [default_bio["host"]]
        # isolation_source #
        line += [default_bio["isolation_source"]]
        # collection_date #
        line += [self.s.info['date'].split(' ')[0]]
        # depth and elevation #
        line += [self.s.info['depth']]
        line += ['-' + self.s.info['depth']]
        # env_biome, env_feature, env_material #
        line += [self.s.info['env_biome']]    #[default_bio["env_biome"]]
        line += [self.s.info['env_feature']]  #[default_bio["env_feature"]]
        line += [self.s.info['env_material']] #[default_bio["env_material"]]
        # geo_loc_name
        line += ["%s: %s" % (self.s.info['country'], ascii(self.s.info['location']))]
        # lat_lon
        coords = (float(self.s.info['latitude'][0]),  # Latitude
                  float(self.s.info['longitude'][0])) # Longitude
        line += ['{:7.6f} N {:7.6f} E'.format(*coords)]
        # Custom extra lines #
        line += [self.s.info.get('primers')['name']]
        line += [self.s.info.get('replicate_id') + self.s.short_name[-1]]
        # Return #
        return line

    @property
    def sra_line(self):
        """Will generate the corresponding entry for SRA submission (second TSV).
        Example usage:
            make_tsv.write_sra_tsv()
        """
        # project accession
        bioproj = self.s.info['bioproject']
        if bioproj and bioproj != "PRJNAXXXXXX":  line = [bioproj]
        else:                                     line = [default_bio['bioproject_id']]
        # sample accession
        line += [self.s.info['biosample']]
        # library id
        line += [self.s.short_name]
        # description
        machine = 'Illumina MiSeq' if 'machine' not in self.s.__dict__ else self.s.machine
        desc = "Sample '%s' (code %s) sampled on %s and run on a %s"
        desc += " -- run number %i, pool number %i, barcode number %i."
        line += [desc % (self.s.long_name, self.s.short_name, self.s.info['date'],
                         machine, self.s.pool.run_num, self.s.pool.num, self.s.num)]
        # library_strategy
        line += [self.s.info['library_strategy']]  # [default_sra['library_strategy']]
        line += [self.s.info['library_source']]    # [default_sra['library_source']]
        line += [self.s.info['library_selection']] # [default_sra['library_selection']]
        line += [default_sra['library_layout']]    # [self.s.info['library_layout']]
        line += [self.s.info['platform']]          # [default_sra['platform']]
        line += [self.s.info['instrument_model']]  # [default_sra['instrument_model']]
        # design_description
        line += [self.s.pool.info['design_description']]
        # reference_genome_assembly
        line += ['', '']
        # forward_read_length
        line += [str(self.s.info['forward_read_length']), str(self.s.info['reverse_read_length'])]
        # forward
        line += [default_sra['forward_filetype'], self.base_name.format("forward")]
        line += [self.s.fastq.fwd.md5]
        # reverse
        line += [default_sra['reverse_filetype'], self.base_name.format("reverse")]
        line += [self.s.fastq.rev.md5]
        # return
        return line

###############################################################################
class LumpSRA(object):
    """A class to generate the spreadsheets required by the NCBI/SRA
    for raw data submission."""

    all_paths = """
    /biosample_creation.tsv
    /sra_submission.tsv
    """

    def __init__(self, lump):
        """You give a lump as input."""
        # Base parameters #
        self.lump    = lump
        self.samples = lump.samples
        # Auto paths #
        self.base_dir = self.lump.p.sra_dir
        self.p        = AutoPaths(self.base_dir, self.all_paths)

    def write_bio_tsv(self, path=None):
        """Will write the TSV required by the NCBI for the creation of 'BioSample' objects
        (first TSV)."""
        header  = '\t'.join(header_bio) + '\n'
        content = '\n'.join('\t'.join(map(str, s.sra.biosample_line)) for s in self.samples)
        if path is None: path = self.p.biosample
        with codecs.open(path, 'w', encoding='utf-8') as handle:
            handle.write(header+content)
        # Message #
        print "Wrote TSV at '%s'" % path

    def write_sra_tsv(self, path=None):
        """Will write the appropriate TSV for the SRA submission in the cluster directory
        (second TSV). Sometimes you need to set the encoding to `windows-1252`."""
        header  = '\t'.join(header_sra) + '\n'
        content = '\n'.join('\t'.join(map(str,s.sra.sra_line)) for s in self.samples)
        if path is None: path = self.p.sra
        with codecs.open(path, 'w', encoding='utf-8') as handle:
            handle.write(header+content)
        # Message #
        print "Wrote TSV at '%s'" % path
