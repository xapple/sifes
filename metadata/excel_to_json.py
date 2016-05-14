#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""
A script that reads an excel file and produces our correctly formated JSON files.

Written by Lucas Sinclair.
Kopimi.

You can use this script from the shell like this:
$ excel_to_json.py ~/repos/illumitag/scripts/runs/run012/run012.xlsx
"""

# Modules #
import sys, os, pandas, codecs, numpy, datetime

# Constants #
home = os.environ['HOME'] + '/'

###############################################################################
template = u"""{
    "contacts": {
        "%(contact_one_function)s": {
            "name": "%(contact_one_name)s",
            "email": "%(contact_one_email)s"
        },
        "%(contact_two_function)s": {
            "name": "%(contact_two_name)s",
            "email": "%(contact_two_email)s"
        }
    },

    "project":       "%(project_short_name)s",
    "project_name":  "%(project_long_name)s",
    "sample":        "%(sample_short_name)s",
    "sample_name":   "%(sample_long_name)s",

    "uppmax_id":     "%(uppmax_proj)s",
    "run_num":       %(run_num)s,
    "run_id":        "%(run_name)s",
    "sample_num":    %(sample_num)s,
    "sample_id":     "%(sample_id)s",
    "forward_reads": "%(fwd_filename)s",
    "reverse_reads": "%(rev_filename)s",
    "forward_md5":   "%(fwd_md5)s",
    "reverse_md5":   "%(rev_md5)s",

    "forward_read_count":   %(fwd_count)s,
    "reverse_read_count":   %(rev_count)s,

    "forward_num":   %(forward_num)s,
    "forward_mid":   "%(forward_mid)s",
    "reverse_num":   %(reverse_num)s,
    "reverse_mid":   "%(reverse_mid)s",
    "barcode_num":   %(barcode_num)s,
    "dna_after":     [%(dna_after)s, "ng/µl"],

    "primers": {
            "name":    "%(primers_name)s",
            "sense":   "5' to 3'",
            "forward": {"name": "%(fwd_primer_name)s", "sequence": "%(fwd_primer_seq)s"},
            "reverse": {"name": "%(rev_primer_name)s", "sequence": "%(rev_primer_seq)s"}
    },

    "library_strategy":     "%(library_strategy)s",
    "library_source":       "%(library_source)s",
    "library_selection":    "%(library_selection)s",
    "library_layout":       "%(library_layout)s",
    "platform":             "%(platform)s",
    "instrument_model":     "%(instrument_model)s",
    "instrument_software":  "%(instrument_software)s",
    "forward_read_length":  %(forward_read_length)s,
    "reverse_read_length":  %(reverse_read_length)s,

    "organism":      "%(organism)s",
    "env_biome":     "%(env_biome)s",
    "env_feature":   "%(env_feature)s",
    "env_material":  "%(env_material)s",

    "date":         "%(date)s",
    "latitude":     [%(latitude)s, "N"],
    "longitude":    [%(longitude)s, "E"],
    "country":      "%(country)s",
    "location":     "%(location)s",

    "design_description": "%(design_description)s",

    "bioproject":   "%(bioproject)s",
    "biosample":    "%(biosample)s"
}"""

###############################################################################
correspondence = {
    u'Uppmax project':                       'uppmax_proj',
    u'Run #':                                'run_num',
    u'Run name':                             'run_name',

    u'Sample #':                             'sample_num',
    u'Sample name':                          'sample_id',
    u'Forward index #':                      'forward_num',
    u'Forward index sequence':               'forward_mid',
    u'Reverse index #':                      'reverse_num',
    u'Reverse index sequence':               'reverse_mid',
    #u'Reverse index (verbatim label)':       None,
    u'Barcode ref.':                         'barcode_num',
    u'DNA con. [ng/µl]':                     'dna_after',
    u'PhiX spiking':                         'lorem',
    u'Forward filename':                     'fwd_filename',
    u'Reverse filename':                     'rev_filename',
    u'Forward reads count':                  'fwd_count',
    u'Reverse reads count':                  'rev_count',
    u'Forward MD5 checksum':                 'fwd_md5',
    u'Reverse MD5 checksum':                 'rev_md5',

    u'Contact 1 function':                   'contact_one_function',
    u'Contact 1 name':                       'contact_one_name',
    u'Contact 1 email':                      'contact_one_email',
    u'Contact 2 function':                   'contact_two_function',
    u'Contact 2 name':                       'contact_two_name',
    u'Contact 2 email':                      'contact_two_email',

    u'Project short name (no spaces and only ascii)': 'project_short_name',
    u'Project long name (free text)':                 'project_long_name',
    u'Sample short name (no spaces and only ascii)':  'sample_short_name',
    u'Sample long name (free text)':                  'sample_long_name',

    u'Primers name':                         'primers_name',
    u'Forward primer name':                  'fwd_primer_name',
    u'Forward primer':                       'fwd_primer_seq',
    u'Reverse primer name':                  'rev_primer_name',
    u'Reverse primer':                       'rev_primer_seq',
    u'Library strategy':                     'library_strategy',
    u'Library source':                       'library_source',
    u'Library selection':                    'library_selection',
    u'Library layout':                       'library_layout',
    u'Platform':                             'platform',
    u'Instrument model':                     'instrument_model',
    u'Instrument software':                  'instrument_software',
    u'Forward read length':                  'forward_read_length',
    u'Reverse read length':                  'reverse_read_length',

    u'Organism':                             'organism',
    u'Environement biome (env_biome)':       'env_biome',
    u'Environement feature (env_feature)':   'env_feature',
    u'Environement material (env_material)': 'env_material',

    u'Sampling Date (YYYY-MM-DD)':           'date',
    u'Latitude (N)':                         'latitude',
    u'Longitude (E)':                        'longitude',
    u'Country':                              'country',
    u'Location (free text)':                 'location',
    u'Design description (free text)':       'design_description',

    u'Bioproject':                           'bioproject',
    u'Biosample':                            'biosample',
}

###############################################################################
extras = {
    u'Depth [m]':                            'lorem',
    u'pH':                                   'lorem',
    u'Temperature [℃]':                      'lorem',
    u'Oxygen (O2) [mg/L]':                   'lorem',
    u'Carbon dioxide (CO2) [µM]':            'lorem',
    u'Methane (CH4) [µM]':                   'lorem',
    u'Iron II [µM]':                         'lorem',
    u'Iron III [µM]':                        'lorem',
    u'TOC [mg/L]':                           'lorem',
    u'SUVA [mg/L*m]':                        'lorem',
    u'total P [µg/L]':                       'lorem',
    u'total N [µg/L]':                       'lorem',
    u'Conductivity [μS/cm]':                 'lorem',
    u'Cell counts [cells/mL]':               'lorem',
    u'Secchi depth [cm]':                    'lorem',
    u'TOP [µg/L]':                           'lorem',
}

###############################################################################
# Get the shell arguments #
if len(sys.argv) < 2: sys.exit(sys.modules[__name__].__doc__)
xlsx_path = sys.argv[1]

# Check that the path is valid #
if not os.path.exists(xlsx_path): raise Exception("No file at %s." % xlsx_path)

# Load data #
xlsx = pandas.ExcelFile(xlsx_path)
df = xlsx.parse(xlsx.sheet_names[0])

# Use the first row (of the dataframe) as the column index and then delete the first row #
df.columns = df.loc[0]
df = df.drop(0)

# Main loop #
for i, row in df.iterrows():
    # Using the column names, make a dict with the ascii codes as keys instead #
    data = dict((correspondence[x], row[x]) for x in row.index if x in correspondence)
    # Skip samples missing short names or project names #
    if data.get('sample_short_name') is numpy.nan: continue
    if data.get('project_short_name') is numpy.nan: continue
    # Sometimes time is introduced in the date #
    if isinstance(data['date'], datetime.datetime): data['date'] = data['date'].date().isoformat()
    # Template the text #
    text = template % data
    text = text.replace('"nan",', 'null,')
    text = text.replace('nan,',   'null,')
    # Figure out the path #
    path = home + "repos/illumitag/json/presamples/run%03d/run%03d-sample%03d.json"
    path = path % (int(data['run_num']), int(data['run_num']), int(data['sample_num']))
    # Create directory if it doens't exist #
    dir_path = os.path.dirname(path)
    if not os.path.exists(dir_path): os.mkdir(dir_path)
    # Write #
    with codecs.open(path, 'w', encoding='utf-8') as handle: handle.write(text)
