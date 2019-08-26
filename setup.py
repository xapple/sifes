from setuptools import setup, find_packages

setup(
    name             = 'sifes',
    version          = '2.0.1',
    description      = 'sifes can analyze DNA data and performs species identification from environmental sequencing.',
    long_description = open('README.md').read(),
    license          = 'Proprietary software, all rights reserved.',
    url              = 'http://github.com/xapple/sifes/',
    author           = 'Lucas Sinclair',
    author_email     = 'lucas.sinclair@me.com',
    packages         = find_packages(),
    classifiers      = ['Topic :: Scientific/Engineering :: Bio-Informatics'],
    requires         = ['plumbing', 'fasta', 'pymarktex', 'sh', 'biopython', 'matplotlib',
                        'threadpool', 'patsy', 'pandas', 'statsmodels', 'rpy2',
                        'scikit-learn', 'rpy2', 'brewer2mpl', 'regex', 'ftputil',
                        'names', 'shell_command', 'pystache', 'tabulate', 'tqdm',
                        'scikit-bio', 'seqenv'],
)