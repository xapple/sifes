from distutils.core import setup

setup(
      name             = 'sifes',
      version          = '2.0.1',
      description      = 'Species Identification From Environmental Sequencing.',
      long_description = open('README.md').read(),
      license          = 'Proprietary software, all rights reserved.',
      url              = 'http://github.com/xapple/sifes/',
      author           = 'Lucas Sinclair',
      author_email     = 'lucas.sinclair@me.com',
      classifiers      = ['Topic :: Scientific/Engineering :: Bio-Informatics'],
      requires         = ['plumbing', 'fasta', 'pymarktex', 'sh', 'biopython', 'matplotlib',
                          'threadpool', 'patsy', 'pandas', 'statsmodels', 'rpy2',
                          'scikit-learn', 'rpy2', 'brewer2mpl', 'regex', 'ftputil',
                          'names', 'shell_command', 'pystache', 'tabulate', 'tqdm',
                          'scikit-bio', 'seqenv'],
      packages         = ['sifes',
                          'sifes.groups',
                          ],
)