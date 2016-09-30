from setuptools import setup, find_packages

__version__ = eval(open('mitty/version.py').read().split('=')[1])
setup(
  name='mitty',
  version=__version__,
  description='Simulator for genomic data',
  author='Seven Bridges Genomics',
  author_email='kaushik.ghose@sbgenomics.com',
  packages=find_packages(include=['mitty*']),
  include_package_data=True,
  entry_points={
    # Command line scripts
    'console_scripts': [
      # Main programs
      'mitty = mitty.cli:cli'
    ],
    # Register the built in plugins
    'mitty.plugins.sfs': ['double_exp = mitty.plugins.site_frequency.double_exp'],
    'mitty.plugins.variants': ['snp = mitty.plugins.variants.snp_plugin',
                               'delete = mitty.plugins.variants.delete_plugin',
                               'uniformdel = mitty.plugins.variants.uniform_deletions',
                               'uniformins = mitty.plugins.variants.uniform_insertions',
                               'insert = mitty.plugins.variants.insert_plugin',
                               #'inversion = mitty.plugins.variants.inversion_plugin',
                               #'low_entropy_insert = mitty.plugins.variants.low_entropy_insert_plugin'
                               ],
    'mitty.plugins.variants.hotspot': ['uniform = mitty.plugins.variants.hotspots.uniformhot'],
    'mitty.plugins.population': ['standard = mitty.plugins.population.standard',
                                 'vn = mitty.plugins.population.vn'],
    'mitty.plugins.reads': ['simple_sequential = mitty.plugins.reads.simple_sequential_plugin',
                            'simple_se = mitty.plugins.reads.simple_se',
                            'simple_illumina = mitty.plugins.reads.simple_illumina_plugin'],
  },
  install_requires=[
    'sphinx',
    'sphinx_rtd_theme',
    'cython',
    'setuptools>=11.0.0',
    'numpy>=1.9.0',
    'docopt>=0.6.2',
    'click>=3.3',
    'pysam==0.8.4',  # 0.9.0 gives a StringIO error. 0.8.3 is confused by some set_tag operations
    'h5py>=2.5.0',
    'matplotlib>=1.3.0',
    'scipy',
    'pandas>=0.18.1',
    'tables>=3.2.2'
  ],
  test_suite='nose.collector',
  tests_require=['nose'],
)