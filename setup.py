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
  entry_points={'console_scripts': ['mitty = mitty.cli:cli']},
  install_requires=[
    'sphinx',
    'sphinx_rtd_theme',
    'cython',
    'setuptools>=11.0.0',
    'numpy>=1.9.0',
    'docopt>=0.6.2',
    'click>=3.3',
    'pysam',
    'h5py>=2.5.0',
    'matplotlib>=1.3.0',
    'scipy',
    'pandas>=0.18.1',
    'tables>=3.2.2'
    'nose'
  ],
  test_suite='nose.collector',
  tests_require=['nose'],
)