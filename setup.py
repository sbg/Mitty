from setuptools import setup, find_packages

__version__ = eval(open('mitty/version.py').read().split('=')[1])
setup(
  name='mitty',
  version=__version__,
  description='Simulator for genomic data',
  author='Seven Bridges Genomics',
  author_email='kaushik.ghose@sbgenomics.com',
  python_requires='>=3.4',
  packages=find_packages(include=['mitty*']),
  include_package_data=True,
  entry_points={'console_scripts': ['mitty = mitty.cli:cli']},
  install_requires=[
    'cython',
    'setuptools>=24.3.0',
    'numpy>=1.11.0',
    'click>=3.3',
    'pysam>=0.10.0',
    'matplotlib>=1.3.0',
    'scipy',
    'nose'
  ],
  test_suite='nose.collector',
  tests_require=['nose'],
)