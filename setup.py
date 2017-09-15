from setuptools import setup, find_packages

with open('long_description.rst') as f:
  ld = f.read()

__version__ = eval(open('mitty/version.py').read().split('=')[1])
setup(
  name='mitty',
  version=__version__,
  description='Simulator for genomic data',
  long_description=ld,
  author='Seven Bridges Genomics',
  author_email='kaushik.ghose@sbgdinc.com',
  url='https://github.com/sbg/Mitty',
  download_url='https://github.com/sbg/Mitty/archive/development.zip',
  keywords=['simulator', 'genomics', 'ngs', 'read mapper', 'aligner', 'variant caller'],
  classifiers=[
    'Development Status :: 5 - Production/Stable',
    'Intended Audience :: Science/Research',
    'Topic :: Scientific/Engineering :: Bio-Informatics',
    'License :: OSI Approved :: Apache Software License',
    'Programming Language :: Python :: 3.4',
  ],
  python_requires='>=3.4',
  packages=find_packages(include=['mitty*']),
  include_package_data=True,
  entry_points={'console_scripts': ['mitty = mitty.cli:cli']},
  install_requires=[
    'cython',
    'setuptools>=24.3.0',
    'numpy>=1.11.0',
    'click>=3.3',
    'pysam==0.10.0',
    'matplotlib>=1.3.0',
    'scipy',
    'pandas',
    'cytoolz',
    'nose'
  ],
  test_suite='nose.collector',
  tests_require=['nose'],
)