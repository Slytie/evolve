#!/usr/bin/env python

from setuptools import setup
from setuptools import find_packages

__minimum_numpy_version__ = '1.9.0'

setup_requires = ['numpy>=' + __minimum_numpy_version__]
#'tensorflow>='+__minimum_tensorflow_version__]

setup(name='evolve',
      version='0.0.1',
      description='Evolution',
      author=['Tyler Clark'],
      author_email=['tyler.wp.clark@gmail.com'],
    setup_requires=setup_requires,  
    tests_require=[
        'pytest>=2.8',
    ],
    package_dir = {'':'src'},
    packages=find_packages('src')
     )

