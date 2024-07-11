#!/usr/bin/env python3

from setuptools import setup

setup(name = 'OilClassification',
      version = '1',
      packages = ['OilClassification',],
      install_requires = ['numpy', 'matplotlib', 'PyYAML',
                          'spectral', 'gdal', 'netCDF4',
                          'jupyter', 'jupyterlab', 'scipy',
                          'scikit-image'])
