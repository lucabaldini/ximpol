#!/usr/bin/env python
#
# Copyright (C) 2015, the ximpol team.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU GengReral Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.


from setuptools import setup, find_packages
import os
import sys
import glob

from ximpol import PACKAGE_NAME
from ximpol.__version__ import TAG


_AUTHOR = 'The ximpol team'
_DESCRIPTION = 'An X-ray polarimetry simulation framework'
_LICENSE = 'GNU General Public License v3 or later'
_PACKAGES = find_packages(exclude='tests')
_URL = 'https://github.com/lucabaldini/ximpol'
_CLASSIFIERS = [
    'Intended Audience :: Science/Research',
    'License :: OSI Approved :: '
    'GNU General Public License v3 or later (GPLv3+)',
    'Operating System :: OS Independent',
    'Programming Language :: Python',
    'Intended Audience :: Science/Research',
    'Programming Language :: Python :: Implementation :: CPython',
    'Topic :: Scientific/Engineering :: Astronomy',
    'Development Status :: 2 - Pre-Alpha'
]
_SCRIPTS = glob.glob('./ximpol/bin/*.py')
_DEPENDENCIES =[
    'numpy',
    'matplotlib',
    'astropy',
    'wcsaxes',
    'scipy'
]


_KWARGS = dict(name=PACKAGE_NAME,
               version=TAG,
               author=_AUTHOR,
               description=_DESCRIPTION,
               license=_LICENSE,
               packages=_PACKAGES,
               include_package_data=True,
               url=_URL,
               classifiers=_CLASSIFIERS,
               scripts=_SCRIPTS,
               install_requires=_DEPENDENCIES)


setup(**_KWARGS)
