#!/usr/bin/env python
#
# Copyright (C) 2016, the ximpol team.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
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


import os

from ximpol import XIMPOL_CONFIG, XIMPOL_DATA, XIMPOL_DOC
from ximpol.core.pipeline import xPipeline
from ximpol.utils.logging_ import logger
from ximpol.config.single_point_source import PL_NORM, PL_INDEX


"""Script-wide simulation and analysis settings.
"""
CFG_FILE = os.path.join(XIMPOL_CONFIG, 'single_point_source.py')
OUT_FILE_PATH_BASE = os.path.join(XIMPOL_DATA, 'single_point_source')
EVT_FILE_PATH = '%s.fits' % OUT_FILE_PATH_BASE
SIM_DURATION = 10000.
OUTPUT_FOLDER = os.path.join(XIMPOL_DOC, 'figures', 'showcase')


"""Main pipeline object.
"""
PIPELINE = xPipeline(clobber=False)


def run():
    PIPELINE.xpobssim(configfile=CFG_FILE, duration=SIM_DURATION,
                      outfile=EVT_FILE_PATH)
    pha1_file_path = PIPELINE.xpbin(EVT_FILE_PATH, algorithm='PHA1')
    spec_fitter = PIPELINE.xpxspec(pha1_file_path)
    (index, index_err), (norm, norm_err) = spec_fitter.fit_parameters()
    logger.info('Fitted PL norm = %.4f +- %4f (input = %.4f)' %\
                (norm, norm_err, PL_NORM))
    logger.info('Fitted PL index = %.4f +- %4f (input = %.4f)' %\
                (index, index_err, PL_INDEX))


if __name__ == '__main__':
    run()
