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


__description__ = 'Quick viewer for an output xpobbsim event file'


import os
import numpy
from astropy.io import fits

from ximpol.utils.matplotlib_ import pyplot as plt
from ximpol.utils.matplotlib_ import context_two_by_two
from ximpol.utils.logging_ import logger, startmsg


def xpevtview(file_path):
    """Run quick event-file viewer.
    """
    logger.info('Opening input file %s...' % file_path)
    hdu_list = fits.open(file_path)
    hdu_list.info()
    evtdata = hdu_list['EVENTS'].data

    with context_two_by_two():
        fig = plt.figure()

    ax = plt.subplot(2, 2, 1)
    plt.hist(evtdata['ENERGY'], bins=numpy.logspace(0,1,50), histtype='step')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Energy [keV]')

    ax = plt.subplot(2, 2, 2)
    plt.hist2d(evtdata['RA'], evtdata['DEC'], 100)
    plt.xlabel('RA [deg]')
    plt.ylabel('Dec [deg]')

    ax = plt.subplot(2, 2, 3)
    plt.plot(evtdata['TIME'], evtdata['ENERGY'], 'o')
    plt.yscale('log')
    plt.xlabel('Time [s]')
    plt.ylabel('Energy [keV]')

    ax = plt.subplot(2, 2, 4)
    pe_angle = evtdata['PE_ANGLE']*(180./numpy.pi)
    le_mask = evtdata['ENERGY'] < 5
    he_mask = evtdata['ENERGY'] > 5
    binning = numpy.linspace(0, 360, 50)
    plt.hist(pe_angle[le_mask], bins=binning, histtype='step')
    plt.hist(pe_angle[he_mask], bins=binning, histtype='step')
    plt.xlabel('PE emission angle [deg]')

    plt.show()


if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__description__)
    parser.add_argument('file_path', type=str,
                        help='the path to the input file')
    args = parser.parse_args()
    startmsg()
    xpevtview(args.file_path)
