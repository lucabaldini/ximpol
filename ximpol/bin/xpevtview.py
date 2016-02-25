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

import matplotlib
from ximpol.utils.matplotlib_ import pyplot as plt
from ximpol.utils.matplotlib_ import context_two_by_two
from ximpol.utils.logging_ import logger, startmsg
from ximpol.irf.mrf import xAzimuthalResponseGenerator
from ximpol.evt.event import xEventFile


def xpevtview(file_path):
    """Run quick event-file viewer.
    """
    event_file = xEventFile(file_path)

    with context_two_by_two():
        fig = plt.figure()

    ax = plt.subplot(2, 2, 1)
    energy = event_file.event_data['ENERGY']
    mcenergy = event_file.event_data['MC_ENERGY']
    emin, emax = energy.min(), energy.max()
    nbins = 100
    logx = False
    if logx:
        binning = numpy.logspace(numpy.log10(emin), numpy.log10(emax), nbins)
    else:
        binning = numpy.linspace(emin, emax, nbins)
    plt.hist(energy, binning, histtype='step', label='Reconstructed E')
    plt.hist(mcenergy, binning, histtype='step', label='True E')
    if logx:
        plt.xscale('log')
    plt.yscale('log')
    plt.axis([emin, emax, None, None])
    plt.xlabel('Energy [keV]')
    plt.ylabel('Counts')
    plt.legend(bbox_to_anchor=(1.025, 0.99))

    ax = plt.subplot(2, 2, 2)
    ra0, dec0 = event_file.roi_center()
    ra = event_file.event_data['RA']
    dec = event_file.event_data['DEC']
    sidex = ra.max() - ra.min()
    sidey = dec.max() - dec.min()
    side = 0.75*max(sidex, sidey)
    nbins = 100
    ra_binning = numpy.linspace(ra0 - side, ra0 + side, nbins)
    dec_binning = numpy.linspace(dec0 - side, dec0 + side, nbins)
    binning = (ra_binning, dec_binning)
    plt.hist2d(ra, dec, binning, norm=matplotlib.colors.LogNorm())
    plt.xlabel('RA [deg]')
    plt.ylabel('Dec [deg]')
    plt.colorbar()

    ax = plt.subplot(2, 2, 3)
    time_ = event_file.event_data['TIME']
    tmin = event_file.min_good_time()
    tmax = event_file.max_good_time()
    nbins = 100
    binning = numpy.linspace(tmin, tmax, nbins)
    plt.hist(time_, binning, histtype='step')
    plt.xlabel('Time [s]')
    plt.ylabel('Energy [keV]')
    plt.axis([tmin, tmax, None, None])

    ax = plt.subplot(2, 2, 4)
    pe_angle = event_file.event_data['PE_ANGLE']
    nbins = 100
    binning = numpy.linspace(0., 2*numpy.pi, nbins)
    hist = plt.hist(pe_angle, bins=binning, histtype='step')
    #fit_results = xAzimuthalResponseGenerator.fit_histogram(hist)
    #fit_results.plot()
    plt.xlabel('PE emission angle [deg]')
    plt.axis([0., 2*numpy.pi, 0., None])

    plt.show()


if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__description__)
    parser.add_argument('file_path', type=str,
                        help='the path to the input file')
    args = parser.parse_args()
    startmsg()
    xpevtview(args.file_path)
