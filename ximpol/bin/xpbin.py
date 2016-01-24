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


__description__ = 'Bin event data'


import os
import numpy
from astropy.io import fits
from astropy import wcs

from ximpol.utils.matplotlib_ import pyplot as plt
from ximpol.utils.logging_ import logger, startmsg, abort
from ximpol.irf.base import xPrimaryHDU, update_header
from ximpol.evt.binning import xColDefsSpectrum, SPECTRUM_HEADER_SPECS
from ximpol.irf.mrf import xAzimuthalResponseGenerator


def pha1(event_data, output_file_path):
    """
    """
    # Horrible---need to put the EBOUND information in the event file.
    _num_chans = 256
    # And also the GTI
    _elapsed_time = event_data['TIME'][-1] - event_data['TIME'][0]
    # Moving on...
    _binning = numpy.linspace(-0.5, _num_chans - 0.5, _num_chans)
    n, bins, patches = plt.hist(event_data['PHA'], bins=_binning)
    primary_hdu = xPrimaryHDU()
    data = [numpy.arange(_num_chans),
            n/_elapsed_time,
            numpy.sqrt(n)/_elapsed_time
    ]
    cols = xColDefsSpectrum(data)
    spectrum_hdu = fits.BinTableHDU.from_columns(cols)
    update_header(spectrum_hdu, SPECTRUM_HEADER_SPECS)
    hdulist = fits.HDUList([primary_hdu, spectrum_hdu])
    hdulist.info()
    hdulist.writeto(output_file_path, clobber=True)
    logger.info('Binned (PHA1) written to %s...' % output_file_path)


def pha2(event_data, output_file_path):
    """
    """
    abort('Not implemented, yet.')


def lc(event_data, output_file_path):
    """
    """
    abort('Not implemented, yet.')


def cmap(event_data, output_file_path, nside=256, mc=False):
    """
    """
    if mc:
        ra = event_data['MC_RA']
        dec = event_data['MC_DEC']
    else:
        ra = event_data['RA']
        dec = event_data['DEC']
    # Horrible: this should be written in the event file at simulation time.
    ra0 = 0.5*(ra.max() + ra.min())
    dec0 = 0.5*(dec.max() + dec.min())
    side = 5./60
    # Moving on...
    delta = side/nside
    binning = numpy.linspace(0, nside, nside + 1)
    # Build the WCS object
    w = wcs.WCS(naxis=2)
    w.wcs.crpix = [nside, 0.]
    w.wcs.cdelt = [-delta, delta]
    w.wcs.crval = [ra0 - 0.5*side, dec0 - 0.5*side]
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.equinox = 2000
    header = w.to_header()
    pix = w.wcs_world2pix(zip(ra, dec), 1)
    n, x, y = numpy.histogram2d(pix[:,1], pix[:,0], bins=(binning, binning))
    hdu = fits.PrimaryHDU(n, header=header)
    hdu.writeto(output_file_path, clobber=True)
    logger.info('Binned (CMAP) written to %s...' % output_file_path)

def ccube(event_data, output_file_path):
    """
    """
    abort('Not implemented, yet.')


def mcube(event_data, output_file_path):
    """
    """
    energy = event_data['MC_ENERGY']
    phi = event_data['PE_ANGLE']
    ebinning = numpy.linspace(1, 10, 10)
    phi_binning = numpy.linspace(0, 2*numpy.pi, 100)
    fit_results = []
    emean = []
    for i, (_emin, _emax) in enumerate(zip(ebinning[:-1], ebinning[1:])):
        #_emean = 0.5*(_emin + _emax)

        _mask = (energy > _emin)*(energy < _emax)
        _energy = energy[_mask]
        _phi = phi[_mask]
        _emean = _energy.sum()/len(_energy)
        emean.append(_emean)
        _hist = plt.hist(_phi, bins=phi_binning, histtype='step')
        _fr = xAzimuthalResponseGenerator.fit_histogram(_hist)
        _fr.emean = _emean
        fit_results.append(_fr)
        _fr.plot(label='Energy: %.2f--%.2f keV' % (_emin, _emax))
        plt.axis([0., 2*numpy.pi, 0., 1.2*_hist[0].max()])
        #plt.show()
        plt.savefig('polarization_fit%d.png' % i)
        plt.clf()

    from ximpol.irf.mrf import xModulationFactor
    from ximpol import XIMPOL_IRF
    file_path = os.path.join(XIMPOL_IRF,'fits','xipe_baseline.mrf')
    modf = xModulationFactor(file_path)

    emean = numpy.array(emean)
    degree = [fr.visibility for fr in fit_results]
    degree = numpy.array(degree)
    degree /= modf(emean)
    err = [fr.visibility_error for fr in fit_results]
    plt.errorbar(emean, degree, yerr=err, fmt='o')
    plt.xlabel('Energy [keV]')
    plt.ylabel('Polarization degree')
    x = numpy.linspace(1, 10, 2)
    y = numpy.array([0.157]*2)
    plt.plot(x, y)
    plt.savefig('polarization_degree.png')
    plt.show()


    angle = [fr.phase for fr in fit_results]
    angle = numpy.degrees(numpy.array(angle))
    err = [fr.phase_error for fr in fit_results]
    err = numpy.degrees(numpy.array(err))
    plt.errorbar(emean, angle, yerr=err, fmt='o')
    plt.xlabel('Energy [keV]')
    plt.ylabel('Polarization angle [$^\\circ$]')
    x = numpy.linspace(1, 10, 2)
    y = numpy.array([161.1]*2)
    plt.plot(x, y)
    plt.savefig('polarization_angle.png')
    plt.show()



BIN_MODE_DICT = {'PHA1' : pha1,
                 'PHA2' : pha2,
                 'LC'   : lc,
                 'CMAP' : cmap,
                 'CCUBE': ccube,
                 'MCUBE': mcube
}
BIN_MODES = BIN_MODE_DICT.keys()
BIN_MODES.sort()


def xpbin(file_path, mode, output_file_path):
    """Application to bin the data.
    """
    assert(file_path.endswith('.fits'))
    logger.info('Opening input file %s...' % file_path)
    hdu_list = fits.open(file_path)
    hdu_list.info()
    event_data = hdu_list['EVENTS'].data
    if output_file_path is None:
        output_file_path = file_path.replace('.fits', '_%s.fits' % mode)
    BIN_MODE_DICT[mode](event_data, output_file_path)



if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__description__)
    parser.add_argument('file_path', type=str,
                        help='the path to the input event file')
    parser.add_argument('-o', '--output-file', type=str, default=None,
                        help='the output binned FITS file')
    parser.add_argument('-m', '--mode', choices=BIN_MODES, required=True,
                        help='the binning mode')
    args = parser.parse_args()
    startmsg()
    xpbin(args.file_path, args.mode, args.output_file)
