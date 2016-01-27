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


__description__ = 'Bin event data in different flavors'


import os
import numpy
from astropy.io import fits

from ximpol.utils.matplotlib_ import pyplot as plt
from ximpol.utils.logging_ import logger, startmsg, abort
from ximpol.irf.mrf import xAzimuthalResponseGenerator


def pha1(file_path, **kwargs):
    """
    """
    from ximpol.evt.binning import xEventBinningPHA1
    xEventBinningPHA1(file_path, **kwargs).bin_()

def pha2(event_data, output_file_path):
    """
    """
    abort('Not implemented, yet.')


def lc(event_data, output_file_path):
    """
    """
    abort('Not implemented, yet.')


def cmap(file_path, **kwargs):
    """
    """
    from ximpol.evt.binning import xEventBinningCMAP
    xEventBinningCMAP(file_path, **kwargs).bin_()

def ccube(event_data, output_file_path):
    """
    """
    abort('Not implemented, yet.')


def mcube(event_data, output_file_path):
    """
    """
    abort('Not implemented, yet.')
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
    plt.errorbar(emean, degree, yerr=err, fmt='o', label='Simulation output')
    plt.xlabel('Energy [keV]')
    plt.ylabel('Polarization degree')
    x = numpy.linspace(1, 10, 2)
    y = numpy.array([0.157]*2)
    plt.plot(x, y, label='Model input')
    plt.legend(bbox_to_anchor=(0.45, 0.95))
    plt.savefig('polarization_degree.png')
    plt.show()
    angle = [fr.phase for fr in fit_results]
    angle = numpy.degrees(numpy.array(angle))
    err = [fr.phase_error for fr in fit_results]
    err = numpy.degrees(numpy.array(err))
    plt.errorbar(emean, angle, yerr=err, fmt='o', label='Simulation output')
    plt.xlabel('Energy [keV]')
    plt.ylabel('Polarization angle [$^\\circ$]')
    x = numpy.linspace(1, 10, 2)
    y = numpy.array([161.1]*2)
    plt.plot(x, y, label='Model input')
    plt.legend(bbox_to_anchor=(0.45, 0.45))
    plt.savefig('polarization_angle.png')
    plt.show()
    """


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
    BIN_MODE_DICT[mode](file_path)



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
