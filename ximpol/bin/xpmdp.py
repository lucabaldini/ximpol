#!/usr/bin/env python
#
# Copyright (C) 2016, the ximpol team.
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


__description__ = 'Calculate the MPD for a given source model'


import os
import numpy
import imp

from ximpol.irf import load_arf, load_mrf
from ximpol.srcmodel.spectrum import xCountSpectrum
from ximpol.srcmodel.roi import xPeriodicPointSource
from ximpol.utils.logging_ import logger, startmsg, abort
from ximpol.core.spline import xInterpolatedUnivariateSplineLinear
from ximpol.core.spline import xInterpolatedBivariateSplineLinear


EBIN_ALGS = ['FILE', 'LIN', 'LOG', 'LIST']


"""Command-line switches.
"""
import ast
import argparse

formatter = argparse.ArgumentDefaultsHelpFormatter
PARSER = argparse.ArgumentParser(description=__description__,
                                 formatter_class=formatter)
#PARSER.add_argument('--outfile', type=str, default=None,
#                    help='the output FITS event file')
PARSER.add_argument('--configfile', type=str, required=True,
                    help='the input configuration file')
PARSER.add_argument('--irfname', type=str, default='xipe_baseline',
                    help='the input configuration file')
PARSER.add_argument('--duration', type=float, default=10,
                    help='the duration (in s) of the simulation')
PARSER.add_argument('--tstart', type=float, default=0.,
                    help='the start time (MET in s) of the simulation')
PARSER.add_argument('--phasemin', type=float, default=0.,
                    help='minimum phase')
PARSER.add_argument('--phasemax', type=float, default=1.,
                    help='maximum phase')
PARSER.add_argument('--ebinalg', choices=EBIN_ALGS, default='LIN',
                    help='energy binning specification')
PARSER.add_argument('--emin', type=float, default=1.,
                    help='minimum energy for LIN/LOG energy binning')
PARSER.add_argument('--emax', type=float, default=10.,
                    help='maximum energy for LIN/LOG energy binning')
PARSER.add_argument('--ebins', type=int, default=5,
                    help='number of bins for LIN/LOG energy binning')
PARSER.add_argument('--ebinfile', type=str, default=None,
                    help='path to the optional energy bin definition file')
PARSER.add_argument('--ebinning', type=ast.literal_eval, default=None,
                    help='the list containing the bin edges')
#PARSER.add_argument('--clobber', type=ast.literal_eval, choices=[True, False],
#                    default=True,
#                    help='overwrite or do not overwrite existing output files')



def _make_binning(ebinalg, emin=1., emax=10., ebins=5, ebinning=None):
    """Build the modulation cube binning.

    Warning
    -------
    This is copied from evt/binning.py and should be factored out.
    """
    if ebinalg == 'LIN':
        ebinning = numpy.linspace(emin, emax, ebins + 1)
    elif ebinalg == 'LOG':
        ebinning = numpy.linspace(numpy.log10(emin), numpy.log10(emax),
                                  ebins + 1)
    #elif ebinalg == 'FILE':
    #    ebinfile = self.get('ebinfile')
    #    assert ebinfile is not None
    #    ebinning = self.read_binning(ebinfile)
    elif ebinalg == 'LIST':
        assert isinstance(ebinning, list)
        ebinning = numpy.array(ebinning, 'd')
    else:
        abort('ebinalg %s not implemented yet' % ebinalg)
    return ebinning

def xpmdp(**kwargs):
    """Calculate the MDP.
    """
    logger.info('Loading the instrument response functions...')
    aeff = load_arf(kwargs['irfname'])
    modf = load_mrf(kwargs['irfname'])
    module_name = os.path.basename(kwargs['configfile']).replace('.py', '')
    ROI_MODEL = imp.load_source(module_name, kwargs['configfile']).ROI_MODEL
    logger.info(ROI_MODEL)

    # This is copied from xpobbsim and should probably be factored out.
    # Actually, this should be a method of the ROI class. TBD
    if kwargs['tstart'] < ROI_MODEL.min_validity_time():
        kwargs['tstart'] = ROI_MODEL.min_validity_time()
        logger.info('Simulation start time set to %s...' % kwargs['tstart'])
    tstop = kwargs['tstart'] + kwargs['duration']
    if tstop > ROI_MODEL.max_validity_time():
        tstop = ROI_MODEL.max_validity_time()
        logger.info('Simulation stop time set to %s...' % tstop)
    kwargs['tstop'] = tstop
    observation_time = kwargs['tstop'] - kwargs['tstart']

    # This is copied from roi.py and should probably be factored out.
    # Again, the ROI class should be able to sum the count spectra of all the
    # component and expose the result.
    sources = ROI_MODEL.values()
    if len(sources) > 1:
        abort('Multiple sources not implemented, yet.')
    source = sources[0]
    if isinstance(source, xPeriodicPointSource):
        psamples = numpy.linspace(kwargs['phasemin'], kwargs['phasemax'], 100)
        logger.info('Sampling phases: %s' % psamples)
        scale_fact = observation_time/source.ephemeris.period(kwargs['tstart'])
        count_spectrum = xCountSpectrum(source.energy_spectrum, aeff, psamples,
                                        scale=scale_fact)
        time_integrated_spectrum = count_spectrum.build_time_integral()
    else:
        tsamples = source.sampling_time(kwargs['tstart'], kwargs['tstop'])
        logger.info('Sampling times: %s' % tsamples)
        count_spectrum = xCountSpectrum(source.energy_spectrum, aeff, tsamples)
        time_integrated_spectrum = count_spectrum.build_time_integral()

    # Thuis should be a callable method in the binning module.
    ebinning =_make_binning(kwargs['ebinalg'], kwargs['emin'], kwargs['emax'],
                            kwargs['ebins'], kwargs['ebinning'])

    # And this might be implemented in the irf.mrf module.
    _x = time_integrated_spectrum.x
    _y = time_integrated_spectrum.y*modf(_x)
    mu_spectrum = xInterpolatedUnivariateSplineLinear(_x, _y)

    for _emin, _emax in zip(ebinning[:-1], ebinning[1:]) +\
        [(ebinning[0], ebinning[-1])]:
        num_counts = count_spectrum.num_expected_counts(emin=_emin, emax=_emax)
        mu_average = mu_spectrum.integral(_emin, _emax)/num_counts
        mdp = 4.29/mu_average/numpy.sqrt(num_counts)
        logger.info('%.2f--%.2f keV: %d counts in %d s, MDP %.2f%%' %\
                    (_emin, _emax, num_counts, observation_time, 100*mdp))




if __name__=='__main__':
    args = PARSER.parse_args()
    startmsg()
    xpmdp(**args.__dict__)
