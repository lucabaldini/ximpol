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
PARSER.add_argument('--clobber', type=ast.literal_eval, choices=[True, False],
                    default=True,
                    help='overwrite or do not overwrite existing output files')



def _make_energy_binning(**kwargs):
    """Build the energy binning for the MDP calculation.

    While there's surely some overlap with the code in ximpol.evt.binning
    module, none of the binning methods implemented there is exactly what
    we need here---the closest being that for the modulation cube, which
    is also supporting the extra EQP binning mode that really does not make
    a lot of sense here (we don't have an event file with the column energy).
    I guess we'll just go along with some code duplication. Too bad.
    """
    ebinalg = kwargs['ebinalg']
    emin = kwargs['emin']
    emax = kwargs['emax']
    ebins = kwargs['ebins']
    ebinning = kwargs['ebinning']
    if ebinalg == 'LIN':
        ebinning = numpy.linspace(emin, emax, ebins + 1)
    elif ebinalg == 'LOG':
        ebinning = numpy.linspace(numpy.log10(emin), numpy.log10(emax),
                                  ebins + 1)
    elif ebinalg == 'FILE':
        ebinfile = self.get('ebinfile')
        assert ebinfile is not None
        ebinning = numpy.loadtxt(ebinfile)
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

    # This is copied from roi.py and should probably be factored out.
    # Again, the ROI class should be able to sum the count spectra of all the
    # component and expose the result.
    sources = ROI_MODEL.values()
    if len(sources) > 1:
        abort('Multiple sources not implemented, yet.')
    source = sources[0]
    if isinstance(source, xPeriodicPointSource):
        observation_time = kwargs['tstop'] - kwargs['tstart']
        psamples = numpy.linspace(kwargs['phasemin'], kwargs['phasemax'], 100)
        logger.info('Sampling phases: %s' % psamples)
        count_spectrum = xCountSpectrum(source.energy_spectrum, aeff, psamples,
                                        scale_factor=observation_time,
                                        column_density=source.column_density,
                                        redshift=source.redshift)
    else:
        tsamples = source.sampling_time(kwargs['tstart'], kwargs['tstop'])
        logger.info('Sampling times: %s' % tsamples)
        count_spectrum = xCountSpectrum(source.energy_spectrum, aeff, tsamples,
                                        column_density=source.column_density,
                                        redshift=source.redshift)

    # Do the actual work.
    ebinning =_make_energy_binning(**kwargs)
    mdp_table = count_spectrum.build_mdp_table(ebinning, modf)
    logger.info(mdp_table)
    return mdp_table



if __name__=='__main__':
    args = PARSER.parse_args()
    startmsg()
    xpmdp(**args.__dict__)
