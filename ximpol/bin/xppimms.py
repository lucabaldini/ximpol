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


__description__ = \
"""Calculate the minimum detectable polarization (MDP) according to the source
parameters provided through command-line switches. The program calculates the
MDP by direct numerical integration of the input spectrum (i.e., there are no
event lists involved). As a consequence, part of the richness of the detector
response (most notably, the energy dispersion and the effective area
vignetting) is not captured here. For use cases beyond simple stationary point
sources, the use of xpobbsim and xpbin in MCUBE mode are recommended, as this
approach offers the maximum flexibility.
"""

import os
import numpy
import imp

from ximpol.irf import load_arf, load_mrf, DEFAULT_IRF_NAME
from ximpol.srcmodel.spectrum import xCountSpectrum
from ximpol.srcmodel.roi import xPointSource, xROIModel
from ximpol.utils.logging_ import logger, startmsg, abort
from ximpol.core.spline import xInterpolatedUnivariateSplineLinear
from ximpol.core.spline import xInterpolatedBivariateSplineLinear
from ximpol.srcmodel.spectrum import power_law
from ximpol.srcmodel.polarization import constant
from xpmdp import _make_energy_binning

EBIN_ALGS = ['FILE', 'LIN', 'LOG', 'LIST']


"""Command-line switches.
"""
import ast
import argparse

formatter = argparse.ArgumentDefaultsHelpFormatter
PARSER = argparse.ArgumentParser(description=__description__,
                                 formatter_class=formatter)
PARSER.add_argument('--outfile', type=str, default=None,
                    help='the output text file')
PARSER.add_argument('--irfname', type=str, default=DEFAULT_IRF_NAME,
                    help='the input configuration file')
PARSER.add_argument('--norm', type=float, default=10.,
                    help='the power-law normalization for energy spectrum')
PARSER.add_argument('--index', type=float, default=2.,
                    help='the power-law index for energy spectrum')
PARSER.add_argument('--nH', type=float, default=0.,
                    help='the hydrogen column density (cm^{-2})')
PARSER.add_argument('--redshift', type=float, default=0.,
                    help='the cosmological redshift')
PARSER.add_argument('--duration', type=float, default=10000,
                    help='the duration (in s) of the simulation')
PARSER.add_argument('--tstart', type=float, default=0.,
                    help='the start time (MET in s) of the simulation')
PARSER.add_argument('--ebinalg', choices=EBIN_ALGS, default='LIN',
                    help='energy binning specification')
PARSER.add_argument('--emin', type=float, default=2.,
                    help='minimum energy for LIN/LOG energy binning')
PARSER.add_argument('--emax', type=float, default=8.,
                    help='maximum energy for LIN/LOG energy binning')
PARSER.add_argument('--ebins', type=int, default=4,
                    help='number of bins for LIN/LOG energy binning')
PARSER.add_argument('--ebinfile', type=str, default=None,
                    help='path to the optional energy bin definition file')
PARSER.add_argument('--ebinning', type=ast.literal_eval, default=None,
                    help='the list containing the bin edges')
PARSER.add_argument('--clobber', type=ast.literal_eval, choices=[True, False],
                    default=True,
                    help='overwrite or do not overwrite existing output files')

def _build_source(**kwargs):
    """Build the source roi for the MDP calculation.
    """
    energy_spectrum = power_law(kwargs['norm'], kwargs['index'])
    polarization_degree = constant(0.)
    polarization_angle = constant(0.)
    column_density = kwargs['nH']
    redshift = kwargs['redshift']
    source = xPointSource('Point source', 0., 0., energy_spectrum,
                          polarization_degree, polarization_angle,
                          column_density, redshift)
    return source

def xppimms(**kwargs):
    """Calculate the MDP.
    """
    logger.info('Loading the instrument response functions...')
    aeff = load_arf(kwargs['irfname'])
    modf = load_mrf(kwargs['irfname'])
    
    # Build the source
    source = _build_source(**kwargs)
    logger.info(source)
    
    # This is copied from xpobbsim and should probably be factored out.
    # Actually, this should be a method of the ROI class. TBD
    if kwargs['tstart'] < source.min_validity_time:
        kwargs['tstart'] = source.min_validity_time
    logger.info('Simulation start time set to %s...' % kwargs['tstart'])
    tstop = kwargs['tstart'] + kwargs['duration']
    if tstop > source.max_validity_time:
        tstop = source.max_validity_time
        logger.info('Simulation stop time set to %s...' % tstop)
    kwargs['tstop'] = tstop
    
    tsamples = source.sampling_time(kwargs['tstart'], tstop)
    logger.info('Sampling times: %s' % tsamples)
    count_spectrum = xCountSpectrum(source.energy_spectrum, aeff, tsamples,
                                    source.column_density, source.redshift)

    # Do the actual work.
    ebinning =_make_energy_binning(**kwargs)
    mdp_table = count_spectrum.build_mdp_table(ebinning, modf)
    logger.info(mdp_table)
    file_path = kwargs['outfile']
    if file_path is not None:
        logger.info('Writing output file path %s...' % file_path)
        open(file_path, 'w').write('%s\n\n%s' % (kwargs, mdp_table))
        logger.info('Done.')
    return mdp_table

if __name__=='__main__':
    args = PARSER.parse_args()
    startmsg()
    xppimms(**args.__dict__)
