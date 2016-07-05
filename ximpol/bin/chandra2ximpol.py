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


__description__ = 'Chandra-to-ximpol converter'


import os
import numpy
import imp
import pyregion
import pyregion._region_filter as filter

from astropy import wcs
from astropy.io import fits
from astropy.coordinates import SkyCoord
from ximpol.irf import load_irfs, DEFAULT_IRF_NAME
from ximpol.utils.os_ import mkdir
from ximpol.utils.profile import xChrono
from ximpol import XIMPOL_IRF, XIMPOL_DATA
from ximpol.srcmodel.roi import xROIModel
from ximpol.srcmodel.polarization import constant
from ximpol.evt.event import xMonteCarloEventList
from ximpol.utils.logging_ import logger, startmsg, abort
from ximpol.core.spline import xInterpolatedUnivariateSplineLinear, \
                               xInterpolatedBivariateSplineLinear


"""Command-line switches.
"""
import ast
import argparse

formatter = argparse.ArgumentDefaultsHelpFormatter
PARSER = argparse.ArgumentParser(description=__description__,
                                 formatter_class=formatter)
PARSER.add_argument('evfile', type=str,
                    help='path to the input FITS event file')
PARSER.add_argument('--outfile', type=str, default=None,
                    help='the output FITS event file')
PARSER.add_argument('--configfile', type=str, default=None,
                    help='the input configuration file')
PARSER.add_argument('--regfile', type=str, default=None,
                    help='the input region file')
PARSER.add_argument('--mc', type=ast.literal_eval, choices=[True, False],
                    default=True,
                    help='use chandra coordinates for source id definition')
PARSER.add_argument('--duration', type=float, default=numpy.nan,
                    help='the duration (in s) of the simulation')
PARSER.add_argument('--irfname', type=str, default=DEFAULT_IRF_NAME,
                    help='the input configuration file')
PARSER.add_argument('--acis', type=str, default='i', choices=['i', 's'],
                    help='the chandra acis detector name')
PARSER.add_argument('--seed', type=int, default=0,
                    help='the random seed for the simulation')
PARSER.add_argument('--clobber', type=ast.literal_eval, choices=[True, False],
                    default=True,
                    help='overwrite or do not overwrite existing output files')

class xSimulationInfo:

    """Empty container to pass along all the relevant information about the
    simulation.
    """

    pass

def filter_region(region, col_mc_ra, col_mc_dec):
    """Check which events are inside the defined region and return the
    corresponding array mask.
    """
    radec = numpy.column_stack((col_mc_ra, col_mc_dec))
    if region.name is 'circle':
        ra, dec, rad = region.coord_list
        reg_filter = filter.Circle(ra, dec, rad)
    elif region.name is 'polygon':
        list_coord = region.coord_list
        reg_filter = filter.Polygon(list_coord[::2], list_coord[1::2])
    elif region.name is 'box':
        ra, dec, width, height = region.coord_list
        reg_filter = filter.Box(ra, dec, width, height)
    else:
        logger.info('Region shape not implemented yet. Ignoring it...')
        return numpy.zeros(len(col_mc_ra), dtype=bool)
    return reg_filter.inside(radec)

def _get_radec(hdr, tbdata):
    """Make the conversion from pixel coordinates to Ra-Dec.
    """
    # Read wcs parameters from header.
    index = hdr.values().index('EQPOS')
    wcs_list = hdr[index+2:index+11].values()

    # Create and set up a WCS object.
    w = wcs.WCS(naxis=2)
    w.wcs.crpix = wcs_list[2:8:5]
    w.wcs.cdelt = numpy.array(wcs_list[3:9:5])
    w.wcs.crval = wcs_list[1:7:5]
    w.wcs.ctype = wcs_list[0:6:5]
    w.wcs.equinox = 2000.
    w.wcs.radesys = 'ICRS'

    # Read the pixel coordinates to be converted.
    x = tbdata['x']
    y = tbdata['y']
    pixcrd = numpy.column_stack((x,y))

    # Convert pixel coordinates to world coordinates (Ra, Dec).
    ra_dec = w.wcs_pix2world(pixcrd, 1)
    return ra_dec[:,0], ra_dec[:,1]

def _load_chandra_arf(arf_file):
    """Load the Chandra effective area data from file.
    """
    hdu_list = fits.open(arf_file)
    tbdata = hdu_list['SPECRESP'].data
    _x = 0.5*(tbdata.field('ENERG_LO') + tbdata.field('ENERG_HI'))
    _y = tbdata.field('SPECRESP')
    fmt = dict(xname='Energy', xunits='keV', yname='Effective area',
               yunits='cm$^2$')
    return xInterpolatedUnivariateSplineLinear(_x, _y, **fmt)

def _load_chandra_vign(vign_file):
    """Load the Chandra vignetting data from file.
    """
    hdu_list = fits.open(vign_file)
    tbdata = hdu_list[1].data
    _x = 0.5*(tbdata.field('ENERG_LO') + tbdata.field('ENERG_HI'))[0,:]
    _y = tbdata.field('THETA')[0,:]
    _vignet = tbdata.field('VIGNET')
    _z = numpy.mean(_vignet[0,:,:,:], axis=0)
    _mask = _y <= 10 #arcmin
    _y = _y[_mask]
    _z = _z[_mask, :]
    fmt = dict(xname='Energy', xunits='keV', yname='Off-axis angle',
               yunits='arcmin', zname='Vignetting')
    return xInterpolatedBivariateSplineLinear(_x, _y, _z.T, **fmt)

def _make_aeff_ratio(aeff, c_aeff, c_vign, emin=0.3, emax=10.):
    """Make the ratio between ximpol and Chandra effective area taking in
    account the vignetting.
    """
    _x = aeff.x[(aeff.x >= emin)*(aeff.x <= emax)]
    _aeff_ratio = aeff(_x)/c_aeff(_x)
    _y = aeff.vignetting.y
    _x_mesh, _y_mesh = numpy.meshgrid(_x, _y)
    _vign_ratio = (aeff.vignetting(_x_mesh, _y_mesh)/c_vign(_x_mesh, _y_mesh)).T
    _z = (_vign_ratio.T*_aeff_ratio).T
    fmt = dict(xname='Energy', xunits='keV', yname='Off-axis angle',
                                yunits='arcmin', zname='Effective area ratio')
    return xInterpolatedBivariateSplineLinear(_x, _y, _z, **fmt)

def _time_scaling(scale, col_mc_energy, col_time, col_mc_ra, col_mc_dec):
    """Make the scaling from Chandra observation time to an arbitrary time
    provided with duration parameter.
    """
    mc_energy = numpy.array([])
    time = numpy.array([])
    mc_ra = numpy.array([])
    mc_dec = numpy.array([])
    delta_time = col_time[-1] - col_time[0]
    for i in range(0, int(scale[1])):
        mc_energy = numpy.append(mc_energy, col_mc_energy)
        time = numpy.append(time, i*delta_time + col_time)
        mc_ra = numpy.append(mc_ra, col_mc_ra)
        mc_dec = numpy.append(mc_dec, col_mc_dec)
    index = 1 + int(scale[0]*len(col_mc_energy))
    mc_energy = numpy.append(mc_energy, col_mc_energy[:index])
    time = numpy.append(time, (i+1)*delta_time + col_time[:index])
    mc_ra = numpy.append(mc_ra, col_mc_ra[:index])
    mc_dec = numpy.append(mc_dec, col_mc_dec[:index])
    return mc_energy, time, mc_ra, mc_dec

def chandra2ximpol(file_path, **kwargs):
    """Make the conversion from Chandra to ximpol.
    """
    assert(file_path.endswith('.fits'))
    if kwargs['outfile'] is None:
        outfile = os.path.basename(file_path).replace('.fits','_xipe.fits')
        mkdir(XIMPOL_DATA)
        kwargs['outfile'] = os.path.join(XIMPOL_DATA, outfile)
        logger.info('Setting output file path to %s...' % kwargs['outfile'])
    if os.path.exists(kwargs['outfile']) and not kwargs['clobber']:
        logger.info('Output file %s already exists.' % kwargs['outfile'])
        logger.info('Remove the file or set "clobber = True" to overwite it.')
        return kwargs['outfile']

    chrono = xChrono()
    logger.info('Setting the random seed to %d...' % kwargs['seed'])
    numpy.random.seed(kwargs['seed'])
    logger.info('Loading the instrument response functions...')
    aeff, psf, modf, edisp = load_irfs(kwargs['irfname'])
    c_aeff_name = 'chandra_acis_%s.arf' % kwargs['acis']
    c_aeff_file = os.path.join(XIMPOL_IRF, 'fits', c_aeff_name)
    logger.info('Reading Chandra effective area data from %s...' % c_aeff_file)
    c_aeff = _load_chandra_arf(c_aeff_file)
    c_vign_name = 'chandra_vignet.fits'
    c_vign_file = os.path.join(XIMPOL_IRF, 'fits', c_vign_name)
    logger.info('Reading Chandra vignetiing data from %s...' % c_vign_file)
    c_vign = _load_chandra_vign(c_vign_file)
    aeff_ratio = _make_aeff_ratio(aeff, c_aeff, c_vign)
    logger.info('Done %s.' % chrono)

    logger.info('Loading the input FITS event file...')
    hdu_list = fits.open(file_path)
    hdr = hdu_list[1].header
    tbdata = hdu_list[1].data
    tstart = hdr['TSTART']
    tstop = hdr['TSTOP']
    obs_time = tstop - tstart
    logger.info('Chandra observation time: %i s.' % obs_time)
    ra_pnt = hdr['RA_PNT']
    dec_pnt = hdr['DEC_PNT']
    ROI_MODEL = xROIModel(ra_pnt, dec_pnt)

    logger.info('Reading Chandra data...')
    col_mc_energy = tbdata['energy']*0.001 # eV -> keV
    col_mc_ra, col_mc_dec = _get_radec(hdr, tbdata)
    ref_skyccord = SkyCoord(ra_pnt, dec_pnt, unit='deg')
    evt_skycoord = SkyCoord(col_mc_ra, col_mc_dec, unit='deg')
    separation = evt_skycoord.separation(ref_skyccord).arcmin

    logger.info('Converting from Chandra to ximpol...')
    rnd_ratio = numpy.random.random(len(col_mc_energy))
    # The condition col_mc_energy < 10. is needed to avoid to take the bunch of
    # events with energy > 10 keV included into the Chandra photon list (we
    # actually don't know the reason).
    _mask = (rnd_ratio < aeff_ratio(col_mc_energy, separation))*\
                                                            (col_mc_energy<10.)
    # This is needed for over-sample the events in case of effective area ratio
    # greater than 1.
    _mask_ratio = (rnd_ratio < (aeff_ratio(col_mc_energy, separation)-1.))*\
                                                            (col_mc_energy<10.)
    col_mc_energy = numpy.append(col_mc_energy[_mask],
                                                    col_mc_energy[_mask_ratio])
    col_mc_ra = numpy.append(col_mc_ra[_mask], col_mc_ra[_mask_ratio])
    col_mc_dec = numpy.append(col_mc_dec[_mask], col_mc_dec[_mask_ratio])
    col_time = numpy.append(tbdata['time'][_mask],tbdata['time'][_mask_ratio])

    # If duration parameter is provided the counts are down- or oversampled.
    duration = kwargs['duration']
    if not numpy.isnan(duration):
        logger.info('Setting the observation time to %d s...' % duration)
        scale = numpy.modf(duration/obs_time)
        tstop = tstart + duration
        logger.info('Scaling counts according to observation time...')
        col_mc_energy, col_time, col_mc_ra, col_mc_dec = _time_scaling(scale,
                                col_mc_energy, col_time, col_mc_ra, col_mc_dec)

    # The Ra Dec coordinates are calculated here because they are needed in
    # source id definition (in case of kwargs['chandra']==False)
    col_ra, col_dec = psf.smear(col_mc_ra, col_mc_dec)
    # The default source id in case of no regfile and for regions not selected
    # is zero. For regions selected in the regfile the source id is determined
    # by the order of definition (starting with 1).
    logger.info('Defining source id...')
    src_id = numpy.zeros(len(col_mc_dec))
    n_reg = -1
    if kwargs['regfile'] is not None:
        regions = pyregion.open(kwargs['regfile'])
        for n_reg, region in enumerate(regions):
            if kwargs['mc'] is True:
                mask = filter_region(region, col_mc_ra, col_mc_dec)
            else:
                mask = filter_region(region, col_ra, col_dec)
            src_id[mask]= n_reg + 1

    # If configuration file is not provided we assume a non-polarized source.
    # In case of configfile with polarization model defined in different
    # regions the photoelectron distribution is generated according to them.
    if kwargs['configfile'] is None:
        logger.info('Configuration file not provided.')
        logger.info('Setting polarization angle and degree to zero...')
        polarization_degree = constant(0.)
        polarization_angle = constant(0.)
        pol_dict = {0: [polarization_degree, polarization_angle]}
        for k in range(1,n_reg+2):
            pol_dict[k] = [polarization_degree, polarization_angle]
    else:
        logger.info('Setting up the polarization source model...')
        module_name = os.path.basename(kwargs['configfile']).replace('.py', '')
        pol_dict = imp.load_source(module_name,
                                        kwargs['configfile']).POLARIZATION_DICT
    logger.info('Generating photoelectron azimuthal distribution...')
    col_pe_angle = numpy.empty(len(col_mc_dec))
    for key, pol_list in pol_dict.items():
        _mask_src = src_id == key
        pol_degree = pol_list[0](col_mc_energy[_mask_src], col_time[_mask_src],
                                 col_mc_ra[_mask_src], col_mc_dec[_mask_src])
        pol_angle = pol_list[1](col_mc_energy[_mask_src], col_time[_mask_src],
                                col_mc_ra[_mask_src], col_mc_dec[_mask_src])
        col_pe_angle[_mask_src] = modf.rvs_phi(col_mc_energy[_mask_src],
                                               pol_degree, pol_angle)

    gti_list = [(tstart, tstop)]
    event_list = xMonteCarloEventList()
    event_list.set_column('TIME', col_time)
    event_list.set_column('MC_ENERGY', col_mc_energy)
    col_pha = edisp.matrix.rvs(col_mc_energy)
    event_list.set_column('PHA', col_pha)
    col_energy = edisp.ebounds(col_pha)
    event_list.set_column('ENERGY',col_energy)
    event_list.set_column('MC_RA', col_mc_ra)
    event_list.set_column('MC_DEC', col_mc_dec)
    event_list.set_column('RA', col_ra)
    event_list.set_column('DEC', col_dec)
    event_list.set_column('MC_SRC_ID', src_id)
    event_list.set_column('PE_ANGLE', col_pe_angle)
    # Set the phase to rnd [0-1] for all non-periodic sources.
    phase=numpy.random.uniform(0,1,len(col_dec))
    event_list.set_column('PHASE', phase)
    logger.info('Done %s.' % chrono)

    simulation_info = xSimulationInfo()
    simulation_info.gti_list = gti_list
    simulation_info.roi_model = ROI_MODEL
    simulation_info.irf_name = kwargs['irfname']
    simulation_info.aeff = aeff
    simulation_info.psf = psf
    simulation_info.modf = modf
    simulation_info.edisp = edisp
    event_list.sort()
    event_list.write_fits(kwargs['outfile'], simulation_info)

    logger.info('All done %s!' % chrono)
    return kwargs['outfile']

if __name__=='__main__':
    args = PARSER.parse_args()
    startmsg()
    chandra2ximpol(args.evfile, **args.__dict__)
