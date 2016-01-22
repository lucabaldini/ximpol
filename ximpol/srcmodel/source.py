#!/urs/bin/env python
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


import numpy

from ximpol.srcmodel.img import xFitsImage
from ximpol.srcmodel.spectrum import xCountSpectrum
from ximpol.evt.event import xMonteCarloEventList

from collections import OrderedDict


class xSourceBase:

    """
        Base class for the source object
    """

    def __init__(self, name):
        """
        """
        self.name = name
        self.identifier = None
        
    def rvs_sky_coordinates(self, size=1):
        """
        """
        pass
    
    def __str__(self):
        return '%s (id:%d)' %(self.name, self.identifier)
    
    def rvs_event_list(self, aeff, psf, modf, edisp, sampling_time):
        count_spectrum = xCountSpectrum(self.spectrum, aeff, sampling_time)
        num_events = numpy.random.poisson(count_spectrum.light_curve.norm())
        col_time = count_spectrum.light_curve.rvs(num_events)
        col_time.sort()
        event_list = xMonteCarloEventList()
        event_list.set_column('TIME', col_time)
        col_mc_energy = count_spectrum.rvs(col_time)
    
        event_list.set_column('MC_ENERGY', col_mc_energy)
        col_pha = edisp.matrix.rvs(col_mc_energy)
        event_list.set_column('PHA', col_pha)
        event_list.set_column('ENERGY', edisp.ebounds(col_pha))
        col_mc_ra, col_mc_dec = self.rvs_sky_coordinates(num_events)
        event_list.set_column('MC_RA', col_mc_ra)
        event_list.set_column('MC_DEC', col_mc_dec)
        col_ra, col_dec = psf.smear(col_mc_ra, col_mc_dec)
        event_list.set_column('RA', col_ra)
        event_list.set_column('DEC', col_dec)
        polarization_degree = self.polarization_degree(col_mc_energy, col_time)
        polarization_angle = self.polarization_angle(col_mc_energy, col_time)
        col_pe_angle = modf.rvs_phi(col_mc_energy, polarization_degree,
                                    polarization_angle)
        event_list.set_column('PE_ANGLE', col_pe_angle)
        event_list.set_column('MC_SRC_ID', self.identifier)
        return event_list

    
class xPointSource(xSourceBase):

    """
    """

    def __init__(self, name, ra, dec):
        """
        """
        xSourceBase.__init__(self, name)
        self.ra = ra
        self.dec = dec

    def rvs_sky_coordinates(self, size=1):
        """
        """
        _ra = numpy.zeros(size)
        _ra.fill(self.ra)
        _dec = numpy.zeros(size)
        _dec.fill(self.dec)
        return (_ra, _dec)


class xExtendedSource(xSourceBase):

    """
    """

    def __init__(self, name, file_path):
        """
        """
        xSourceBase.__init__(self, name)
        self.image = xFitsImage(file_path)

    def rvs_sky_coordinates(self, size=1):
        """
        """
        return self.image.rvs_coordinates(size)

    
class xROIModel(OrderedDict):
    '''
    '''

    def __init__(self, roi_ra, roi_dec, roi_radius):
        '''
        '''
        OrderedDict.__init__(self)
        self.ra=roi_ra
        self.dec=roi_dec
        self.radius=roi_radius
        pass

    def add_source(self, source):
        source.identifier=len(self)
        self[source.name]=source
        pass

    def __str__(self):
        txt=""
        for s in self.values():
            txt+="%s\n" % s
            pass
        return txt

    def rvs_event_list(self, aeff, psf, modf, edisp, sampling_time):
        event_list = xMonteCarloEventList()
        for source in self.values():
            event_list += source.rvs_event_list(aeff, psf, modf, edisp, sampling_time)
            pass
        event_list.sort()
        return event_list
    pass
