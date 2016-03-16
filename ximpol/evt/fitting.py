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


import xspec

from ximpol.evt.binning import xBinnedCountSpectrum
from ximpol.irf import irf_file_path


class xSpectralFitter:

    """Interface to XSPEC spectral fitting.
    """

    def __init__(self, file_path, **kwargs):
        """
        """
        xspec.AllData.clear()
        self.count_spectrum = xBinnedCountSpectrum(file_path)
        irf_name = self.count_spectrum.primary_header_keyword('IRFNAME')
        self.spectrum = xspec.Spectrum(file_path)
        self.spectrum.response = irf_file_path(irf_name, 'rmf')
        self.spectrum.response.arf = irf_file_path(irf_name, 'arf')
        self.spectrum.ignore('**-%f' % kwargs['emin'])
        self.spectrum.ignore('%f-**' % kwargs['emax'])
        self.model = xspec.Model(kwargs['model'])
        xspec.Fit.perform()

    def fit_parameter(self, i):
        """Return the best-fit value of a fit parameter, along with the
        associated error.
        """
        return self.model(i).values[0], self.model(i).sigma

    def fit_parameters(self):
        """Return the fit parameters and the associated errors.
        """
        return [self.fit_parameter(i) for i in \
                range(1, self.model.nParameters + 1)]

    def plot(self):
        """Plot the fit.
        """
        xspec.Fit.show()
        xspec.Plot.device = '/xs'
        xspec.Plot.xAxis = 'keV'
        xspec.Plot('ldata', 'resid')
