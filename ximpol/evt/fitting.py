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
import numpy

from ximpol.evt.binning import xBinnedCountSpectrum
from ximpol.irf import irf_file_path
from ximpol.utils.logging_ import logger


class xSpectralFitter:

    """Interface to XSPEC spectral fitting.
    """

    def __init__(self, file_path, **kwargs):
        """
        """
        xspec.AllData.clear()
        self.count_spectrum = xBinnedCountSpectrum(file_path)
        self.spectrum = xspec.Spectrum(file_path)
        try:
            self.spectrum.response
        except Exception:
            irf_name = self.count_spectrum.primary_header_keyword('IRFNAME')
            self.spectrum.response = irf_file_path(irf_name, 'rmf')
            self.spectrum.response.arf = irf_file_path(irf_name, 'arf')
        self.spectrum.ignore('**-%f' % kwargs['emin'])
        self.spectrum.ignore('%f-**' % kwargs['emax'])
        self.model = xspec.Model(kwargs['model'])

    def set_parameters(self, param_list, fix=False):
        """Set the values of all the model parameters passed through a list. If
        fix is set to True the parameters are frozen.
        """
        self.model.setPars(param_list)
        if fix:
            for par_index in range(1, len(param_list) + 1):
                self.model(par_index).frozen = True

    def set_parameter(self, par_index, par_value, fix=False):
        """Set the model parameter value corresponding to provided index. If fix
        is set to True the parameter is frozen.
        """
        self.model(par_index).values = par_value
        if fix:
            self.model(par_index).frozen = True

    def fit(self):
        """Perform the fit and display the results.
        """
        xspec.Fit.perform()
        xspec.Fit.show()

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

    def plot(self, arg_list=['ldata', 'delchi'], device='/xs', title=None,
             xaxis='keV', logx=True):
        """Plot the data and the input model. With arg_list parameter it's
        possible to change the displayed quantity (e.g. 'ufspec' to show the
        unfolded spectrum (ph/cm^2/s/keV) or 'ldata' for the folded spectrum
        (counts/s/keV).
        """
        xspec.Plot.device = device
        xspec.Plot.xAxis = xaxis
        xspec.Plot.xLog = logx
        if title is not None:
            xspec.Plot.addCommand('label top '+ title)
        args = ', '.join(x for x in arg_list)
        xspec.Plot(args)

    def save_plot(self, outfile, arg_list=[], device='/cps', title=None,
                  xaxis='keV', logx=True):
        """Save the plot to the outfile directory. The default format is '.ps'.
        """
        xspec.Plot.device = outfile + device
        xspec.Plot.xAxis = xaxis
        xspec.Plot.xLog = logx
        if title is not None:
            xspec.Plot.addCommand('label top '+ title)
        args = ', '.join(x for x in arg_list)
        logger.info('Saving the plot to %s...' %outfile)
        xspec.Plot(args)

    def save_spectrum(self, outfile, unfolded=True):
        """Save the data and the model to file. With the unfolded keyword is
        possible to choose to save the unfolded or foleded flux values.
        """
        xspec.Plot.device = '/null'
        xspec.Plot.xAxis = 'keV'
        if unfolded:
            xspec.Plot('ufspec') #ph/cm^2/s/keV
        else:
            xspec.Plot('data')   #counts/s/keV
        x = xspec.Plot.x()
        y = xspec.Plot.y()
        y_err = xspec.Plot.yErr()
        mod_vals = xspec.Plot.model()
        logger.info('Saving the spectrum to %s...' %outfile)
        numpy.savetxt(outfile, numpy.c_[x,y ,y_err,mod_vals])
