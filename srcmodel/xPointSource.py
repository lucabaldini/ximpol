#!/usr/bin/env python
# *********************************************************************
# * Copyright (C) 2015 Luca Baldini (luca.baldini@pi.infn.it)         *
# *                                                                   *
# * For the license terms see the file LICENSE, distributed           *
# * along with this software.                                         *
# *********************************************************************
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



from ximpol.srcmodel.xModelElementBase import xModelElementBase
from ximpol.srcmodel.xModelParameter import xModelParameter
from ximpol.srcmodel.xSpectralComponent import xSpectralComponent



class xPointSource(xModelElementBase):

    """
    """

    REQUIRED_KEYS = ['RA', 'Dec', 'spectrum']

    def __init__(self, name, **kwargs):
        """
        """
        xModelElementBase.__init__(self, name, **kwargs)
        self.RA = xModelParameter('RA', **self.RA)
        self.Dec = xModelParameter('Dec', **self.Dec)
        for key, value in self.spectrum.items():
            self.spectrum[key] = xSpectralComponent(key, **value)

    def spectralComponents(self):
        """ Return the (unordered) list of spectral components.
        """
        return self.spectrum.values()

    def spectralComponent(self, name = 'main'):
        """ Return a given spectral component.

        By default returns the component named "main".
        """
        return self.spectrum[name]

    def __str__(self):
        """ String formatting.
        """
        _str = 'Point source %s at (%s, %s)' % (self.name(), self.RA, self.Dec)
        for comp in self.spectralComponents():
            _str += '\n- %s' % comp
        return _str




def test():
    """ Test code.
    """
    import os
    import yaml
    from ximpol.__package__ import XIMPOL_SRCMODEL
    filePath = os.path.join(XIMPOL_SRCMODEL, 'yaml', 'simple_source.yaml')
    tree = yaml.load(open(filePath))
    srcName = 'test source'
    srcType = 'point source'
    s = xPointSource(srcName, **tree[srcName][srcType])
    print(s)



if __name__ == '__main__':
    test()
