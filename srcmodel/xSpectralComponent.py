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



class xSpectralComponent(xModelElementBase):

    """
    """

    REQUIRED_KEYS = ['shape', 'polarization']
    OPTIONAL_KEYS = []

    def __init__(self, name, **kwargs):
        """
        """
        xModelElementBase.__init__(self, name, **kwargs)
        self.linearPolarization()['angle'] = xModelParameter(
            'polarization angle', **self.linearPolarization()['angle'])
        self.linearPolarization()['degree'] = xModelParameter(
            'polarization degree', **self.linearPolarization()['degree'])

    def linearPolarization(self):
        """
        """
        return self.polarization['linear']

    def linearPolarizationAngle(self):
        """
        """
        return self.linearPolarization()['angle']

    def linearPolarizationDegree(self):
        """
        """
        return self.linearPolarization()['degree']

    def __str__(self):
        """ String formatting.
        """
        _str = 'Spectral shape: %s, %s, %s' %\
               (self.shape, self.linearPolarizationDegree(),
                self.linearPolarizationAngle())
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
    name = 'main'
    c = xSpectralComponent(name, **tree[srcName][srcType]['spectrum'][name])
    print(c)

    

if __name__ == '__main__':
    test()
