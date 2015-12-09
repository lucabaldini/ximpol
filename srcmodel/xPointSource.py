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



class xPointSource(xModelElementBase):

    """
    """

    REQUIRED_KEYS = ['RA', 'Dec', 'spectrum']

    def __init__(self, **kwargs):
        """
        """
        xModelElementBase.__init__(self, **kwargs)
        self.RA = xModelParameter(**self.RA)
        self.Dec = xModelParameter(**self.Dec)

    def __str__(self):
        """ String formatting.
        """
        return 'Point source at (RA = %s, Dec = %s)' % (self.RA, self.Dec)



def test():
    """ Test code.
    """
    kwargs = {'RA': {'value': 23.63875, 'unit': 'deg'},
              'Dec': {'value': 51.736298, 'unit': 'deg'},
              'spectrum': 'test'}
    s = xPointSource(**kwargs)
    print(s)

    

if __name__ == '__main__':
    test()
