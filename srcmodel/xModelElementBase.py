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



from ximpol.__logging__ import abort



class ModelElementKeyMissing(Exception):
    pass



class ModelElementKeyUnknown(Exception):
    pass



class ModelElementKeyCastError(Exception):
    pass



class xModelElementBase:

    """
    """

    REQUIRED_KEYS = []
    OPTIONAL_KEYS = []
    TYPE_DICT = {}

    def __init__(self, **kwargs):
        """
        """
        for key in self.REQUIRED_KEYS:
            if not kwargs.has_key(key):
                raise ModelElementKeyMissing('key "%s" is required by %s.' %\
                                             (key, self.__class__.__name__))
        for key in kwargs.keys():
            if key not in self.REQUIRED_KEYS + self.OPTIONAL_KEYS:
                raise ModelElementKeyUnknown('unknown key "%s" for %s.' %\
                                             (key, self.__class__.__name__))
        for key, typeCast in self.TYPE_DICT.items():
            if kwargs.has_key(key):
                try:
                    kwargs[key] = typeCast(kwargs[key])
                except Exception as e:
                    raise ModelElementKeyCastError('%s in %s' %\
                                                   (e, self.__class__.__name__))
        self.__Kwargs = kwargs

    def __str__(self):
        """ Default string formatting.
        """
        return str(self.__Kwargs)

    def __getattr__(self, name):
        """
        """
        return self.__Kwargs[name]



def test():
    """ Test code.
    """
    pass



if __name__ == '__main__':
    test()
