#!/usr/bin/env python
import scipy as sp
from astropy.io import fits

class xEventList():
    def __init__(self):
        self.energy_array = sp.array([],dtype='double')
        self.time_array   = sp.array([],dtype='double')
        self.ra_array     = sp.array([],dtype='double')
        self.dec_array    = sp.array([],dtype='double')
        self.x_array      = sp.array([],dtype='int')
        self.y_array      = sp.array([],dtype='int')
        self.l_array      = sp.array([],dtype='double')
        self.b_array      = sp.array([],dtype='double')
        self.angle_array  = sp.array([],dtype='double')
        pass
    
    def fill(self,event):
        self.energy_array=sp.append(self.energy_array,event.energy)
        self.time_array=sp.append(self.time_array,event.time)
        self.ra_array=sp.append(self.ra_array,event.ra)
        self.dec_array=sp.append(self.dec_array,event.dec)
        self.x_array=sp.append(self.x_array,event.x)
        self.y_array=sp.append(self.y_array,event.y)
        self.l_array=sp.append(self.l_array,event.l)
        self.b_array=sp.append(self.b_array,event.b)
        self.angle_array=sp.append(self.angle_array,event.angle)
        pass
    
    def write_fits(self,file_path):       
        col_time = fits.Column(name='time', format='E', array=self.time_array)
        col_energy = fits.Column(name='energy', format='E', array=self.energy_array)
        col_ra = fits.Column(name='ra', format='E', array=self.ra_array)
        col_dec = fits.Column(name='dec', format='E', array=self.dec_array)
        col_l = fits.Column(name='l', format='E', array=self.l_array)
        col_b = fits.Column(name='b', format='E', array=self.b_array)
        col_x = fits.Column(name='x', format='D', array=self.x_array)
        col_y = fits.Column(name='y', format='D', array=self.y_array)
        col_angle = fits.Column(name='angle', format='E', array=self.angle_array)
        cols = fits.ColDefs([col_time,col_energy,col_ra,col_dec,col_l,col_b,col_x,col_y,col_angle])
        tbhdu = fits.BinTableHDU.from_columns(cols)
        tbhdu.writeto(file_path,clobber=True)    
