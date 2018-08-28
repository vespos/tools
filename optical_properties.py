# -*- coding: utf-8 -*-
"""
Created on Tue Aug 21 11:09:30 2018

@author: esposito_v
"""

import numpy as np
import matplotlib.pyplot as plt

#plt.close('all')

""" CONSTANTS """
e_0 = 8.85E-12; # [s/(Ohm*m)]
c = 299792458; # [m/s]

""" MATERIAL AND LIGHT DATA """
sigma = 500 # [1/(Ohm*cm)]
sigma = sigma*100 # [1/(Ohm*m)]
e_r = 0.5 # [] real part of epsilon
wvl = 800E-9 # [m] wavelength
f = c/wvl # [1/s]
omega = 2*np.pi*f # [1/s]





n = np.sqrt(e_r+1j*sigma/e_0/omega)
alpha = 2*n.imag*omega/c

str = '\npenetration depth 1/alpha = {:0.2f} [nm]'.format((1/alpha)*10**9)
print(str)

""" REFLECTIVITY """
nr = n.real;
ni = n.imag;
thetai = np.deg2rad(np.arange(0,91,1))
thetat = np.arcsin( np.sin(thetai)/n )

r = (np.cos(thetai)-n*np.cos(thetat)) / (np.cos(thetai)+n*np.cos(thetat))
R = np.absolute( r )**2 # reflectivity at incident angle theta1 (s-pol)
R_N = ( 1-n / (1+n) )**2 # normal incidence reflectivity

""" p-polarized light """
num = -( (nr**2-ni**2+1j*(2*nr*ni)) )*np.cos(thetai) + np.sqrt( (nr**2-ni**2-np.sin(thetai)**2+1j*2*nr*ni) )
den = ( (nr**2-ni**2+1j*(2*nr*ni)) )*np.cos(thetai) + np.sqrt( (nr**2-ni**2-np.sin(thetai)**2+1j*2*nr*ni) )
R_p = np.absolute(num / den)**2

""" s-polarized light """
num = np.cos(thetai) - np.sqrt( (nr**2-ni**2-np.sin(thetai)**2+1j*2*nr*ni) )
den = np.cos(thetai) + np.sqrt( (nr**2-ni**2-np.sin(thetai)**2+1j*2*nr*ni) )
R_s = np.absolute(num / den)**2

thetai = np.rad2deg(thetai)
thetat = np.rad2deg(thetat.real)
plt.figure('Reflectivity')
plt.title('Reflectivity')
plt.plot(thetai, R_p, label='p-polarized')
plt.xlabel('incident angle')
plt.ylabel('reflectance')
plt.plot(thetai, R_s, label='s-polarized')
plt.xlabel('incident angle')
plt.ylabel('reflectance')
plt.legend(loc='upper left')

plt.figure('Refraction')
plt.title('Refraction')
plt.plot(thetai, thetat)
plt.xlabel('theta_i')
plt.ylabel('theta_t')
