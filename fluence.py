# -*- coding: utf-8 -*-
"""
Created on Thu Aug 23 14:30:20 2018

@author: esposito_v
"""

import numpy as np


""" material properties """
""" low T """
z0 = 52 # effective penetration depth [nm]
T = 0.8 # transmission
Cp = 5 # heat capacity [J/K(/mol)]
""" high T """
z0 = 64
T = 0.82
Cp = 10

""" laser properties """
power = 5 # [mW]
rep_rate = 1000 # [Hz]
spotx = 0.046 # FWHM [cm]
spoty = 0.052 # FWHM [cm]
err_spot = 0.001 # [cm]



""" fluence """
fluence = power/rep_rate/spotx/spoty # [mJ/cm^2]
fluence_max = power/rep_rate/(spotx-err_spot)/(spoty-err_spot)
fluence_min = power/rep_rate/(spotx+err_spot)/(spoty+err_spot)
print('\nfluence = {:0.2f} ({:1.2f}, {:2.2f}) mJ/cm^2'.format(fluence, fluence_min-fluence, fluence_max-fluence))

""" energy density """
dlayer = 1 # [nm]
flu_top = fluence*T*np.exp(-0*dlayer/z0)
flu_bottom = fluence*T*np.exp(-1*dlayer/z0)
n0 = (flu_top-flu_bottom)/(dlayer*1E-7)/1000

flu_top = fluence_max*T*np.exp(-0*dlayer/z0)
flu_bottom = fluence_max*T*np.exp(-1*dlayer/z0)
n0_max = (flu_top-flu_bottom)/(dlayer*1E-7)/1000

flu_top = fluence_min*T*np.exp(-0*dlayer/z0)
flu_bottom = fluence_min*T*np.exp(-1*dlayer/z0)
n0_min = (flu_top-flu_bottom)/(dlayer*1E-7)/1000
print('\n n = {:0.2f} ({:1.2f}, {:2.2f}) mJ/cm^2'.format(n0, n0_min-n0, n0_max-n0))



""" Average heating """
heating = power/rep_rate/Cp*T *1000 #[mK]
print('\n heating = {:0.3f} K'.format(heating))