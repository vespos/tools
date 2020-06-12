import numpy as np
import matplotlib.pyplot as plt
import time
from scipy.special import erf

""" --------------- Statistics --------------- """
def moments(x,y):
    mean = np.sum(x*y)/np.sum(y)
    variance = np.sum( (x-mean)**2*y )/y.sum()
    sigma = np.sqrt(variance)
    skew = np.sum( (x-mean)**3*y ) /sigma**3
    return mean, sigma, skew


""" --------------- Exponentials and Gaussian --------------- """
def gaussian(x,A,x0,sigma,norm=False):
    if norm:
        return 1/sigma/np.sqrt(2*np.pi)*np.exp(-(x-x0)**2/2/sigma**2)
    return A*np.exp(-(x-x0)**2/2/sigma**2)


def gaussian2d(x,y,A,x0,sigma, angle=0):
    angle = np.deg2rad(angle)
    a = np.cos(angle)**2/2/sigma[0]**2 + np.sin(angle)**2/2/sigma[1]**2
    b = -np.sin(2*angle)/4/sigma[0]**2 + np.sin(2*angle)/4/sigma[1]**2
    c = np.sin(angle)**2/2/sigma[0]**2 + np.cos(angle)**2/2/sigma[1]**2
    return A*np.exp(-( a*(x-x0[0])**2 + 2*b*(x-x0[0])*(y-x0[1]) + c*(y-x0[1])**2 ) )


def step_fct(x,A,x0,FWHM,offset):
    sigma = FWHM/2./np.sqrt(2.*np.log(2))
    y = 0.5*(erf((x-x0)/sigma) + 1)
    y = A*y + offset
    return y


def stretched_exp(t,A,tau,beta):
    return 1-A*np.exp(-(t/tau)**beta)



""" ---------------- Utilities --------------- """
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx


def get_bin_centers(bin_edges):
    return (bin_edges[:-1] + bin_edges[1:]) / 2


# def get_bin_edges(bin_centers):
#     bins = (bin_centers[:-1] + bin_center[1:]) / 2
#     bins = np.append()


def count_occurences(array):
    a = array.flatten()
    values = np.unique(array)
    occurences = [np.sum(a==value) for value in values]
    return values, np.asarray(occurences)
    



""" ---------------- Special --------------- """
def ackermann(m,n):
     if m == 0:
          return (n + 1)
     elif n == 0:
          return ackermann(m - 1, 1)
     else:
          return ackermann(m - 1, ackermann(m, n - 1))