## This file contains some functions required to fit galaxy spectra with templates to extract the velocity dispersion 
## information of the galaxy. 

## Author : Pritom Mozumdar

#from ppxf.ppxf import ppxf
#import ppxf.ppxf_util as util
#from specim.specfuncs import spec1d
import numpy as np
#import matplotlib.pyplot as plt
#import glob
#from random import sample
#import pandas as pd
#import seaborn as sn
#from collections import Counter

def velocity_scale(lamda_gal):
    '''
    This function calculates and returns the associated velocity scale of the
    galaxy.
    
    Parameters
    ---------------
    lamda_gal: array
        An array containing the wavelengths of the galaxy spectra. 
    
    Returns
    -------------
    vel_scale: float
        Velocity scale of the galaxy.
    '''
        
    c = 299792.458                             # speed of light in km/s
    frac_lamda = lamda_gal[1] / lamda_gal[0]   # Constant lambda fraction per pixel
    vel_scale =  np.log(frac_lamda)*c         # velocity scale in km/s per pixel
    print('Velocity scale = %f km/s' %vel_scale)
    
    return vel_scale

###############################################################################################

def gen_sigma_diff(sig_ins=0, fwhm_temp=0, lam_gal=0, lam_temp=0):
    '''
    This function calculates and returns the differences in sigma per wavelength of 
    the instrumental LSF's used to collect galaxy spectrum and template spectra.
    
    Parameters
    ---------------
    sig_ins: float
        sigma value of the instrumental LSF used to collect galaxy spectra.
    
    fwhm_temp: float
        FWHM value of the template spectra.
        
    lam_gal: array
       An array containing the wavelengths of the galaxy spectra.
       
    lam_temp: array
       An array containing the wavelengths of the template spectra.
    
    Returns
    -------------
    sigma_diff: array
        An array containing the differences in sigma per wavelength.
     
    '''
    
    sigma_instrument = sig_ins                       #sigma of the instrumental LSF
    fwhm_galaxy = 2.355 * sigma_instrument           # FWHM of every pixel in Angstrom
    fwhm_galaxy_spectra  = np.full(len(lam_gal), fwhm_galaxy)
    
    fwhm_interp_gal_spec = np.interp(lam_temp, lam_gal, fwhm_galaxy_spectra)  #interpolated fwhm
    fwhm_diff = np.sqrt(fwhm_interp_gal_spec**2 - fwhm_temp**2)
    sigma_diff = fwhm_diff / 2.355
    
    return sigma_diff