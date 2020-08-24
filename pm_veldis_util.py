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
