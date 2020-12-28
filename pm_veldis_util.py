## This file contains some functions required to fit galaxy spectra with templates to extract the velocity dispersion 
## information of the galaxy. 

## Author : Pritom Mozumdar

#from ppxf.ppxf import ppxf
import ppxf.ppxf_util as util
from specim.specfuncs import spec1d
import numpy as np
import matplotlib.pyplot as plt
import glob
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

def wav_dev(lamda_gal):
    '''
    This function calculates the parameter 'dv' to account for the difference
    of the initial wavelength in the galaxy and template spectra.
    
    Parameters
    ---------------
    lamda_gal: float
        Starting wavelength of the galaxy spectra. 
    
    Returns
    -------------
    dv: float
        The parameter to account for the initial wavelegth difference.
    '''
    
    c = 299792.458               # speed of light in km/s
    lamda_temp = 3465.00         # starting wavelength of the templates
                                 # in the Indo-US library.
    dv = c*np.log(lamda_temp / lamda_gal) 
    print('dv = %f ' %dv)
    
    return dv

###############################################################################################

def gen_sigma_diff(sig_ins=0, fwhm_temp=None, lam_gal=0, lam_temp=0):
    '''
    This function calculates and returns the differences in sigma per wavelength of 
    the instrumental LSF's used to collect galaxy spectrum and template spectra.
    
    Parameters
    ---------------
    infile: string 
        contains a template file location and/or file name.
        
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
    try:
        if sig_ins > 0 :
            sigma_instrument = sig_ins  #sigma of the instrumental LSF
    except:
        print("Error : sigma of the instrumental LSF should be greater than zero.")
        
    if fwhm_temp is None:
        print("\nAs no \'fwhm_temp\' value is provided, FWHM for the Indo-US template"\
              "library will be used as default value")
        fwhm_temp = 1.35                             # for indo-us template library
    else:
        fwhm_temp = fwhm_temp
    
    lam_temp = spec1d.Spec1d('../TEXT/101484.txt', verbose=False)['wav']    
    fwhm_galaxy = 2.355 * sigma_instrument           # FWHM of every pixel in Angstrom
    fwhm_galaxy_spectra  = np.full(len(lam_gal), fwhm_galaxy)
    
    fwhm_interp_gal_spec = np.interp(lam_temp, lam_gal, fwhm_galaxy_spectra)  #interpolated fwhm
    fwhm_diff = np.sqrt(fwhm_interp_gal_spec**2 - fwhm_temp**2)
    sigma_diff = fwhm_diff / 2.355
    
    plt.plot(lam_temp, sigma_diff,'.', label='sigma_diff')
    plt.legend()
    plt.show()
    
    return sigma_diff

####################################################################################################

def masking(pixel_range, log_lamda_gal):
    '''
    This function generate and returns a boolean array with value 'False'
    in the pixel locations which should be excluded from the fit.
    
    Parameters
    ---------------   
    pixel_range: list
        A list of tuples where each tuple contains start and end of the
        pixel range needs to be excluded.
        
    log_lamda_gal: array
        This array contains the values of the logarithmically 
        rebinned wavelengths.
    
    Returns
    -------------
    mask : boolean array
        Boolean array with with value 'False' in the pixel locations 
        which should be excluded from the fit.
        
    '''
    
    mask = np.zeros(len(log_lamda_gal), dtype=bool)
    for i,p in enumerate(pixel_range):
        mask |= (log_lamda_gal>=p[0]) & (log_lamda_gal <= p[1])
    return (~mask)

######################################################################################################

def gen_rebinned_templates(lib_path=0, temp_num=0, temp_array=0, sigma_diff=0, v=0):
    '''
    This function generates and returns an array containing logarithmically 
    rebinned template spectra.
    
    Parameters
    ---------------
    lib_path: string
        path to the directory containing template library.
        
    temp_num: int
        Number of templates that would be logarithmically rebinned. If
        given that amount of template spectra would be fetched from
        library, rebinned and stored in the array.
    
    temp_array: array
        An array containing template file names which would be 
        logarithmically rebinned. If given only those template spectra 
        would be fetched from library, rebinned and stored in the array.
    
    sigma_diff: array
        An array with the differences in sigma per wavelength of 
        the instrumental LSF's used to collect galaxy spectrum and 
        template spectra.
    
    v: float
        Velocity scale of the galaxy.
    
    Returns
    -------------
    templates: array
        An array containging all the logarithmically rebinned and normalized
        template spectra.
    
    '''
    
    filename = []
    templates = []
    
    if(temp_num):
        
        indo_us_library = glob.glob(lib_path)[:temp_num]
        
        lam_temp = spec1d.Spec1d(indo_us_library[0], verbose=False)['wav']
        lam_temp_range = [lam_temp[0], lam_temp[-1]]

        for j, file_name in enumerate(indo_us_library):
        
            filename.append(file_name)
        
            template_data = spec1d.Spec1d(file_name, verbose=False)
            template_spectra = template_data['flux']

            # perform convolution with variable sigma_diff    
            convolved_temp_spectra= util.gaussian_filter1d(template_spectra, sigma_diff)  

            template_spectra_rebinned = util.log_rebin(lam_temp_range, convolved_temp_spectra, 
                                                       velscale=v)[0]
            nor_temp = template_spectra_rebinned / np.median(template_spectra_rebinned)  
            
            templates.append(nor_temp)
            
        templates = np.swapaxes(np.array(templates), 0, 1)
            
        return templates
    
    else:
        
        lam_temp = spec1d.Spec1d(temp_array[0], verbose=False)['wav']
        lam_temp_range = [lam_temp[0], lam_temp[-1]]
        
        for j, file_name in enumerate(temp_array):
        
            filename.append(file_name)
        
            template_data = spec1d.Spec1d(file_name, verbose=False)
            template_spectra = template_data['flux']

            # perform convolution with variable sigma_diff    
            convolved_temp_spectra= util.gaussian_filter1d(template_spectra, sigma_diff)  

            template_spectra_rebinned = util.log_rebin(lam_temp_range, convolved_temp_spectra, 
                                          velscale=v)[0]
            nor_temp = template_spectra_rebinned / np.median(template_spectra_rebinned)
            
            templates.append(nor_temp)
          
        templates = np.swapaxes(np.array(templates), 0, 1)
            
        return templates