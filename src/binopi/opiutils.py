import oivis.oifits.oifits_read as pif
import oivis.oifits.oifits_transform as oitrans
import oivis.oivisfit.fitfunctions as fitfunctions
import os
import glob
import matplotlib.pyplot as plt
import numpy as np
import binopi.jdutil as jd
from scipy import constants

### Global Parameters ### 
max_spfreq = 70 #Mlambda
max_baseline = 200 #m

########################
### USEFUL FUNCTIONS ###
########################

def teff_to_spi(teff, refwave):
    reffreq = constants.c/refwave
    r = (constants.h*reffreq)/(constants.k*teff)
    return (3 - (r*np.exp(r))/(np.exp(r)-1))

def randInterval(low, up):
    return low + (up - low)*np.random.random()

def randGauss(low, up, mid, spread = 4):
    return np.random.normal(mid, (up-low)/spread)

def importOIFits(fileDirectory, recursive = True):
    """ Imports all .fits files in fileDirectory, 
    using the Oifits class defined in the oivis module. 
    
    Parameters:
    -fileDirectory: Str indicating the path where all the oifits files are saved.
    Can be a parent directory with subfolders.
    
    Returns:
    oifits_read.Oifits object"""

    filelist = [name for name in glob.glob(fileDirectory+'/**/*.fits', recursive=recursive)]
    return pif.Oifits(filelist)

def snr_filter(oifiles, snr_threshold):
    """ Filters the oifiles [Oifits] oidata lower than snr_threshold [Int/Float] 
    
    Parameters:
    -oifiles: OI Data in the Oifits class format.
    
    -snr_threshold: Int/Float, SNR threshold below which to remove observations.
    
    Returns:
    oifits_read.Oifits object"""
    newoifiles = pif.Oifits()
    filtereddata = [oi for oi in oifiles.oidata if (oi['TYPE']=='V2' and oi['VIS2DATA']/oi['VIS2ERR'] >= snr_threshold and oi['VIS2DATA'] <= 1) or (oi['TYPE']=='T3' and abs(oi['T3PHI']) <= 180 and oi['T3PHIERR'] <= 180 and abs(oi['T3PHI']/oi['T3PHIERR']) >= 1)]
    #above returns runtime warning for oifiles with vi2err = 0 (division by 0)
    newoifiles.oidata = filtereddata
    newoifiles.get_v2data()
    newoifiles.get_cpdata()
    return newoifiles

def date_filter(oifiles):
    """ Returns a dict of oifiles corresponding to each individual night.

    Parameters:
    -oifiles: OI Data in the Oifits class format.
     
    Returns:
    dict of oifits_read.Oifits objects"""

    dates = np.unique([int(oi['MJD']) for oi in oifiles.oidata]) #get all individual nights
    filtered_oifiles = {} 
    for date in dates:
        oifile = pif.Oifits()    
        oifile.oidata = [oi for oi in oifiles.oidata if int(oi['MJD']) == date]
        oifile.get_cpdata()
        oifile.get_v2data()
        filtered_oifiles[str(date)] = oifile

    return filtered_oifiles

def getTelescopes(oifile):
    """ Returns a string of the telescope locations used.
    
    Parameters:
    -oifiles: OI Data in the Oifits class format.

    Returns:
    Str
    """
    tel1 = np.unique(oifile.cpdata['TEL1'])
    tel2 = np.unique(oifile.cpdata['TEL2'])
    tel3 = np.unique(oifile.cpdata['TEL3'])
    tel123 = np.unique(np.concatenate((tel1,tel2,tel3)))
    telescopes = ''
    for tel in tel123:
        telescopes += tel
    return telescopes

def telescope_type(oi):
    """ Returns the type of telescope used for the observation ('UT' or 'AT'). 

    Parameters:
    -oi: single observation (element of oifits_read.Oifits.oidata)

    Returns:
    Str
    """
    config = oi['CONFIG']
    if 'U' in config:
        return 'UT'
    else:
        return 'AT'

def telescopetype_filter(oifiles, telescopetype):
    """ Filters the observation based on the specified telescope type ('UT' or 'AT').
    
    Parameters:
    -oifiles: OI Data in the Oifits class format.
    
    -telescopetype: 'UT' or 'AT'
    
    Returns:
    oifits_read.Oifits object, filtered"""    


    newoifiles = pif.Oifits()
    newoifiles.oidata = [oi for oi in oifiles.oidata if telescope_type(oi) == telescopetype]
    newoifiles.get_cpdata()
    newoifiles.get_v2data()
    return newoifiles

def instrument_filter(oifiles, instrument):
    """ Filters the observation based on the specified instrument.
    
    Parameters:
    -oifiles: OI Data in the Oifits class format.
    
    -instrument: name of instrument
    
    Returns:
    oifits_read.Oifits object, filtered"""  
    newoifiles = pif.Oifits()
    newoifiles.oidata = [oi for oi in oifiles.oidata if oi['INSNAME'] == instrument]
    newoifiles.get_cpdata()
    newoifiles.get_v2data()
    return newoifiles
