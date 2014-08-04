"""
pyparticleprobe.dsd_calcs.params
===============================

A grouping of functions to calculate parameters of a drop size distribution.

Adapted by Nick Guy.

"""
#*****************************
#  params.py
#*****************************
# DESCRIPTION::
#   This is a grouping of scripts designed to process DSD data of any type.
#   The impetus for the creation was NOAA P-3 2-D optical precipitaton 
#   imaging probe (PIP) data.
# HISTORY::
#   8 Jan 2014 - Nick Guy. NRC, NOAA/NSSL (nick.guy@noaa.gov)
#                Converted NCL functions below to Python 
# HISTORY::
#  14 May 2014 - Updated routine structure NOAA/NSSL/WRDD, NRC
# NOTES::
#
# FUNCTIONS::
# refl_h - Horizontal equivalent radar reflectivity
# refl_v - Vertical equivalent radar reflectivity
# diff_refl - Differential reflectivity
# spec_diff_phase - Specific differential phase
#-------------------------------------------------------------------
# Load the needed packages
import numpy as np
from muphyspy.dsd.moments import *
from muphyspy.dsd.params import *

from pytmatrix.tmatrix import Scatterer
from pytmatrix.psd import PSDIntegrator, GammaPSD
from pytmatrix import orientation, radar, tmatrix_aux, refractive
#-------------------------------------------------------------------
# Define various constants that may be used for calculations
rho0=1.2 # reference density of air at the surface [kg m^-3]
rhoL=1000. # density of water [kg m^-3]
#
#===============================================================
# BEGIN FUNCTIONS
#===============================================================
def refl_h(Nd,Diam,sum=False,axis=None):
    """Horizontal equivalent reflectivity factor computed given drop size distribution 
 INPUT::
  Nd              = Drop concentraton as a function of drop size [m^-3]
  Diam            = Drop size diameter [in mm]
       OPTIONAL
  sum             = True results in a summed variable additionally returned
  axis            = The axis to sum over is sum set to True
 OUTPUT::
  Ze_h              = Reflectivity [mm^6 * m^-3]
  Ze_h_sum          = Optional returned array, Ze summed along an axis,masked = 0
 USAGE::
  Ze_h = refl_h(Nd,Diam)
 NOTES::
  Input variables must have same dimensionality

  Calculation comes from integrating (summing) over the 
  Ze = Summation(Integral)[i=0,n] N(Di) * Di^6

      n
     ---
     \              
     /   N(Di) * Di^6
     ---
     i=0
     
  As in Brandes et al. (1995) "A study of thunderstorm microphysics with multiparameter
  radar and aircraft observations"
    """
#---------------------------------------
# Calculate horizontal equivalent reflectivity
    Ze_h = Nd * Diam**6.

    if sum:
        Ze_h_sum = np.sum(Ze_h,axis=axis)
        Ze_h_sum = np.ma.masked_equal(Ze_h_sum,0)
        
        return Ze_h,Ze_h_sum
    else:
        return Ze_h
#**====================================================

def refl_v(Nd,Diam,sum=False,axis=None):
    """Vertical equivalent reflectivity factor computed given drop size distribution 
 INPUT::
  Nd              = Drop concentraton as a function of drop size [m^-3]
  Diam            = Drop size diameter [in mm]
       OPTIONAL
  sum             = True results in a summed variable additionally returned
  axis            = The axis to sum over is sum set to True
 OUTPUT::
  Ze_v              = Reflectivity [mm^6 * m^-3]
  Ze_v_sum          = Optional returned array, Ze summed along an axis,masked = 0
 USAGE::
  Ze_v = refl_v(Nd,Diam)
 NOTES::
  Input variables must have same dimensionality

  Calculation comes from integrating (summing) over the 
  Ze = Summation(Integral)[i=0,n] N(Di) * Di^6

      n
     ---
     \              
     /   N(Di) * Di^6 * r^(7/3)
     ---
     i=0
     
  As in Brandes et al. (1995) "A study of thunderstorm microphysics with multiparameter
  radar and aircraft observations"
  The vertical polarization value is reduced using the axis ratio relationship computed
  for wind-tunnel data of Pruppacher and Beard (1970).
    """
#---------------------------------------
    # Calculate vertical equivalent reflectivity by reducing via by axis ratio
    Ze_v = Nd * Diam**6. * (1.03 - 0.062*Diam)**(7./3.)

    if sum:
        Ze_v_sum = np.sum(Ze_v,axis=axis)
        Ze_v_sum = np.ma.masked_equal(Ze_v_sum,0)
        
        return Ze_v,Ze_v_sum
    else:
        return Ze_v
#**====================================================

def diff_refl(Nd,Diam,sum=False,axis=None):
    """Vertical equivalent reflectivity factor computed given drop size distribution 
 INPUT::
  Nd              = Drop concentraton as a function of drop size [m^-3]
  Diam            = Drop size diameter [in mm]
       OPTIONAL
  sum             = True results in a summed variable additionally returned
  axis            = The axis to sum over is sum set to True
 OUTPUT::
  ZDR              = Differential Reflectivity [mm^6 * m^-3]
  ZDR_sum          = Optional returned array, Ze summed along an axis,masked = 0
 USAGE::
  ZDR = refl(Nd,Diam)
 NOTES::
  Input variables must have same dimensionality

  Calculation comes from the ratio of horizontal reflectivity over vertical reflectivity
  
            Zh
     ZDR = ----
            Zv
    """
#---------------------------------------
    # Calculate horizontal and vertical equivalent reflectivity
    Ze_h = refl_v(Nd,Diam,sum=sum,axis=axis)
    Ze_v = refl_v(Nd,Diam,sum=sum,axis=axis)
    
    ZDR = Ze_h/Ze_v

    if sum:
        Ze_v_sum = np.sum(Ze_v,axis=axis)
        Ze_v_sum = np.ma.masked_equal(Ze_v_sum,0)
        
        return ZDR,ZDR_sum
    else:
        return ZDR
#**====================================================

def spec_diff_phase(Nd,Diam,Dm,sum=False,axis=None):
    """Differential reflectivity factor computed given drop size distribution 
 INPUT::
  Nd              = Drop concentraton as a function of drop size [m^-3]
  Diam            = Drop size diameter [mm]
  Dm              = Mass-weighted mean diameter [mm]
  Lam             = Wavelength of radar [m]
       OPTIONAL
  sum             = True results in a summed variable additionally returned
  axis            = The axis to sum over is sum set to True
 OUTPUT::
  KDP             = Specific differential phase [deg km^-1]
  KDP_sum          = Optional returned array, KDP summed along an axis,masked = 0
 USAGE::
  KDP = spec_diff_phase(Nd,Diam,Dm,**[args])
 NOTES::
  Calculation as in Bringi and Chandrasekar (2001) "Polarimetric Doppler Weather Radar" 
  KDP = (180/Lam) * 10^-3 * C * LWC * (0.062 * Dm) [Eq. 7.18]
  
  Uses the equilibrium drop shape from the linear fit model of Pruppacher and Beard (1970):
  b/a = 1.03 - 0.062D [Eq. 7.2 in BC01]  range 1 <= D <= 9 mm

    """
#---------------------------------------
    # Set the constant C
    C = 3.75
    
    # Check the sum and axis
    if axis is None:
        sumaxis=0
    else:
        sumaxis=axis

    # Compute the water mass and masked summed time series (LWC) (masked water = 0)
    Mwater, LWC = watermass(Nd,Diam,sum=True,axis=sumaxis)
    
    # Calculate specific differential phase
    KDP = (180./Lam) * 10**-3 * C * LWC * (0.062 * Dm)

    if sum:
        KDP_sum = np.sum(KDP,axis=axis)
        KDP_sum = np.ma.masked_equal(KDP_sum,0)
        
        return KDP,KDP_sum
    else:
        return KDP
#**====================================================