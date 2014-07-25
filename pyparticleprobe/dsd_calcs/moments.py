"""
pyparticleprobe.dsd_calcs.moments
===============================

A grouping of functions to calculate moments of a drop size distribution.

Adapted by Nick Guy.

"""
# HISTORY::
#   8 Jan 2014 - Nick Guy. NRC, NOAA/NSSL (nick.guy@noaa.gov)
#                Converted NCL functions below to Python 
#   9 Mar 2014 - NG: Added option to sum along axis and return additional variable
# NOTES::
#
# FUNCTIONS::
# moment_nth - The nth moment of distribution
# moment_3rd - The 3rd moment of distribution
# moment_4th - The 4th moment of distribution
# moment_6th - The 6th moment of distribution
#-------------------------------------------------------------------
# Load the needed packages
import numpy as np
#-------------------------------------------------------------------
# Define various constants that may be used for calculations
rho0=1.2 # reference density of air at the surface [kg m^-3]
rhoL=1000. # density of water [kg m^-3]
#===============================================================
# BEGIN FUNCTIONS
#===============================================================

def moment_nth(Nd,Diam,X,Uout,sum=False,axis=None):
    """Nth moment of drop size distribution computed
 INPUT::
  Nd              = 2-D Drop concentraton as a function of drop size [m^-3]
  Diam            = 1-D Drop size diameter [in m]
  X               = Moment (integer) to compute
  Uout            = String indicating desired output units
       OPTIONAL
  sum             = True results in a summed variable additionally returned
  axis            = The axis to sum over is sum set to True
 OUTPUT::
  Mx              = X moment of DSD [unit dependent upon Uout]
  Mx_sum          = Optional returned array, Mx summed along an axis,masked = 0
 USAGE::
  Mx = moment_nth(Nd,Diam,X,Uout)
 NOTES::
  Input variables must have same dimensionality
      n
     ---
     \              
     /   N(Di) * Di^M
     ---
     i=0
    """
#---------------------------------------
# Check the units desired
  # Set the unit
    if Uout == 'm':
        Fact=1.
    elif Uout == 'cm':
        Fact=1E2
    elif Uout == 'mm':
        Fact=1E3
    else:
        print('Choose m, cm, or mm units only')
        Fact=np.nan

  # Calculate the third moment
    Mx = (Nd * (Diam * Fact)**X)

    if sum:
        Mx_sum = np.sum(Mx,axis=axis)
        Mx_sum = np.ma.masked_equal(Mx_sum,0)
        
        return Mx,Mx_sum
    else:
        return Mx
#    return Mx
#**====================================================

def moment_3rd(Nd,Diam,sum=False,axis=None):
    """Third moment of drop size distribution computed
 INPUT::
  Nd              = Drop concentraton as a def of drop size [m^-3]
  Diam            = Drop size diameter [in m]
       OPTIONAL
  sum             = True results in a summed variable additionally returned
  axis            = The axis to sum over is sum set to True
 OUTPUT::
  M3              = Third moment of DSD [unitless]
  M3_sum          = Optional returned array, M3 summed along an axis,masked = 0
 USAGE::
  M3 = moment_3rd(Nd,Diam)
 NOTES::
  Input variables must have same dimensionality
  
  Can be used to compute Volume (or Mass)-weighted mean diameter [g/m^3], 
    Also known as mean volume diameter.
    """
#---------------------------------------
# Calculate the third moment
    M3 = (Nd * Diam**3.)

    if sum:
        M3_sum = np.sum(M3,axis=axis)
        M3_sum = np.ma.masked_equal(M3_sum,0)
        
        return M3,M3_sum
    else:
        return M3
#**====================================================

def moment_4th(Nd,Diam,sum=False,axis=None):
    """Fourth moment of drop size distribution computed
 INPUT::
  Nd              = Drop concentraton as a def of drop size [m^-3]
  Diam            = Drop size diameter [in m]
       OPTIONAL
  sum             = True results in a summed variable additionally returned
  axis            = The axis to sum over is sum set to True
 OUTPUT::
  M4              = Fourth moment of DSD [m]
  M4_sum          = Optional returned array, M4 summed along an axis,masked = 0
 USAGE::
  M4 = moment_4th(Nd,Diam)
 NOTES::
  Input variables must have same dimensionality
  
  Can be used to compute Volume (or Mass)-weighted mean diameter [g/m^3], 
    Also known as mean volume diameter.
    """
#---------------------------------------
# Calculate the fourth moment
    M4 = (Nd * Diam**4.)

    if sum:
        M4_sum = np.sum(M4,axis=axis)
        M4_sum = np.ma.masked_equal(M4_sum,0)
        
        return M4,M4_sum
    else:
        return M4
#**====================================================

def moment_6th(Nd,Diam,sum=False,axis=None):
    """Third moment of drop size distribution computed
 INPUT::
  Nd              = Drop concentraton as a def of drop size [m^-3]
  Diam            = Drop size diameter [in m]
       OPTIONAL
  sum             = True results in a summed variable additionally returned
  axis            = The axis to sum over is sum set to True
 OUTPUT::
  M6              = Sixth moment of DSD [m^3]
  M6_sum          = Optional returned array, M6 summed along an axis,masked = 0
 USAGE::
  M6 = moment_6th(Nd,Diam)
 NOTES::
  Input variables must have same dimensionality
  
  Can be used to compute Reflectivity (def above) or the G parameter
   in the method of moments (as discussed in Tokay and Short (JAM# 1996).
    """
#---------------------------------------
# Calculate the sixth moment
    M6 = (Nd * Diam**6.)

    if sum:
        M6_sum = np.sum(M6,axis=axis)
        M6_sum = np.ma.masked_equal(M6_sum,0)
        
        return M6,M6_sum
    else:
        return M6
#====================================================

