"""
pyparticleprobe.dsd_calcs.exp_dsd 
===============================

A grouping of functions to calculate parameters of a drop size distribution
 using an exponential model distribution.  Originally proposed by 
 Marshall and Palmer (1948).

Per Walvogel (1974).
      N(D) = N0 * exp(-Lamdba * D) (0 <= D <= Dmax)

 where N0 is the intercept parameter (m^-3 cm^-1) and Lambda is
 the slope parameter (cm^-1).

Adapted by Nick Guy.

"""
# HISTORY::
#  28 Jan 2014 - Nick Guy.   NRC, NOAA/NSSL (nick.guy@noaa.gov)   
#                Adapted from scripts coded by Terry Schuur (OU CIMMS)
#-------------------------------------------------------------------

def intercept(M3,M6):
    """Compute the intercept parameter (N0) using the method of moments as in 
  Waldvogel (1974)
 INPUT::
  M3              = Third moment of the DSD [unitless]
  M6              = Sixth moment of the DSD [m^3]
 OUTPUT::
  N0              = Intercept [m^-1]
 USAGE::
  N0 = intercept(M3,M6)
 NOTES::
  This method usd the 3rd and 6th moments to calculate the intercept.
    """
#---------------------------------------
    N0 = 446.437 * ((M3/M6)**(4./3.)) * M3

    return N0
#====================================================

def slope(M3,M6):
    """Compute the slope (Lambda) using the method of moments as in 
  Waldvogel (1974)
 INPUT::
  M3              = Third moment of the DSD [unitless]
  M6              = Sixth moment of the DSD [m^3]
 OUTPUT::
  Lambda              = Slope [m^-1]
 USAGE::
  Lambda = slope(M3,M6)
 NOTES::
  This method usd the 3rd and 6th moments to calculate the slope.
    """
#---------------------------------------
    Lambda = 6.12 * ((M3/M6)**(1./3.))

    return Lambda
#====================================================
