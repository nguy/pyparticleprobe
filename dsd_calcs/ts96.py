"""
plotparticleprobe.dsd_calcs.ts96
===============================

A grouping of functions to calculate parameters of a drop size distribution
using the methodology of Tokay and Short (1996). 
That study assumes a gamma distributionof the DSD 
      N(D) = N0 * D^m * exp(-Lambda*D)

 and uses the method of moments, with the xth moment expressed:
                 Gamma(m + x + 1)
      Mx = N0 * -----------------
                Lambda^(m + x + 1)

 where N0 is the intercept parameter (mm^(-1-m) m^-3), 
  Lambda is the slope parameter (mm^-1), and m is the shape 
  parameter (unitless) following Kozu and Nakamura (1991).
 
            (11G - 8 + [G(G+8)]^(1/2)
       m =  -------------------------
                     2(1-G)
 where
               M4^3
       G = -----------
             M3^2 * M6

 and 
               (m+4)M3
      Lambda =  --------
                  M4

 and 
            Lambda^(m+4)*M4
      N0 = ----------------
              Gamma(M+4)

Adapted by Nick Guy.

"""
# HISTORY::
#   8 Jan 2014 - Nick Guy.   NRC, NOAA/NSSL (nick.guy@noaa.gov)   
#                Converted NCL functions below to Python
# FUNCTIONS::
#  shape - Calculate the shape parameter (mu)
#  slope - Calculate the slope parameter (lambda)
#  intercept - Calculate the intercept parameter (N0)
#  d0 - Calculate the median volume diameter (D0)
#-------------------------------------------------------------------
# Load the needed packages
import numpy as np
import scipy.special as scifunct
#-------------------------------------------------------------------
# Define various constants that may be used for calculations
rho0=1.2 # reference density at the surface [kg m^-3]
rhoL=1000. # density of water [kg m^-3]
#===============================================================
# BEGIN FUNCTIONS
#**===============================================================
def shape(M3,M4,M6):
    """Compute the shape parameter using the method of moments as in 
  Tokay and Short (JAM 1996), Eqn 3
 INPUT::
  M3              = Third moment of the DSD [unitless]
  M4              = Fourth moment of the DSD [m]
  M6              = Sixth moment of the DSD [m^3]
 OUTPUT::
  m               = Shape parameter of gamma DSD model [unitless]
 USAGE::
  m = shape(M3,M4,M6)
 NOTES::
  This particular methodology uses the 3rd, 4th, and 6th moments.
  First the 3rd moment normalized by Dm (M4/M3) is calculated.
  Next the shape paramter is computed.
    """
#---------------------------------------
    # Calculate G parameter
    G = (M4**3)/(M3**2 * M6)
    # Mask values equal to 1 will result in undefined values of m (below)
    G = np.ma.masked_equal(G,1.)

    # Now calculate the shape parameter
    m = (11. * G - 8 + sqrt(G * (G + 8)))/(2. * (1 - G))

    return m
#**====================================================

def slope(M3,M4,m):
    """Compute the slope (Lambda) using the method of moments as in 

  Tokay and Short (JAM 1996), Eqn 5
 INPUT::
  M3              = Third moment of the DSD [unitless]
  M4              = Fourth moment of the DSD [m]
  m               = Shape parameter of gamma DSD model [m]
 OUTPUT::
  Lambda           = Slope [m^-1]
 USAGE::
  Lambda = slope(M3,M4,m)
 NOTES::
  This particular methodology uses the 3rd and 4th moments and shape
  paramter to compute the slope is computed.
    """
#---------------------------------------
    Lambda = (m + 4) * M3/M4

    return Lambda
#====================================================

def intercept(M3,m,Lambda):
    """Compute the intercept parameter (N0) using the method of moments as in 
  Tokay and Short (JAM 1996), Eqn 4
 INPUT::
  M3              = Third moment of the DSD [unitless]
  m               = Shape parameter of gamma DSD model [m]
  Lambda           = Slope [m^-1]
 OUTPUT::
  N0              = Intercept parameter [mm^(-1-shape) m^-3]
 USAGE::
  m = calc_intercept(M3,M4,M6)
    """
#---------------------------------------
    N0 = (Lambda**(m + 4) * M3)/scifunct.gamma(m + 4)

    return N0
#**====================================================

def d0(m,Lambda):
    """Compute the median volume diameter (D0) using the method of moments as in 
  Tokay and Short (JAM 1996) - which assumes a gamma DSD model so that
  the xth moment is expressed:
                 Gamma(m + x + 1)
      Mx = N0 * -----------------
                Lambda^(m + x + 1)
 INPUT::
  m               = Shape parameter of gamma DSD model [m]
  Lambda           = Slope parameter of gamma DSD model [m^-1]
 OUTPUT::
  D0              = Median volume diameter [mm^(-1-shape) m^-3]
 USAGE::
  D0 = d0(M3,M4,M6)
 NOTES::
  This particular methodology uses the 3rd, 4th, and 6th moments.
  First the 3rd moment normalized by Dm (M4/M3) is calculated.
  Next the shape paramter is computed.
    """
#---------------------------------------
    # Calculate D0 [Appears following eq. 9]
    D0 = (3.67 + m)/Lambda

    return D0
#====================================================
