"""
pyparticleprobe.dsd_calcs.ua98
===============================

A grouping of functions to calculate parameters of a drop size distribution 
 using the methodology of  Ulbrich and Atlas (JAM 1998, JAMC 2007).  
That study assumes a gamma distribution
 (using the complete gamma function) of the DSD 
      N(D) = N0 * D^mu * exp(-Lambda*D)

 and uses the method of moments, with the nth moment expressed:
                               Gamma(n + mu + 1)
      Mn = D^n*N(D)*dD = N0 * -------------------
                              Lambda^(n + mu + 1)

 where N0 is the intercept parameter (mm^(-1-m) m^-3), mu is the shape
  parameter (unitless) and Lambda is the slope paramter (m).
 
            (7-11eta) - [(7-11eta)^2 -4(eta-1)(30eta-12)]^(1/2)
      mu =  -----------------------------------------------------
                                 2(eta-1)
 where
               M4^2
      eta = -----------
             M2 * M6

 and 
               ( (4+mu)(3+mu)M2 )^(1/2)
      Lambda =  (----------------)
               (       M4       )


 and N0 can be found using the 6th moment (M6 or Ze) and solving for N0
           M6 * Lambda^(7+mu)
      N0 = -----------------
              Gamma(7+mu)

Adapted by Nick Guy.

"""
# HISTORY::
#   8 Jan 2014 - Nick Guy. Converted NCL functions below to Python
# FUNCTIONS::
# eta_ratio - Eta ratio
# shape - Shape parameter (mu)
# slope - Slope parameter (Lambda)
# intercept - Intercept parameter (N0)
# d0 - Median volume diameter (D0)
# zr_a - a coefficient (prefactor) in Z-R relationship of form Z = aR^b
# zr_b - b coefficient (exponent) in Z-R relationship of form Z = aR^b
# norm_intercept - Normalized intercept parameter (Nw)
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
#===============================================================

def eta_ratio(M2,M4,M6):
    """Compute the ratio eta using the method of moments

  Ulbrich and Atlas (JAM 1998), Eqn 8
 INPUT::
  M2              = Second moment of the DSD [m^-1]
  M4              = Fourth moment of the DSD [m]
  M6              = Sixth moment of the DSD [m^3]
 OUTPUT::
  Eta             = Ratio of M4^2/(M2*M6) [unitless]
 USAGE::
  m = eta_ratio(M3,M4,M6)
 NOTES::
  This particular methodology uses the 2nd, 4th, and 6th moments.
    """
#---------------------------------------
    Eta = M4**2 / (M2 * M6)

    return Eta
#**====================================================

def shape(M2,M4,M6):
    """Compute the shape parameter mu using the method of moments 

  Ulbrich and Atlas (JAM 1998), Eqns 8 and 9
 INPUT::
  M2              = Second moment of the DSD [m^-1]
  M4              = Fourth moment of the DSD [m]
  M6              = Sixth moment of the DSD [m^3]
 OUTPUT::
  mu              = Shape parameter of gamma DSD model [unitless]
 USAGE::
  mu = shape(M2,M4,M6)
 NOTES::
  This particular methodology uses the 2nd, 4th, and 6th moments.
  First the 3rd moment normalized by Dm (M4/M3) is calculated.
  Next the shape paramter is computed.
  
  The denominator [2*(Eta-1)] values are filtered between -1E-1 and +1E-1
  to ensure that unrealistically large shape parameter values 
  are not calculated.
    """
#---------------------------------------
    # Calculate the eta ratio [eq. 8]
    Eta = M4**2/(M2 * M6)

    # Now calculate the shape parameter [eq. 3]
    muNumer = (7. - 11. * Eta) - np.sqrt((7. - 11. * Eta)**2 - 4. * (Eta - 1.)
                                    * (30. * Eta - 12.))
    muDenom = 2. * (Eta - 1.)

    # Mask any zero or unrealistically low values in denominator
    muDenom = np.ma.masked_inside(muDenom,-1E-1,1E-1)

    mu = muNumer / muDenom

    return mu
#**====================================================

def slope(M2,M4,mu):
    """Compute the slope (Lambda) using the method of moments

  Ulbrich and Atlas (JAM 1998), Eqn 10
 INPUT::
  M2              = Second moment of the DSD [m^-1]
  M4              = Fourth moment of the DSD [m]
  mu              = Shape parameter of gamma DSD model [unitless]
 OUTPUT::
  Lambda           = Slope [m^-1]
 USAGE::
  Lambda = slope(M2,M4,mu)
 NOTES::
  This particular methodology uses the 3rd and 4th moments and shape
  paramter to compute the slope is computed.
    """
#---------------------------------------
    Lambda = np.ma.sqrt((4 + mu) * (3 + mu) * M2 / M4)

    return Lambda
#**====================================================

def intercept(M6,mu,Lambda):
    """Compute the intercept parameter (N0) using the method of moments 

  Ulbrich and Atlas (JAM 1998), Eqn 6 solved for N0
 INPUT::
  M6              = Sixth moment of the DSD [m^3]
  mu              = Shape parameter of gamma DSD model [unitless]
  Lambda           = Slope [m^-1]
 OUTPUT::
  N0              = Intercept parameter [m^(-1-shape) m^-3]
 USAGE::
  N0 = intercept(M6,mu,Lambda)
 NOTES::
  The exponents work out in the following way:
   M6*Lambda^(7+mu) = m^-3*m^-1(7+mu) = m^-3*m^-7*m^-mu = m^(-1-mu) * m^-3
    """
#---------------------------------------
    # Calculate numerator
    IntNumer = M6 * Lambda**(7. + mu)

    # Calculate denominator
    # Mask values near 0., otherwise gamma function returns "inf"
    mucopy = np.ma.masked_inside(mu,-1E-1,1E-1)
    IntDenom = scifunct.gamma(7 + mucopy)
    # Mask any invalid values
    np.ma.masked_invalid(IntDenom)

    # Mask any zero values from Denom (-999. a problem)
    IntDenom = np.ma.masked_equal(IntDenom,0.)
    #print(mu(:,0)+"  "+Lambda(:,0)+"  "+IntNumer(:,0)+"  "+IntDenom(:,0))

    # Calculate N0
    N0 = IntNumer / IntDenom

    return N0
#**====================================================

def d0(mu,Lambda):
    """Compute the median volume diameter (D0) using the method of moments  

  Ulbrich and Atlas (JAM 1998), Eqn 11
 INPUT::
  mu               = Shape parameter of gamma DSD model [unitless]
  Lambda            = Slope paramter of gamma DSD model [m^-1]
 OUTPUT::
  D0              = Median volume diameter [m]
 USAGE::
  D0 = d0(mu,Lambda)
    """
#---------------------------------------
    D0 = (3.67 + mu) / Lambda

    return D0
#**====================================================

def zr_a(mu,N0):
    """Assuming a rainfall parameter relationship of Z=AR^b,
   Compute the A prefactor using gamma distribution. 
 
  Ulbrich and Atlas (JAMC 2007), Eqn T5
 INPUT::
  mu               = Shape parameter of gamma DSD model [unitless]
  N0               = Intercept parameter [m^(-1-shape) m^-3]
 OUTPUT::
  A                = Z-R prefactor (see description)
 USAGE::
  A = zr_a(mu)
    """
#---------------------------------------
    # Mask 0. values, otherwise gamma function returns "inf"
    mucopy = np.ma.masked_equal(mu,0.)

    # gamma_fix is a patch for versions earlier than NCL v6.2,
    # which cannot handle missing data in the gamma function
    ANumer = 10E6 * scifunct.gamma(7 + mucopy) * N0**(-2.33/(4.67 + mu))
    ADenom = (33.31 * scifunct.gamma(4.67 + mucopy))**((7 + mu)/(4.67 + mu))

    # Mask any zero values from Denom
    ADenom = np.ma.masked_equal(ADenom,0.)

    A = ANumer / ADenom

    return A
#**====================================================

def zr_b(mu):
    """Assuming a rainfall parameter relationship of Z=AR^b,
  Compute the b exponent using gamma distribution.  

  Ulbrich and Atlas (JAMC 2007), Eqn T5
 INPUT::
  mu               = Shape parameter of gamma DSD model [unitless]
 OUTPUT::
  b                = Z-R prefactor (see description)
 USAGE::
  b = zr_b(mu)
    """
#---------------------------------------
    b = (7 + mu) / (4.67 + mu)

    return b
#**====================================================

def norm_intercept(LWC,Dm):
    """Calculates the normalized intercept parameter, which
  is more physically meaningful than N0.  

  Ulbrich and Atlas (JAMC 2007), Eqn T9
 INPUT::
  LWC              = Liquid water content [g m^-3]
  Dm               = Volume weigthed mean diameter [m]
 OUTPUT::
  Nw                = Normalized intercept parameter
 USAGE::
  Nw = norm_intercept(mu)
    """
#---------------------------------------
    # The factor of 1000 converts the water density from kg/m^3 to g/m^3
    Nw = (256 / (np.pi * rhoL * 1000.)) * (LWC / Dm**4)

    return Nw
#====================================================
