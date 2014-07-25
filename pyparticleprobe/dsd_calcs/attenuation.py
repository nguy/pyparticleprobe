"""
pyparticleprobe.dsd_calcs.attenuation
=========================

A grouping of functions that calcuates attenuation characteristics

7 Feb 2014 - Adapted by Nick Guy NOAA/NSSL/WRDD, NRC

"""
# HISTORY::
#   3 Feb 2014 - Nick Guy. NOAA/NSSL, NRC (nick.guy@noaa.gov) 
#
# NOTES::
#   Arrays seem to be able to be passed, but make sure they are float arrays
#    (e.g. created with numpy) and not lists
#
# FUNCTIONS::
# abs_coeff - Absorption coefficient
# scat_coeff - Scattering coefficient
# ext_coeff - Extinction coefficient
# spec_atten - Specific attenuation
#-------------------------------------------------------------------
# Load the needed packages
import numpy as np
#===============================================================
# BEGIN FUNCTIONS
#===============================================================

def abs_coeff(D,lam,m):
    """Absorption coefficient of a spherical particle
  From Doviak and Zrnic (1993), Eqn 3.14a or Battan (1973), Eqn 6.6
 INPUT::
  D             = Particle diameter [m]
  lam           = Radar wavelength [m]
  m             = Complex refractive index [unitless]
 OUTPUT::
  Qa            = Absorption coefficient
 USAGE::
  Qa = abs_coeff(D,lam,m)
 NOTES::
  The default is for a dielectric factor value for water.  This can be 
   changed by the user, e.g. K=0.208 for particle sizes of equivalent melted
   diameters or K=0.176 for particle sizes of equivalent ice spheres.
    """
#---------------------------------------
    Km = (m**2 - 1) / (m**2 + 2)
    Qa = (np.pi**2 * D**3 / lam) * np.imag(-1 * Km)

    return Qa
#====================================================

def scat_coeff(D,lam,m):
    """Scattering coefficient of a spherical particle
  From Doviak and Zrnic (1993), Eqn 3.14b or Battan (1973), Eqn 6.5
 INPUT::
  D             = Particle diameter [m]
  lam           = Radar wavelength [m]
  m             = Complex refractive index [unitless]
 OUTPUT::
  Qs            = Scattering coefficient
 USAGE::
  Qs = scat_coeff(D,lam,m)
    """
#---------------------------------------
    Km = (m**2 - 1) / (m**2 + 2)
    Qs = (2 * np.pi**5 * D**6 / (3 * lam**4) * (np.absolute(Km))**2)

    return Qs
#====================================================

def ext_coeff(D,lam,m):
    """Extinction coefficient of a spherical particle
  From Doviak and Zrnic (1993), Eqn 3.14b or Battan (1973), Eqn 6.5
 INPUT::
  D             = Particle diameter [m]
  lam           = Radar wavelength [m]
  m             = Complex refractive index [unitless]
 OUTPUT::
  Qe            = Scattering coefficient
 USAGE::
  Qe = ext_coeff(D,lam,m)
 NOTES::
  The default is for a dielectric factor value for water.  This can be 
   changed by the user, e.g. K=0.208 for particle sizes of equivalent melted
   diameters or K=0.176 for particle sizes of equivalent ice spheres.
    """
#---------------------------------------
    Qa = abs_coeff(D,lam,m)
    Qs = scat+coeff(D,lam,m)
    Qe = Qa + Qs

    return Qe
#====================================================

def spec_atten(Nd,Diam,lam,m):
    """Extinction coefficient of a spherical particle
  From Doviak and Zrnic (1993), Eqn 3.15
 INPUT::
  Nd            = Drop concentration as a function of drop size [m^-3]
  Diam          = Drop size diameter [mm]
  lam           = Radar wavelength [m]
  m             = Complex refractive index [unitless]
 OUTPUT::
  K            = Specific attenuation [dB/km]
 USAGE::
  K = spec_atten(Nd,Diam,lam,m)
 NOTES::
  The default is for a dielectric factor value for water.  This can be 
   changed by the user, e.g. K=0.208 for particle sizes of equivalent melted
   diameters or K=0.176 for particle sizes of equivalent ice spheres.
    """
#---------------------------------------
    Qa = abs_coeff(Diam,lam,m)
    Qs = scat_coeff(Diam,lam,m)
    Qe = Qa + Qs
    
    # Calculate specific attenuation
    K = 4.34e3 * Nd * Qe

    return Qe
#====================================================
