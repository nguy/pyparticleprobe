"""
pyparticleprobe.dsd_calcs.axis_ratio
===============================

A grouping of functions to calculate the axis ratio given
#    a drop diameter.  A number of developed relationships is provided.

Adapted by Nick Guy.

"""
#*****************************
#  axis_ratio.py
#*****************************
# HISTORY::
#   12 Feb 2014 - Nick Guy. NRC, NOAA/NSSL (nick.guy@noaa.gov)
#                Converted NCL functions below to Python 
# FUNCTIONS::
# axis_ratio_PB70 - Axis ratio, Pruppacher and Beard (1970)
# axis_ratio_BC87 - Axis ratio, Beard and Chuang (1987)
# axis_ratio_ABL99 - Axis ratio, Andsager et al. (1999)
# axis_ratio_GET94 - Axis ratio, Goddard et al. (1994)
# axis_ratio_BZV02 - Axis ratio, Brandes et al. (2002)
# axis_ratio_TB05 - Axis ratio, Thurai and Bringi (2005)
# axis_ratio_THBRS07 - Axis ratio, Thurai et al. (2007)
#-------------------------------------------------------------------
# Load the needed packages
import numpy as np
#===============================================================
# BEGIN FUNCTIONS
#===============================================================

def axis_ratio_PB70(Diam):
    """Calculates the axis ratio of a droplet (semiminor/semimajor; b/a)
  From Pruppacher and Beard (1970,QJRMS)
 INPUT::
  Diam            = Drop size diameter [in m]
 OUTPUT::
  r               = Axis ratio [mm]
 USAGE::
  r = axis_ratio_PB70(Diam)
 NOTES::
  Linear fit to wind tunnel data, valid (1 <= D <= 9 mm) 

  Falls below bridge experiment data at sizes below 4 mm from 
   Thurai and Bringi (2005, JAOT). Deviates from sphericity attributed to 
   drop oscillations.  Good estimate above.
    """
#---------------------------------------
# Convert the input diameter to mm
    Dmm = Diam * 1000.

# Check the size and calculate axis ratio
    if Dmm >= 1.:# and Dmm <= 9.:
        r = 1.03 - 0.062 * Dmm
    elif Dmm < 1.:
        r = 1.

    return r
#====================================================

def axis_ratio_BC87(Diam):
    """Calculates the axis ratio of a droplet (semiminor/semimajor; b/a)
  From Beard and Chuang (1987, JAS)
 INPUT::
  Diam            = Drop size diameter [in m]
 OUTPUT::
  r               = Axis ratio [mm]
 USAGE::
  r = axis_ratio_BC87(Diam)
 NOTES::
  4th order polynomial fit to numerical model, valid (1 <= D <= 7 mm)

  Upper bound of confidence reasonably close to  bridge experiment data from 
   Thurai and Bringi (2005, JAOT).
    """
#---------------------------------------
# Convert the input diameter to mm
    Dmm = Diam * 1000.

# Check the size and calculate axis ratio
    if Dmm >= 1.:# and Dmm <= 7.:
        r = 1.0048 + 5.7E-4 * Dmm - 2.628E-2 * Dmm**2 + 3.682E-3 * Dmm**3 -1.677E-4 * D**4
    elif Dmm < 1.:
        r = 1.

    return r
#====================================================

def axis_ratio_ABL99(Diam):
    """Calculates the axis ratio of a droplet (semiminor/semimajor; b/a)
  From Andsager et al. (1999, JAS)
 INPUT::
  Diam            = Drop size diameter [in m]
 OUTPUT::
  r               = Axis ratio [mm]
 USAGE::
  r = axis_ratio_ABL99(Diam)
 NOTES::
  2nd order polynomial fit to aircraft measurements, valid (1 <= D <= 4 mm)

  Overestimates in comparison to  bridge experiment data from 
   Thurai and Bringi (2005, JAOT) from D = 2-3.5 mm.
    """
#---------------------------------------
# Convert the input diameter to mm
    Dmm = Diam * 1000.

# Check the size and calculate axis ratio
    if Dmm >= 1. and Dmm <= 4.:
        r = 1.012 - 0.01445 * Dmm - 1.028E-2 * Dmm**2
    elif Dmm < 1.:
        r = 1.
    elif Dmm > 4.:
        r = np.nan

    return r
#====================================================

def axis_ratio_GET94(Diam):
    """Calculates the axis ratio of a droplet (semiminor/semimajor; b/a)
  From Goddard et al (1994, IEEE)
 INPUT::
  Diam            = Drop size diameter [in m]
 OUTPUT::
  r               = Axis ratio [mm]
 USAGE::
  r = axis_ratio_GET94(Diam)
 NOTES::
  3rd order polynomial fit inferred from dual-polarimetric radar and
   Joss disdrometer data, valid (1 <= D <= 5 mm)

  Closest to bridge experiment data from Thurai and Bringi (2005, JAOT)
    """
#---------------------------------------
# Convert the input diameter to mm
    Dmm = Diam * 1000.

# Check the size and calculate axis ratio
    if Dmm >= 1.:# and Dmm <= 5.:
        r = 1.075 - 6.5E-2 * Dmm - 3.6E-3 * Dmm**2 + 4.0E-3 * Dmm**3
    elif Dmm < 1.:
        r = 1.

    return r
#====================================================

def axis_ratio_BZV02(Diam):
    """Calculates the axis ratio of a droplet (semiminor/semimajor; b/a)
  From Brandes et al (2002, JAM)
 INPUT::
  Diam            = Drop size diameter [in m]
 OUTPUT::
  r               = Axis ratio [mm]
 USAGE::
  r = axis_ratio_BZV02(Diam)
 NOTES::
  4th order polynomial fit based upon many previous experimental data,
    valid (1 <= D <= 7 mm)

  Close to bridge experiment data from Thurai and Bringi (2005, JAOT) over 
   valid range mentioned above.
    """
#---------------------------------------
# Convert the input diameter to mm
    Dmm = Diam * 1000.

# Check the size and calculate axis ratio
    if Dmm >= 1.:# and Dmm <= 7.:
        r = 0.9951 + 2.51E-2 * Dmm - 3.644E-2 * Dmm**2 + 5.303E-3 * Dmm**3 - 2.492E-4 * Dmm**4
    elif Dmm < 1.:
        r = 1.

    return r
#====================================================

def axis_ratio_TB05(Diam):
    """Calculates the axis ratio of a droplet (semiminor/semimajor; b/a)
  From Thurai and Bringi (2005, JAOT)
 INPUT::
  Diam            = Drop size diameter [in m]
 OUTPUT::
  r               = Axis ratio [mm]
 USAGE::
  r = axis_ratio_TB05(Diam)
 NOTES::
  4th order polynomial fit based upon 80-m bridge experiment, 
   valid (1.5 <= D <= 8 mm)
 
  This offers similar performance as the Brandes et al 2002 formulation:
   axis_ratio_BZV05
    """
#---------------------------------------
# Convert the input diameter to mm
    Dmm = Diam * 1000.

# Check the size and calculate axis ratio
    if Dmm >= 1.5:# and Dmm <= 7.:
        r = 0.9707 + 4.26E-2 * Dmm - 4.29E-2 * Dmm**2 + 6.5E-3 * Dmm**3 - 3.0E-4 * Dmm**4
    elif Dmm < 1.5:
        r = 1.

    return r
#====================================================

def axis_ratio_THBRS07(Diam):
    """Calculates the axis ratio of a droplet (semiminor/semimajor; b/a)
  From Thurai et al. (2007, JAOT)
 INPUT::
  Diam            = Drop size diameter [in m]
 OUTPUT::
  r               = Axis ratio [mm]
 USAGE::
  r = axis_ratio_BZV02(Diam)
 NOTES::
  4th order polynomial fits based upon 80-m bridge experiment, 
   valid (1.5 <= D <= 9 mm) and laboratory measurements of Beard and
   Kubesh (1991) valid (0.7 <= D <=1.5 mm)
  Below 0.7 mm is assumed spherical
 
  This offers very good performance.
    """
#---------------------------------------
# Convert the input diameter to mm
    Dmm = Diam * 1000.

# Check the size and calculate axis ratio
    if Dmm >= 1.5:# and Dmm <= 9.: # Eqn 2
        r = 1.065 - 6.25E-2 * Dmm - 3.99E-3 * Dmm**2 + 7.66E-4 * Dmm**3 - 4.095E-5 * Dmm**4
    elif Dmm < 1.5 and Dmm >= 0.7: # Eqn 3
        r = 1.173 - 0.5165 * Dmm + 0.4698 * Dmm**2 - 0.1317 * Dmm**3 - 8.5E-3 * Dmm**4
    elif Dmm < 0.7:
        r = 1.

    return r
#====================================================

