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
#                Converted NCL functions below to Python originally coded
#                 between Jan - Sep 2013.
#   9 Mar 2014 - NG: Added optional summation along axis
#  28 Mar 2014 - NG: Updated code to make more "pythonic"
#   9 Jul 2014 - NG: Refactored code to take advantage of the power of 
#                    object-oriented code.
#                    Added a couple of functions to pull out different groups
#                     of variables for use in pyparticleprobe package.
# NOTES::
#
# FUNCTIONS::
# term_vel - Hydrometeor terminal velocity
# refl - Reflectivity factor
# rainrate - Rainfall rate
# watermass - Liquid water content (mass)
# norm_intercept - Normalized intercept parameter (input moments)
# norm_intercept2 - Normalized intercept parameter (input Conc, Diam)
# norm_intercept3 - Normalized intercept parameter (Alternate calc)
# get_params - Calculate a number of parameters
#-------------------------------------------------------------------
# Load the needed packages
import numpy as np
from moments import *
import scipy.special as scifunct
#-------------------------------------------------------------------
# Define various constants that may be used for calculations
rho0=1.2 # reference density of air at the surface [kg m^-3]
rhoL=1000. # density of water [kg m^-3]
#**MAY WANT TO USE AtmConst or other???****
#
#===============================================================
# BEGIN FUNCTIONS
#===============================================================

def term_vel(Diam,Type,Opt,Correct,Rhoair,Factor):
    """Hydrometeor terminal velocity is calculated given drop size distribution
    INPUT::
     Diam            = Drop size diameter [in mm]
     Type            = Type or precipitation 0 for rain, 1 for ice (Integer)
     Opt             = Method to use for calculation (Integer)
     Correct         = Apply altitude correction? 0=No, 1=Yes
     Rhoair          = Air density [kg/m^3], same dimension as Diam
                    If Correct=0, then set to any number.
     Factor          = Density correction factor (e.g. 0.40)
                    If Correct=0, then set to any number.
    OUTPUT::
     Vt              = Terminal velocity calculations (same dim as Diam) [m/s]
    USAGE::
     Vt = term_vel(Diam,Type,Opt,Correct,Rhoair,Factor)
    NOTES::
     If rain (Type=0) is chosen there are calculation options:
        Opt = 1  ==> Based upon empirical formulae from Brandes et al 
                     (JAM# 2002), which found a polynomial fit to the
                     Gunn and Kinzer (J. Meteor.# 1949) and Pruppacher and
                     Pitter (JAS# 1971) laboratory measurements.
            = 2  ==> Exponential fit derived by Lhermitte (GRL# 1988), 
                     accurate to within 3 cm/s between 0.5 - 6.0 mm.
            = 3  ==> Atlas and Ulbrich (# 1977)
            = 4  ==> Atlas et al. (Rev. Geophys.# 1973)
     If ice (Type=1) is chosen these are the calculation options:
        Opt = 1  ==> Khain et al. (Atmos. Res.# 2000),Kumjian et al. (JAS# 2012)
    """
#---------------------------------------
# Calculate terminal velocity values
    if Type == 0: # Rain
        if Opt == 1:
            Vt = -0.1021 + 4.932 * Diam - 0.9551 * Diam**2. + 0.07934 * Diam**3.
            - 0.002362 * Diam**4.
        elif Opt == 2:
            Vt = 9.25 * (1 - np.exp(-6.8 * (Diam/10.)**2 - 4.88 * (Diam/10.)))
        elif Opt == 3:
            Vt = 3.78 * Diam**0.67
        elif Opt == 4:
            Vt = 9.65 - 10.3 * np.exp(-6 * Diam)
    elif Type == 1: # Ice
        if Opt == 1:
            Vt = 0.2259 + 1.5954 * Diam - 0.0405 * Diam**2.

# Apply altitude correction if desired
# First introduced in Foote and duToit (JAM# 1969) and expanded in Beard (JTech# 1985).  This correction applies a density correction factor for finding terminal velocity of raindrops aloft and not just at sea level.  The exponential factor given in Beard (1985) ranges from 0.41 - 0.46, with suggested values of 0.42 for the calculation of rainfall rate and 0.45 for mean

# Calculate Corrected velocity?
    if Correct == 0:
        Vtf = Vt
    elif Correct == 1:
        Vtf = Vt * (rho0 / Rhoair)**Factor

    return Vtf
#**====================================================

def refl(Conc,Diam,sum=False,axis=None):
    """Reflectivity factor computed given drop size distribution 
    INPUT::
     Conc            = Drop concentraton as a function of drop size [m^-3]
     Diam            = Drop size diameter [in mm]
       OPTIONAL
     sum             = True results in a summed variable additionally returned
     axis            = The axis to sum over is sum set to True
    OUTPUT::
     Ze              = Reflectivity [mm^6 * m^-3]
     Ze_sum          = Optional returned array, Ze summed along an axis,masked = 0
    USAGE::
     Ze = refl(Conc,Diam)
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
    """
#---------------------------------------
# Calculate equivalent reflectivity
    Ze = Conc * Diam**6.

    if sum:
        Ze_sum = np.sum(Ze,axis=axis)
        Ze_sum = np.ma.masked_equal(Ze_sum,0)
        
        return Ze,Ze_sum
    else:
        return Ze
#**====================================================

def rainrate(Conc,Diam,Vt,sum=False,axis=None):
    """Rainfall rate computed given drop size distribution
    INPUT::
     Conc            = Drop concentraton as a function of drop size [m^-3]
     Diam            = Drop size diameter [mm]
     Vt              = Terminal velocity [m/s]
       OPTIONAL
     sum             = True results in a summed variable additionally returned
     axis            = The axis to sum over is sum set to True
    OUTPUT::
     RR              = Rainfall rate [mm/hr]
     RR_sum          = Optional returned array, RR summed along an axis,masked = 0
    USAGE::
     RR = rainrate(Conc,Diam,Vt)
    NOTES::
     Input variables must have same dimensionality

     Third moment calculation of the drop size distribution  
     Calculation comes from integrating (summing) over the 
     RR = (pi/6) *Summation(Integral)[i=0,n] N(Di) * Di^3 * Vt(Di)
     This form yields m/s.  Because DSD is in mm and it is common to 
     to have the form mm/hr, a conversion factor of 3600E-6 is included.
       n
      ---
      \              
      /   (pi/6) *N(Di) * Di^3 * Vt(Di)
      ---
      i=0
    """
#---------------------------------------
# Calculate rainfall rate
# 3600E-6 is conversion factor given D in mm and Vt m/s to get to mm/hr
    RR = (np.pi/6.)  *3600E-6 * (Conc * Diam**3.) * Vt

    if sum:
        RR_sum = np.sum(RR,axis=axis)
        RR_sum = np.ma.masked_equal(RR_sum,0)
        
        return RR,RR_sum
    else:
        return RR
#**====================================================

def watermass(Conc,Diam,sum=False,axis=None):
    """Liquid water content (mass) computed given drop size distribution
    INPUT::
     Conc            = Drop concentraton as a function of drop size [m^-3]
     Diam            = Drop size diameter [in m]
       OPTIONAL
     sum             = True results in a summed variable additionally returned
     axis            = The axis to sum over is sum set to True
    OUTPUT::
     W               = Liquid water mass [g/m^3] at each bin!
     LWC             = Optional returned array, W summed along an axis, masked = 0
    USAGE::
     W = watermass(Conc,Diam)
    NOTES::
     Input variables must have same dimensionality
  
     Third moment calculation of the drop size distribution
     Calculation comes from integrating (summing) over the 
     W = (pi * rhoL/6) *Summation(Integral)[i=0,n] N(Di) * Di^3
  
       n
      ---
      \              
      /   (pi *rhoL/6) * N(Di) * Di^3
      ---
      i=0

     LWC can be retrieved by summing along the bin (diameter) axis
    """
#---------------------------------------
# Calculate water mass
# Multiply rhoL by 1000 to get in g/m^3
    W = (np.pi / 6.) * (rhoL * 1000) * (Conc * Diam**3.)

    if sum:
        LWC = np.sum(W,axis=axis)
        LWC = np.ma.masked_equal(LWC,0)
        
        return W,LWC
    else:
        return W,LWC
#**====================================================

def norm_intercept(M3,M4):
    """Calculates the normalized intercept parameter, which
    is more physically meaningful than N0.  
    as in Testud et al. (JAM# 2001), Eqn 11 
    INPUT::
     M3              = Third moment of DSD [m]
     M4              = Fourth moment of DSD [m]
    OUTPUT::
     Nw              = Normalized intercept parameter [m^-1 m^-3]
    USAGE::
     Nw = norm_intercept(M3,M4)
    """
# HISTORY::
#  9 Sep 2013 - Nick Guy NOAA/NSSL/WRDD, NRC
#---------------------------------------
  # Replace any zero number with missing data in denominator
    M4Denom = M4**4
    M4Denom = np.ma.masked_equal(M4Denom,0.)

  # Calculate Nw [eq. 15]
    Nw = (4.**4/scifunct.gamma(4)) * (M3**5/M4Denom)

    return Nw
#**====================================================

def norm_intercept2(Conc,Diam):
    """Calculates the normalized intercept parameter, which
    is more physically meaningful than N0.  
    as in Testud et al. (JAM# 2001) , Eqn 11
    INPUT::
     Conc            = Drop concentraton as a def of drop size [m^-3]
     Diam            = Drop size diameter [in m]
    OUTPUT::
     Nw              = Normalized intercept parameter [m^-1 m^-3]
    USAGE::
     Nw = norm_intercept2(Conc,Diam)
    NOTES::
     This routine calculates at each diameter if the input is such, in other
      words the values are not summed over the diameter range as should be 
      done for many calculations.
    """
#---------------------------------------
# Calculate the 3rd and 4th order moments
    M3 = moment_nth(Conc,Diam,3,'m')
    M4 = moment_nth(Conc,Diam,4,'m')

# Replace any zero number with missing data in denominator
    M4Denom = M4**4
    M4Denom = np.ma.masked_values(M4Denom,0.)

# Calculate Nw [eq. 15]
    Nw = (4.**4/(scifunct.gamma(4))) * (M3**5/M4Denom)

    return Nw
#**====================================================

def norm_intercept3(Conc,Diam):
    """Calculates the normalized intercept parameter, which
    is more physically meaningful than N0.  
    as in Testud et al. (JAM 2001), Eqn 8
    used by Bringi et al. (JAS 2003), Eqn 4
    INPUT::
      Conc            = Drop concentraton as a def of drop size [m^-3]
      Diam            = Drop size diameter [in m]
      OUTPUT::
      Nw              = Normalized intercept parameter [m^-1 m^-3]
    USAGE::
      Nw = norm_intercept3(Conc,Diam)
    NOTES::
     This routine calculates at each diameter if the input is such, in other
      words the values are not summed over the diameter range as should be 
      done for many calculations.
    """
#---------------------------------------
# Calculate the water mass
    W = watermass(Conc,Diam)
    W = np.ma.masked_values(W, 0.)
# Calculate the 4th order moments
    M3 = moment_nth(Conc,Diam,3,'m')
    M4 = moment_nth(Conc,Diam,4,'m')
    
# Calculate Dm
    Dm = M4/M3

# Replace any zero number with missing data in denominator
    M4Denom = Dm**4
    M4Denom = np.ma.masked_values(M4Denom,0.)

# Calculate Nw (the factor of 1000 is to convert water density units
    Nw = (4.**4/(np.pi * rhoL * 1000.)) * (W/M4Denom)

    return Nw
#**====================================================

def get_drops_params(Conc,Diam,RhoAir):
    """Calculate a series of parameters using DSD data.  This function can be used to drive
    many other functions in this set of routines.  
    INPUT::
     Conc            = Drop concentraton as a def of drop size [m^-3]
     Diam            = Drop size diameter [in m]
     RhoAir          = Air density [kg/m^3], same dimension as Diam
    OUTPUT::
     Data            = Dictionary containing the following values:
  
        Ze              = Reflectivity [mm^6 * m^-3]
        dBZe            = Reflectivity [dBZ]
        RR              = Rainfall rate [mm/hr]
        MWater          = Mass of droplet [g]
    USAGE::
     Data = get_params_drops(Conc,Diam,RhoAir,W)
    NOTES::
     The arrays returned are masked arrays.
     The Conc,Diam,RhoAir should all be the same shape
    """
#--------------------------------------------------
    # Using moment calculations, estimate Z and R via distribution (N,D)
    #--------------------------------------------
    # Compute equivalent reflectivity and masked summed time series (masked Z = 0)
    Ze = refl(Conc,Diam*1000.,sum=False)

    # Compute terminal vel for R (run 2x for water and ice-note not robust ice method)
    Vt = term_vel(Diam*1000.,0,1,1,RhoAir,0.42)

    # Compute rainrate and masked summed time series (masked RR = 0)
    RR = rainrate(Nd,Diam*1000.,Vt,sum=False)

    #---Now estimate gamma parameters---
    # Compute moments of DSD aConc a masked summed time series (masked M = 0)
    M2 = moment_nth(Conc,Diam,2,"m",sum=False)
    M3 = moment_nth(Conc,Diam,3,"m",sum=False)
    M4 = moment_nth(Conc,Diam,4,"m",sum=False)
    M6 = moment_nth(Conc,Diam,6,"m",sum=False)

    # Compute the water mass and masked summed time series (LWC) (masked water = 0)
    Mwater = watermass(Conc,Diam,sum=False)

    # Delete some variables to clean up
    del M2,M3,M4,M6

    # Calculate a reflectivity in dBZ
    dBZe = 10. * np.log10(Ze)
        
    # Create a dictionary to transfer the data
    data = {'Ze' : Ze,
            'dBZe' : dBZe,
            'rain_rate' : RR,
            'water_mass' : MWater}
    
    return data
#====================================================

def get_bins_params(Conc,sumaxis=0):
    """Calculate a series of parameters using DSD data.  This function can be used to drive
    many other functions in this set of routines.  
    INPUT::
     Conc            = Drop concentraton as a def of drop size [m^-3]
     sumaxis         = Index of axis to sum over for 
    OUTPUT::
     Data            = Dictionary containing the following values:
  
        Nd              = Total concentration [m^-3] summed at each bin (masked = 0)
        Nd_Avg          = Average concentration [m^-3] summed at each bin (masked = 0)
        Nd_Med          = Average concentration [m^-3] summed at each bin (masked = 0)
    USAGE::
     data = get_param_bins(Conc,sumaxis=1)
    NOTES::
     The arrays returned are masked arrays.
    """
#--------------------------------------------------
# Sum over entire data set time period for total drop concentration
    Nd = np.sum(Conc,axis=sumaxis)

    # Average/Median drop concentration valuesfor the data set time period
    Nd_Avg = np.mean(Conc,axis=sumaxis)
    Nd_Med = np.median(Conc,axis=sumaxis)

    # Mask the zero values
    Nd = np.ma.masked_equal(Nd,0)
    Nd_Avg = np.ma.masked_equal(Nd_Avg,0)
    Nd_Med = np.ma.masked_equal(Nd_Med,0)
        
    # Create a dictionary to transfer the data
    data = {'Nd' : Nd,
            'Nd_Avg' : Nd_Avg,
            'Nd_Med' : Nd_Med}
    
    return data
#====================================================

def get_timeseries_params(Conc,Diam,RhoAir,sumaxis=0):
    """Calculate a series of parameters using DSD data.  This function can be used to drive
    many other functions in this set of routines.  
    INPUT::
     Conc            = Drop concentraton as a def of drop size [m^-3]
     Diam            = Drop size diameter [in m]
     RhoAir          = Air density [kg/m^3], same dimension as Diam
     sumaxis         = Index of axis to sum over for 
    OUTPUT::
     Data            = Dictionary containing the following values:
  
         Ze_sum          = Reflectivity [mm^6 * m^-3] summed along an axis (masked = 0)
         dBZ             = Reflectivity [dBZ] summed along an axis
         RR_sum          = Rainfall rate [mm/hr] summed along an axis (masked = 0)
         LWC             = Liquid water content [g/m^3]
         Nw              = Normalized gamma intercept [mm^-1 m^-3]
         Dm              = Mean mass-weighted diameter [m]
         Nt              = Total concentration per time [m^-3]
    USAGE::
     data = get_param_timeseries(Conc,Diam,RhoAir,sumaxis=1)
    NOTES::
     The arrays returned are masked arrays.
     The Conc,Diam,RhoAir should all be the same shape
    """
#--------------------------------------------------
    # Using moment calculations, estimate Z and R via distribution (N,D)
    #--------------------------------------------
    # Compute equivalent reflectivity and masked summed time series (masked Z = 0)
    Ze, Ze_sum = refl(Conc,Diam*1000.,sum=True,axis=sumaxis)

    # Compute terminal vel for R (run 2x for water and ice-note not robust ice method)
    Vt = term_vel(Diam*1000.,0,1,1,RhoAir,0.42)

    # Compute rainrate and masked summed time series (masked RR = 0)
    RR, RR_sum = rainrate(Conc,Diam*1000.,Vt,sum=True,axis=sumaxis)

    #---Now estimate gamma parameters---
    # Compute moments of DSD and a masked summed time series (masked M = 0)
    M2, M2_sum = moment_nth(Conc,Diam,2,"m",sum=True,axis=sumaxis)
    M3, M3_sum = moment_nth(Conc,Diam,3,"m",sum=True,axis=sumaxis)
    M4, M4_sum = moment_nth(Conc,Diam,4,"m",sum=True,axis=sumaxis)
    M6, M6_sum = moment_nth(Conc,Diam,6,"m",sum=True,axis=sumaxis)

    # Calculate the volume (or mass) -weighted mean diameter Dm = M4/M3
    Dm = M4_sum/M3_sum

    # Compute the water mass and masked summed time series (LWC) (masked water = 0)
    Mwater, LWC = watermass(Conc,Diam,sum=True,axis=sumaxis)

    # Compute the normalized gamma intercept, Nw
    Nw = norm_intercept(M3_sum,M4_sum)

    # Delete some variables to clean up
    del M2,M3,M4,M6,M2_sum,M3_sum,M4_sum,M6_sum

    # Calculate a reflectivity in dBZ
    dBZ_sum = 10. * np.ma.log10(Ze_sum)
    
    # Sum over all bins for a concentration per time step
    Nt_sum = np.sum(Conc,axis=sumaxis)
        
    # Create a dictionary to transfer the data
    data = {'Ze' : Ze_sum,
            'dBZ' : dBZ_sum,
            'RR' : RR_sum,
            'Dm' : Dm,
            'LWC' : LWC,
            'Nw' : Nw,
            'Nt' : Nt_sum}
    
    return data
#====================================================

def get_model_params(Conc,Diam,RhoAir,modelName=None,sumaxis=0):
    """Calculate a series of parameters using DSD data.  This function can be used to drive
    many other functions in this set of routines.  
    INPUT::
     Conc            = Drop concentraton as a def of drop size [m^-3]
     Diam            = Drop size diameter [in m]
     RhoAir          = Air density [kg/m^3], same dimension as Diam
     W               = Vertical Velocity
     modelName       = String of fit model to use for DSD.  Choices are:
                        'gamma_ua98' - Ulbrich and Atlas 1998 gamma dist method
                        'gamma_ts96' - Tokay and Short 1996 gamma dist method
     sumaxis         = Index of axis to sum over for 
    OUTPUT::
     Data            = Dictionary containing the following values:
  
         D0              = Median volume diameter [mm]
         shape           = Model fit shape parameter
         slope           = Model fit slope parameter
         N0              = Model fit intercept parameter
    USAGE::
     data = get_model_params(Conc,Diam,RhoAir,model='gamma_ua98'sumaxis=1)
    NOTES::
     The arrays returned are masked arrays.
     The Conc,Diam,RhoAir should all be the same shape
    """
#--------------------------------------------------
    #---Now estimate parameters---
    # Compute moments of DSD and a masked summed time series (masked M = 0)
    M2, M2_sum = moment_nth(Conc,Diam,2,"m",sum=True,axis=sumaxis)
    M3, M3_sum = moment_nth(Conc,Diam,3,"m",sum=True,axis=sumaxis)
    M4, M4_sum = moment_nth(Conc,Diam,4,"m",sum=True,axis=sumaxis)
    M6, M6_sum = moment_nth(Conc,Diam,6,"m",sum=True,axis=sumaxis)
    
    # Check which model to use to fit data for calculations
    if modelName == 'gamma_ua98':
        import ua98 as model
        shape = model.shape(M2_sum,M4_sum,M6_sum) # Shape parameter
        slope = model.slope(M2_sum,M4_sum,shape)  # Slope parameter
        D0 = model.d0(shape,slope)                # Median volume diameter
        N0 = model.intercept(M6_sum,shape,slope)  # Intercept parameter
        
    elif modelName == 'gamma_ts96':
        import ts96 as model
        shape = model.shape(M3_sum,M4_sum,M6_sum) # Shape parameter
        slope = model.slope(M3_sum,M4_sum,shape)  # Slope parameter
        D0 = model.d0(shape,slope)                # Median volume diameter
        N0 = model.intercept(M3_sum,shape,slope)  # Intercept parameter
        
    elif modelName == 'exp':
        import exp_dsd as model
        slope = model.slope(M3_sum,M6_sum)
        N0 = model.intercept(M3_sum,M6_sum)
        shape = []
        D0 = []
        
    else:
        print "Need to choose a model, choices are: "
        print "'gamma_ua98'"
        print "'gamma_ts96'"
        print "'exp'"
        return

    # Get rid of "inf" and "nan" values from computations
    N0 = np.ma.masked_invalid(N0)
    slope = np.ma.masked_invalid(slope)
        
    # Delete some variables to clean up
    del M2,M3,M4,M6,M2_sum,M3_sum,M4_sum,M6_sum

    # Create a dictionary to transfer the data
    data = {'shape' : shape,
            'slope' : slope,
            'N0' : N0,
            'D0' : D0}
    
    return data
#====================================================