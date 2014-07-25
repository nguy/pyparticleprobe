"""
pyparticleprobe.dsd_calcs.partition
=========================

A grouping of functions that paritions time series of data into Convective-Stratiform-Transitional rain components

Adapted by Nick Guy.

"""
# HISTORY::
#   19 Feb 2014 - Nick Guy. NOAA/NSSL, NRC (nick.guy@noaa.gov)
# FUNCTIONS::
# Rmag_cs - Testud et al 2001 Con-Sf, Rain rate magnitude
# Rsd_cs - Bringi et al 2003 Con-Sf, Rain rate Std Dev.
# RmagSD_cs - Islam et al 2012 Con-Sf, R magnitude and Std Dev.
# Wvel_cs - Atlas et al 2000 Con-Sf, Vert Velocity magnitude
# Nw_D0_cs - Bringi et al 2009, Con-Sf, Nw-D0 space separation
# Nw_D0_cst - Bringi et al2009, Thurai et al 2010, Con-Sf-Trans,Nw-D0 separation
#-------------------------------------------------------------------
# Load the needed packages
import numpy as np
#===============================================================
# BEGIN FUNCTIONS
#===============================================================

def Rmag_cs(R,Thresh=10.,window=3):
    """Methodology used in Testud et al. (JAM 2001). 
  At point R(k), If all points in range R(k-5) -  R(k+5) < 10 mm/h, 
   then R(k) is considered stratiform, otherwise R(k) is convective.  This
   study used aircraft data, so the +/- 5 corresponds to +/- 3.6 km 
   denoted as the "radius of influence" of a convective cell.
 INPUT::
  R             = Rain rate time series [mm/h]
  Thresh        = Rain rate threshold for which above is Convective ID  [mm/h]
  window        = Size of window to analyze in time series
 OUTPUT::
  TS            = Array identifying Conv-Strat Classification
                  1 = Stratiform
                  2 = Convective
 USAGE::
  TS = T01_cs(R,[Thresh=10.],[window=3])
 NOTES::
  Default values:
   The threshold is that used in Testud et al. 2001.
   The window size has been calibrated to match a radius of 2.52 km 
    travelled by the NOAA P-3.  These can be changed when calling the function.
    """
#---------------------------------------
    # Create the variable to pass from function
    TS = np.zeros_like(R)
    # Find the number of points in the time series
    nPts = len(R)

    # Create a padded array and apply the mirroring technique for end points
    #  Caveat should be noted that if variability exists at end points 
    #  (convective cell) it may be  "watered down" or amplified.
    R_pad = np.pad(R,window,'reflect')

    for x in range(nPts):
        # Set the beginning and end indices according to added padding
        Inds = np.arange(2 * window + 1) + x

        # Create an array of threshold values the same size as window
        #  and note the addition of +1 to account for current time step
        ThArr = np.zeros(2 * window + 1) + Thresh

        if np.less([R_pad[Inds]],[ThArr]).all():
            TS[x] = 1  # Stratiform
        else:
            TS[x] = 2  # Convective
        del Inds, ThArr # Clean up the indices just in case

    return TS
#====================================================

def Rsd_cs(R,Thresh=1.5,window=4):
    """Methodology used in Bringi et al. (JAS 2003). 
  The standard deviation of rain rate over 5 consecutive samples is computed.
   Std Dev <= 1.5 mm/h is considered stratiform
   Std Dev >= 1.5 mm/h is considered convective
 INPUT::
  R             = Rain rate time series [mm/h]
  Thresh        = Standard deviation of Rain rate threshold 
                   for which above is Convective ID  [mm/h]
  window        = Size of window to analyze in time series
 OUTPUT::
  TS            = Array identifying Conv-Strat Classification
                  1 = Stratiform
                  2 = Convective
 USAGE::
  TS = B03_cs(R,[Thresh=1.5],[window=4])
    """
    
#---------------------------------------
    # Create the variable to pass from function
    TS = np.zeros_like(R)
    # Find the number of points in the time series
    nPts = len(R)

    # Create a padded array and apply the mirroring technique for end points
    #  Caveat should be noted that if variability exists at end points 
    #  (convective cell) it may be  "watered down" or amplified.
    R_pad = np.pad(R,(0,window),'reflect')

    for x in range(nPts):
        # Set the beginning and end indices according to added padding
        Inds = np.arange(window) + x
        
        if R_pad[Inds].std() <= Thresh:
            TS[x] = 1  # Stratiform
        else:
            TS[x] = 2  # Convective
        del Inds # Clean up the indices just in case

    return TS
#====================================================

def RmagSD_cs(R,ThreshR=10.,ThreshSD=1.5,window=4):
    """Methodology used in Islam et al. (Atmos Res 2012).  The methodologies of
   Testud et al. (2001) and Bringi et al. (2003) are combined.
 INPUT::
  R             = Rain rate time series [mm/h]
  ThreshR       = Rain rate threshold for which above is Convective ID  [mm/h]
  ThreshSD      = Standard deviation of Rain rate threshold 
                   for which above is Convective ID  [mm/h]
  window        = Size of window to analyze in time series
 OUTPUT::
  TS            = Array identifying Conv-Strat Classification
                  1 = Stratiform
                  2 = Convective
 USAGE::
  TS = I12_cs(R,[ThreshR=10.],[ThreshSD=1.5],[window=4])
    """
    
#---------------------------------------
    # Create the variable to pass from function
    TS = np.zeros_like(R)
    # Find the number of points in the time series
    nPts = len(R)

    # Create a padded array and apply the mirroring technique for end points
    #  Caveat should be noted that if variability exists at end points 
    #  (convective cell) it may be  "watered down" or amplified.
    R_pad = np.pad(R,window,'reflect')

    for x in range(nPts):
        # Set the beginning and end indices according to added padding
        Inds = np.arange(2 * window + 1) + x

        # Create an array of threshold values the same size as window
        #  and note the addition of +1 to account for current time step
        ThArr = np.zeros(2 * window + 1) + ThreshR

        if np.less([R_pad[Inds]],[ThArr]).all() and R_pad[Inds].std() <= ThreshSD:
            
            TS[x] = 1  # Stratiform
        else:
            TS[x] = 2  # Convective
        del Inds, ThArr # Clean up the indices just in case

    return TS

#====================================================

def Wvel_cs(VertVel,Wthresh=1.0):
    """Methodology used in Atlas et al. (JGR 2000).
  This study used aircraft data.  Separation was based upon vertical
   velocity measurements obtained by aircraft instrumentation. 
   In theory, these thresholds could be adjusted to use alternate data,
   though care should be taken.
 INPUT::
  VertVel       = Vertical velocity time series [m/s]
  Wthresh       = Vertical velocity  threshold for which 
                   above is Convective ID  [m/s]
 OUTPUT::
  TS            = Array identifying Conv-Strat Classification
                  1 = Stratiform
                  2 = Convective
 USAGE::
  TS = A00_cs(VertVel,[Thresh=1.3])
 NOTES::
  Ensure that time series and threshold values have same units!
  Default values:
   The A00 paper uses a Wthresh = 1.0 m/s.  This was found by looking at
    6-s average vertical velocity measurements as a function of A coefficients
    of a Z-R relationship.
   Wthresh = 1.3 m/s used for DYNAMO project comparison.
    """
#---------------------------------------
    # Create the variable to pass from function
    TS = np.zeros_like(VertVel)
    # Find the number of points in the time series
    nPts = len(VertVel)

    TS = np.where(abs(VertVel) > Wthresh, 2, 1)

    return TS
#====================================================

def Nw_D0_cs(Nw,D0,slope=-1.65,intercept=6.5):
    """Methodology used in Bringi et al. (JAOT 2009).
  This techniques separates convective and stratiform based upon 
   Normalized gamma intercept parameter - Mean Drop Diameter, 
   log10(Nw)-D0 space.
   This is a physically based separation technique based upon DSD retrieved
   from dual-frequency profiler data.  
  A fit line is found in the form log10(Nw) = C1 * D0 + C2.
  Values of log10(Nw) above the line are considered convective.
  Values of log10(Nw) below the line are considered stratiform.
 INPUT::
  Nw            = Normalized gamma intercept parameter [mm^-1 m^-3]
  D0            = Volume-weighted Mean diameter [mm]
  slope         = Slope, C1 constant as defined above [unitless]
  intercept     = Intercept, C2 constant as define above [unitless]
 OUTPUT::
  TS            = Array identifying Conv-Strat Classification
                  1 = Stratiform
                  2 = Convective
 USAGE::
  TS = Nw_D0_cs(Nw,D0,[slope=-1.65],[intercept=6.5])
 NOTES::
  Default values:
   slope     = -1.65 Midpoint of range given in Thurai et al. (2010)
   intercept = 6.5   Midpoint of range given in Thurai et al. (2010)
    """
#---------------------------------------
    # Create the variable to pass from function
    TS = np.zeros_like(Nw)
    # Find the number of points in the time series
    nPts = len(Nw)

    TS = np.where(np.log10(Nw) >= (slope * D0 + intercept),2,1)

    return TS
#====================================================

def Nw_D0_cst(Nw,D0,slope=-1.65,intercept=6.5,Cthresh=0.1,Sthresh=-0.1):
    """Methodology put forth in Bringi et al. (JAOT 2009), 
   discussed in Thurai et al (JAOT 2010).
  This techniques separates convective and stratiform based upon 
   Normalized gamma intercept parameter - Mean Drop Diameter, 
   log10(Nw)-D0 space.

  This techniue calculates an index: Index = log10(NwData) - log10(NwLine).

  A fit line is found in the form log10(Nw) = C1 * D0 + C2.
   Values in a range about the fit line, specified by thresholds are 
    considered transition, while outside these thresholds are either 
    convective or stratiform as below .
  
  Values above log10(Nw)+Cthresh are considered convective.
  Values below log10(Nw)-Sthresh are considered stratiform.
 INPUT::
  Nw            = Normalized gamma intercept parameter [mm^-1 m^-3]
  D0            = Volume-weighted Mean diameter [mm]
  slope         = Slope, C1 constant as defined above [unitless]
  intercept     = Intercept, C2 constant as define above [unitless]
  Cthresh       = Threshold added to fit line, above which considered Con
  Sthresh       = Threshold added to fit line, below which considered Sf
 OUTPUT::
  TS            = Array identifying Conv-Strat Classification
                  1 = Stratiform
                  2 = Convective
                  3 = Transition
 USAGE::
  TS = Nw_D0_cs(Nw,D0,[slope=-1.65],[intercept=6.5],[Cthresh=0.2],[Sthresh=-1.])
 NOTES::
  Default values:
   slope     = -1.65 Midpoint of range given in Thurai et al. (2010)
   intercept = 6.5   Midpoint of range given in Thurai et al. (2010)
   Cthresh   = 0.1   Used in Thurai et al. (2010)
   Sthresh   = 0.1   Used in Thurai et al. (2010)
    """
#---------------------------------------
    # Create the variable to pass from function
    TS = np.zeros_like(Nw)
    # Find the number of points in the time series
    nPts = len(Nw)
    
#    # Mask any zero values of Nw to limit math errors
#    np.ma.masked_equal(Nw,0.)

    # Create index line
    ND0_Index = np.ma.log10(Nw) - (slope * D0 + intercept)

    # Search for stratiform
    TS[ND0_Index <= Sthresh] = 1
# = np.where(ND0_Index <= Sthresh,1,0)

    # Search for convective
    TS[ND0_Index >= Cthresh] = 2

    # Search for transition
    TS[(ND0_Index > Sthresh) & (ND0_Index < Cthresh)] =  3
#    TS[ = np.where(np.log10(Nw) >= (slope * D0 + intercept),2,1)

    return TS
#====================================================
