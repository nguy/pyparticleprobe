"""
pyparticleprobe.dsd_calcs.zr
===============================

A grouping of functions for calculations of a Z-R relationship from a drop
size distribution.

Adapted by Nick Guy.

"""
# HISTORY::
#   28 Feb 2014 - Nick Guy. NRC, NOAA/NSSL (nick.guy@noaa.gov)
#                 Converted NCL functions below to Python 
# FUNCTIONS::
# linreg - Least squares linear regression fit
# regfit_powerlaw - Create a fit line from the regression data
# SD_filter - Filter a variable given a standard deviation
# regfit_abcd - Calculate a,b,c,d power law coefficients from regression data
# save_stats - Save a number of statistics into a text file for documentation
# get_zr_linreg - Solve for a,b,c,d coefficient and exponents in ZR relationship
#-------------------------------------------------------------------
# Load the needed packages
import numpy as np
from scipy import stats
#-------------------------------------------------------------------
# Define various constants that may be used for calculations
#
#===============================================================
# BEGIN FUNCTIONS
#**===============================================================

def linreg(Z,R):
    """Calculate the linear regression of two variables
 INPUT::
  Z               = Reflectivity [mm^6 m^-3]
  R               = Rainfall rate [mm h^-1]
 OUTPUT::
  slope           = Slope of the regression line
  intercept       = Intercept of the regression line
  rVal            = Correlation coefficient
  pVal            = Two-sided P-value test (null hypothesis slope = 0)
  std_err         = Standard Error of estimate
 USAGE::
  slope, intercept, r_value, p_value, std_err = linreg(x,y)
 NOTES::
  The masked stats function is used, otherwise the regression is performed on
    all data.
  Note that both variables are put in log-base10 space, to keep linear and because
    the vast majority of data is generally weighted to lower values.
    """
#---------------------------------------
# Use the Scipy linear regression algorithm to calculate
    slope, intercept, rVal, pVal, std_err = stats.mstats.linregress(np.log10(Z),np.log10(R))

    return slope, intercept, rVal, pVal, std_err
#**====================================================

def regfit_powerlaw(Z,slope,intercept,limLo=1E-2,limHi=1E6,Rtest=False):
    """Calculate a fit line to the linearly regressed data by 
  Create a linear array of independent reflectivity values.  Calculate the dependent
   Rainfall rate array using a power law distribution.

  Optionally return a "test" array calculating the rainfall from the actual
   reflectivity values via power law distribution.  
   This can then be used (with the SD_filter function).
   to remove outlier data. 
 INPUT::
  Z               = Reflectivity [mm^6 m^-3]
  slope           = Slope of regression line
  intercept       = Intercept of the regression line
       OPTIONAL
  limLo           = Lower limit of line
  limHi           = Upper limit of line
  Rtest           = Set True to return Rcalc "test" array
 OUTPUT::
  Zfit            = Fit line of Z as the independent variable
  Rfit            = Fit line of R as the dependent variable
  Rcalc           = Optionally returned array of rainfall rate based upon input Z
 USAGE::
  Zfit, Rfit, [Rtest] = linreg(x,y)
 NOTES::
  The masked stats function is used, otherwise the regression is performed on
    all data.
    """
#---------------------------------------
    # Create the fit lines
    Zfit = np.linspace(limLo,limHi,len(Z))
    Rfit = (10.**intercept) * (Zfit**slope)
    
    if Rtest:
        Rcalc = (10.**intercept) * (Z**slope)
        return Zfit,Rfit,Rcalc
    else:
        return Zfit,Rfit
#**====================================================

def SD_filter(Var,R,Rfit,Xsd):
    """Applies a filter to data at each point using the standard deviation
   of the entire data series.
 INPUT::
  Var             = Variable to be filtered
  R               = Rainfall rate [mm h^-1]
  Xsd             = Multiplication factor of Std Dev
 OUTPUT::
  VarFilt         = Filtered input variable (masked)
 USAGE::
  VarOut = SD_filter(VarIn,R,Xsd)
 NOTES::
  This module assumes that the variable to be filtered is the same
    dimensionality as the Rainfall rate variable.
    """
#---------------------------------------
    # Find the standard deviation of scatter to remove outliers
    sigDev = Xsd * R.std()
    
    # Create the array for the filtered data
    VarFilt = Var.copy()
    
    # Apply the condition to mask the data
    VarFilt = np.ma.masked_where((R <= (Rfit-sigDev)) | (R >= (Rfit+sigDev)),VarFilt,
                                   copy=False)
                                   
    return VarFilt
#**====================================================

def regfit_abcd(slope,intercept):
    """Calculate the a, b, c, d coefficients give the interept and slope, 
   assuming that reflectivity is the dependent variable 
 INPUT::
  slope           = Slope of regression line
  intercept       = Intercept of the regression line
  
 OUTPUT::
  a               = a coefficient in Z = aR^b power law
  b               = b coefficient in Z = aR^b power law
  c               = c coefficient in R = cZ^d power law
  d               = d coefficient in R = cZ^d power law
 USAGE::
  a,b,c,d = fit_abcd(slope,intercept)
 NOTES::
   This method assumes that reflectivity (Z) was the independent variable and
    Rainfall rate (R) was the dependent variable during linear regression.
    """
#---------------------------------------
    # Calculate a, b, c, and d coefficients
    c = 10.**intercept
    d = slope
    a = (1./c)**(1./d)
    b = (1./d)
    
    return a,b,c,d
#**====================================================

def get_zr_linreg(Z,R,filter=False,SDmult=1.,limLo=1E-2,limHi=1E6):
    """Use linear regression to solve find a Z-R relationship
 INPUT::
  Z               = Reflectivity [mm^6 m^-3]
  R               = Rainfall rate [mm h^-1]
     OPTIONAL::
  filter          = Set True to also return filtered a,b,c,d
  SDmult          = Multiplier for Standard deviation filter (if 3; then 3 * Std Dev)
  limLo           = Lower limit of line
  limHi           = Upper limit of line
  See the printout for details of inputs
 OUTPUT::
  a               = a coefficient in Z = aR^b power law
  b               = b exponent in Z = aR^b power law
  c               = c coefficient in R = cZ^d power law
  d               = d exponent in R = cZ^d power law
 USAGE::
  zr.save_stats(fname,[**args])
    """
#---------------------------------------
    # Calculate a least squares linear regression fit of log-log distribution
    Regslp, RegInt, rVal, pVal, stdErr = linreg(Z,R)
    
    # Assign values from linear regression for coefficients
    a_all,b_all,c_all,d_all = regfit_abcd(Regslp,RegInt)

    # Apply a filter if requested
    if filter:
        # Line fits for independent Z (linear array) and dependent R via power law relationship
        # The Rtest = true returns the test array (power law) for next filtering step 
        Zfit, RegFit, RRtest = regfit_powerlaw(Z,Regslp,RegInt,limLo=limLo,limHi=limHi,Rtest=True)

        #mask = (R >= (RRtest-sigDev)) & (R <= (RRtest+sigDev))
        # Filter the arrays within specified std deviation
        Z_filt = SD_filter(Z,R,RRtest,SDmult)
        R_filt = SD_filter(R,R,RRtest,SDmult)

        # Calculate the least squares linear regression fit of log-log distribution of filtered data
        RegslpFilt, RegIntFilt, rValFilt, pValFilt, stdErrFilt = linreg(Z_filt,R_filt)

        # Create an array for line fit using power law for filtered data
        ZfitFilt, RegFitFilt = regfit_powerlaw(Z,RegslpFilt,RegIntFilt,Rtest=False)
        del ZfitFilt

        a_filt,b_filt,c_filt,d_filt = regfit_abcd(RegslpFilt,RegIntFilt)
        
        # Find the number of elements in the filtered array
        nPts = R_filt.count()
        
        Info = 'Filtered'
        a, b, c, d, nPts = a_filt, b_filt, c_filt, d_filt, nPts      
        
        return a_filt,b_filt,c_filt,d_filt,nPts
    else:
        # Find the number of elements in the array
        nPts = len(R)
        
        Info = 'NonFiltered'
        a, b, c, d, nPts = a_all, b_all, c_all, d_all, nPts
        
    # Create a dictionary to transfer the data
    data = {'Mode' : Info,
            'a' : a,
            'b' : b,
            'c' : c,
            'd' : c,
            'number_pts' : nPts
            }
    return data
#**====================================================

def save_stats(fname,title=None,Conc=None,nPtsAll=None,cFactAll=None,dFactAll=None,aFactAll=None,
               bFactAll=None,nPtsFilt=None,cFactFilt=None,dFactFilt=None,aFactFilt=None,
               bFactFilt=None,rValAll=None,pValAll=None,stdErrAll=None,rValFilt=None,
               pValFilt=None,stdErrFilt=None,Nw=None,
               D0=None,Nw_D0_cst=None,W=None):
    """Save a text file with output stats calculated from Z-R relationship calculations
 INPUT::
  fname           = Name out output file
  title           = Title information to identify statistics
     OPTIONAL::
  See the printout for details of inputs
 OUTPUT::
  fname           = Text file
 USAGE::
  zr.save_stats(fname,[**args])
    """
#---------------------------------------
    # Create a single element needed to save file
    empty = [0.]

    ZRstatsTex = "=======================================\n"
    ZRstatsTex += "**** "+fname+" ****\n"
    ZRstatsTex += "DROP SIZE DISTRIBUTION CONCENTRATION\n"
    ZRstatsTex += "min Conc = " + str(Conc.min())+"\n"
    ZRstatsTex += "max Conc = " + str(Conc.max())+"\n"
    ZRstatsTex += "=======================================\n"
    ZRstatsTex += " \n"
    ZRstatsTex += title+"\n"
    ZRstatsTex += "========================================\n"
    ZRstatsTex += "PREFACTOR AND EXPONENT ESTIMATION\n"
    ZRstatsTex += "All Data: R=cZ^d:: c = "+str(cFactAll)+"  , d = "+str(dFactAll)+"\n"
    ZRstatsTex += "          Z=aR^b:: a = "+str(aFactAll)+"  , b = "+str(bFactAll)+"\n"
    ZRstatsTex += "          Correlation = "+str(rValAll)+" p = "+str(pValAll)+" StdErr = "+str(stdErrAll)+"\n"
    ZRstatsTex += "          # Points = "+str(nPtsAll)+"\n"
    ZRstatsTex += "          -----------------\n"
    ZRstatsTex += "Filtered Data: R=cZ^d:: c = "+str(cFactFilt)+"  , d = "+str(dFactFilt)+"\n"
    ZRstatsTex += "               Z=aR^b:: a = "+str(aFactFilt)+"  , b = "+str(bFactFilt)+"\n"
    ZRstatsTex += "               Correlation = "+str(rValFilt)+" p = "+str(pValFilt)+" StdErr = "+str(stdErrFilt)+"\n"
    ZRstatsTex += "          # Points = "+str(nPtsFilt)+"\n"
    ZRstatsTex += "=========================================\n"
    ZRstatsTex += " \n"

    ZRstatsTex += "==============================================\n"
    ZRstatsTex += " \n"
    ZRstatsTex += "Bringi et al. 2009, Conv-Strat-Trans\n"
    ZRstatsTex += "Stratiform: "+str(len(Nw_D0_cst[Nw_D0_cst == 1]))+" points"+"\n"
    ZRstatsTex += "Mean Nw = "+str(np.log10(Nw[Nw_D0_cst == 1].mean()))+",  SD = "+str(np.log10(Nw[Nw_D0_cst == 1].std()))+"\n"
    ZRstatsTex += "Mean D0 = "+str(D0[Nw_D0_cst == 1].mean())+",  SD = "+str(D0[Nw_D0_cst == 1].std())+"\n"
    ZRstatsTex += "================================\n"
    ZRstatsTex += "Convective: "+str(len(Nw_D0_cst[Nw_D0_cst == 2]))+" points"+"\n"
    ZRstatsTex += "Mean Nw = "+str(np.log10(Nw[Nw_D0_cst == 2].mean()))+",  SD = "+str(np.log10(Nw[Nw_D0_cst == 2].std()))+"\n"
    ZRstatsTex += "Mean D0 = "+str(D0[Nw_D0_cst == 2].mean())+",  SD = "+str(D0[Nw_D0_cst == 2].std())+"\n"
    ZRstatsTex += "=================================\n"
    ZRstatsTex += "Transition: "+str(len(Nw_D0_cst[Nw_D0_cst == 3]))+" points"+"\n"
    ZRstatsTex += "Mean Nw = "+str(np.log10(Nw[Nw_D0_cst == 3].mean()))+",  SD = "+str(np.log10(Nw[Nw_D0_cst == 3].std()))+"\n"
    ZRstatsTex += "Mean D0 = "+str(D0[Nw_D0_cst == 3].mean())+",  SD = "+str(D0[Nw_D0_cst == 3].std())+"\n"
    ZRstatsTex += "=================================\n"
#    ZRstatsTex += "Mean W = "+str(W[Nw_D0_cs == 2].mean())+",  SD = "+str(W[Nw_D0_cs == 2].std())+"\n"
    ZRstatsTex += "==============================================\n"
    ZRstatsTex += " \n"
    
    # Save the file
    np.savetxt(fname, empty, header=ZRstatsTex)
#====================================================

