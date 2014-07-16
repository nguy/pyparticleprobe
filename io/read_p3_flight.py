"""
pyaircraft.io.read_p3_flight
=========================

This is a grouping of scripts designed to process NOAA P-3 
 flight level data recorded during flights and put into NetCDF
 format by NOAA AOC.

Created by Nick Guy.

"""
# NOTES:: This has only been tested with DYNAMO data files, versions
#         may change and another function may be needed.
# HISTORY::
#   8 Jan 2014 - Nick Guy.   NRC, NOAA/NSSL (nick.guy@noaa.gov)   
#                Converted NCL functions below to Python
# FUNCTIONS::
#  flight_level_variable - Read in a variable from flight level NetCDF
#  flight_track - Read in data to for flight track
#-------------------------------------------------------------------
# Load the needed packages
from scipy.io import netcdf
from netCDF4 import Dataset,num2date
import json
import numpy as np
import pytz
#-------------------------------------------------------------------
# Define various constants that may be used for calculations
#===============================================================
# BEGIN FUNCTIONS
#===============================================================
def flight_level_variable(fname,Rec):
    """Read in data from NetCDF file containing P3 flight level data created
  by NOAA AOC.  The NetCDF should be read in the main program and passed
  to this function.
  A call such as this can be used in the main program:
      FltncID=addfile(FlightFileStringName,"r")
 INPUT::
  fname           = Filename [string]
  Rec             = Variable name to be pulled out [string]
 OUTPUT::
  VarOut          = Masked array containing variable data
 USAGE::
  Lat = read_flight_level_dynamo('P3.nc','LatGPS.3')
 NOTES::
 Data file structure::
  Available variables (not full list) :
  LonGPS.3      = Novatel GPS Longitude
  LatGPS.3      = Novatel GPS Latitude
  AltGPS.3      = Novatel GPS Altitude [m]
  THdgI-GPS.1   = True heading [deg]
  TRK.1         = Track [deg]
  AltPaADDU.1   = Pressure altitude [m]
  WSZ_DPJ       = Vertical wind via D Jorgensen calculation [m/s]
  TA.1          = Ambient Temperature [C]
  TD.1          = Dewpoint Temperature [C]
  TVIRT.1       = Virtual Temperature [K]
  THETA.1       = Potential Temperature [K]
  THETAE.1      = Equivalent Potential Temperature [K]
  THETAV.1      = Virtual Potential Temperature [K]
  WS.1          = Wind Speed [m/s]
  WD.1          = Wind Direction [deg]
  HUM_REL.1     = Relative Humidity [%]
  HUM_SPEC.1    = Specific Humidity [g/kg]
  MR.1          = Mixing ratio [g] [g/g?]
  EE.1          = Vapor Pressure [hPa]
  EW.1          = Saturated Vapor Pressure [hPa]
      """
# HISTORY::
#  25 Jul 2013 - Nick Guy NOAA/NSSL/WRDD, NRC
#   8 Jan 2014 - With change to python, the procedure of setting the bad
#                 data by user is no longer used, but the numpy masking
#                 technique is instead used to nuke "bad" occurences
#---------------------------------------------------
    # Read the NetCDF
    ncFile = netcdf.netcdf_file(fname,'r')

    # Get the variable of interest
    VarOut = ncFile.variables['Rec'][:]

    # Mask any "_FillValue" or some missing_data type attribute
    #try:
    #    VarOut = np.ma.masked_values(VarOut, VarOut.missing_value)
    #except:
    #    pass
    #try:
    #    VarOut = np.ma.masked_values(VarOut, VarOut.Missing_Value)
    #except:
    #    pass
    #try:
    #    VarOut = np.ma.masked_values(VarOut, VarOut._FillValue)
    #except:
    #    pass

    # Mask any NaN values
    #VarOut = np.ma.masked_values(VarOut, np.isnan(VarOut)

    return VarOut
#**====================================================
def flight_track(fname):
    """Read in data from NetCDF file containing P3 flight level data created
  by NOAA AOC.  Pull out the needed variables for flight track info.
 INPUT::
  fname           = Filename [string]
 OUTPUT::
  Lat             = Aircraft latitude
  Lon             = Aircraft longitude
  Alt             = Aircraft altitude
  PAlt            = Aircraft pressure altitude
  Time            = Aircraft time array
 USAGE::
  Lat,Lon,Alt,PAlt = flight_track(fname)
    """
# HISTORY::
#   7 Mar 2014 - Nick Guy NOAA/NSSL/WRDD, NRC
#---------------------------------------------------
    # Read the NetCDF
    ncFile = Dataset(fname,'r')
    
    # Pull out each variable
    Lat = ncFile.variables['LatGPS.3'][:]
    Lon = ncFile.variables['LonGPS.3'][:]
    Alt = ncFile.variables['AltGPS.3'][:]
    PAlt = ncFile.variables['AltPaADDU.1'][:]
    
    # Pull out the start time
    StartTime = ncFile.StartTime
    
    # Create a time array 
    TimeSec = np.linspace(StartTime,StartTime + len(Lat), len(Lat))
    
    Time_unaware = num2date(TimeSec,'seconds since 1970-01-01 00:00:00+0:00')
    Time = Time_unaware#.replace(tzinfo=pytz.UTC)
    
    # Now mask missing values
    np.ma.masked_invalid(Lat)
    np.ma.masked_invalid(Lon)
    np.ma.masked_invalid(Alt)
    np.ma.masked_invalid(PAlt)
    
    return Lat,Lon,Alt,PAlt,Time

