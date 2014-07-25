"""
muphyspy.io.write_netcdf
=========================

A grouping of scripts designed to create a NetCDF file (NetCDF4).


Created by Nick Guy.

"""
# HISTORY::
#  19 May 2014 - Nick Guy.   NRC, NOAA/NSSL (nick.guy@noaa.gov)
# FUNCTIONS::
# rab_netcdf_vars - Pulls out a number of variables from NetCDF file created
#                   using R Black results.
#-------------------------------------------------------------------
# Load the needed packages
from netCDF4 import Dataset
#===============================================================
# BEGIN FUNCTIONS
#===============================================================
def create_file(FilePathName,format='NetCDF4'):
    """Open a write-enabled NetCDF file.
 INPUT::
  FilePathName    = Long string path to NetCDF file to be written
        OPTIONAL
  format          = File formatDefault to NetCDF4
  StartT          = Start time for subsetting [matlab date number instance]
  EndT            = End time for subsetting [matlab date number instance]
 OUTPUT::
  ncF             = Output file pointer
 USAGE::
  ncF = create_file(filepathString,FillVal)
    """
# HISTORY::
#  19 May 2014 - Nick Guy NOAA/NSSL/WRDD, NRC
# NOTES::
#---------------------------------------
    # Open the file
    ncf = Dataset(FilePathName, 'w', format=format)
    return
#===============================================================
#def set_dims(dims,labs):
#    # Check to see if dimension of dimensions and labels are the same
#    if len(dims) != len(labs):
##        print 'Number of labels must be the same as the number of dimensions'
#        return
#    else:
#        for nn in len(dims):
            
    