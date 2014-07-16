# -*- coding: utf-8 -*-
"""
ReadFile.py - Class for reading Data from file
"""
import io.read_p3_2dimage_probe as ReadData
import os
from DropSizeDistribution import DropSizeDistribution

import numpy as np
########################
## BEGIN MAIN CODE
########################

def read_file(filename, file_format=None, Instr=None, flCDF=None, \
              Subset=False, SubType=None, Smin=None, Smax=None, \
              AvDat=True, RunAv=6,\
              ):
    '''
    Takes a filename pointing to a probe data file and returns
    a drop size distribution object.

    Usage:
    data = read_file(filename)

    Returns:
    FileReader object
    
    '''
    
    reader = FileReader(filename, file_format=file_format, Instr=Instr, flCDF=flCDF, \
                        Subset=Subset, SubType=SubType, Smin=Smin, Smax=Smax, \
                        AvDat=AvDat, RunAv=RunAv)
    
    dsd = DropSizeDistribution(reader)
    
    return dsd



class FileReader(object):

    '''
    FileReader class to process data files.  
    '''

    def __init__(self, filename, file_format=None, Instr=None, flCDF=None, \
                        Subset=False, SubType=None, Smin=None, Smax=None, \
                        AvDat=False, RunAv=None,\
                        ):
        self.filename = filename
        self.fileformat =  file_format
        self.flightlevel_filename = flCDF
        
        yyyy = os.path.basename(filename).split(".")[1][0:4]
        mm = os.path.basename(filename).split(".")[1][4:6]
        dd = os.path.basename(filename).split(".")[1][6:8]
        
        if file_format is None:
            file_format = 'ucsc_netcdf'
        if Instr is None:
            Instr = 'pip'
        
        if file_format is 'rab_netcdf':
            sizebins, Time, Conc, rhoair, Wvel, Alt = ReadData.rab_netcdf_vars(filename,\
              Subset=Subset, SubType=SubType, Smin=Smin, Smax=Smax)
            
            self.sizebins = probe['Sizebins']
            self.time = probe['Time']
            self.conc = probe['Conc']
            self.rhoair = probe['RhoAir']
            self.w_vel = probe['w_vel_air']
        
        if file_format is 'ucsc_netcdf':
            probe = ReadData.ucsc_netcdf_vars(filename, \
              Instr, yyyy, mm, dd, flCDF, \
              Subset=Subset, SubType=SubType, Smin=Smin, Smax=Smax, \
              AvDat=AvDat, RunAv=RunAv)
            
            self.sizebins = probe['Sizebins']
            self.time = probe['Time']
            self.conc = probe['Conc']
            self.rhoair = probe['RhoAir']
            self.w_vel_air = probe['W_vel_air']
            self.alt = probe['Altitude']
            self.bin_edges = probe['Bin_edges']
            
        del probe,yyyy,mm,dd
            
#        del sizebins,Time,Conc,rhoair,Wvel,Alt,yyyy,mm,dd
        
            
        
