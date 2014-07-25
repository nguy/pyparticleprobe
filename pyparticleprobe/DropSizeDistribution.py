# -*- coding: utf-8 -*-
"""
DropSizeDistribution.py - Class for DSD parameter calculations
"""
from dsd_calcs.moments import *
from dsd_calcs.params import get_drops_params, get_timeseries_params, get_model_params, get_bins_params
import dsd_calcs.partition as CSpart
from dsd_calcs import axis_ratio
from dsd_calcs.zr import get_zr_linreg 

import inspect
from copy import deepcopy
import numpy as np

import pytmatrix
from pytmatrix.tmatrix import Scatterer
from pytmatrix.psd import PSDIntegrator
from pytmatrix import orientation, radar, tmatrix_aux, refractive

########################
## BEGIN MAIN CODE
########################

class DropSizeDistribution(object):

    '''
    DropSizeDistribution class to hold DSD's and calculate parameters
    and relationships. Should be returned from the disdrometer*reader objects.
    '''

    def __init__(self, reader):

        """
        DropSizeDistribution class to hold DSD's and calculate parameters
        and relationships. Should be returned from the disdrometer*reader objects.
        
        Parameters returned:
        filename : Filename of probe data processed [string]
        flightlevel_filename : Filename of flight level data processed [string]
        sizebins : Mid-point size of probe bins [micron]
        time : Time [matlab date object, fractional day]
        conc : Concentration of particles [m^-3]
        rhoair : Density of air [kg/m^3]
        w_vel_air : Vertical velocity [m/s]
        alt : Aircraft altitude [m]
        bin_edges : Bin edges [micron], Note one element greater than sizebins
        """
    
        if (inspect.isclass(type(reader)) or inspect.isclass(reader)):
#            dsd = deepcopy(reader)
            self.filename = reader.filename
            self.flightlevel_filename = reader.flightlevel_filename
            self.sizebins = reader.sizebins
            self.time = reader.time
            self.conc = reader.conc
            self.rhoair = reader.rhoair
            self.w_vel_air = reader.w_vel_air
            self.alt = reader.alt
            self.bin_edges = reader.bin_edges
            
            lt = len(reader.time)
            self.Zh = np.zeros(lt)
            self.Zdr = np.zeros(lt)
            self.Kdp = np.zeros(lt)
            self.Ai = np.zeros(lt)
        else:
            #Initializes class instance but leaves it to other methods to
            #populate the class attributes.
            return
            
        # Find which dimensions are the time and bin dimensions
        for ndim in range(len(self.conc.shape)):
            if self.time.shape[0] == self.conc.shape[ndim]:
                self.timeDim = ndim     # Set the time dimension
            if self.sizebins.shape[0] == self.conc.shape[ndim]:
                self.binDim = ndim     # Set the bin dimension

    ##########################################
        
    def get_ts_parameters(self):
        '''
        Calculates the time series parameters by summing droplets at
        each time interval.  These are equivalent to disdrometer observations.
        
        Parameters returned:
        Ze : Reflectivity [mm^6 / m^3]
        dBZ : Reflectivity [dBZ]
        RR : Rainfall rate [mm / h]
        Dm : Mean mass-weighted diameter [m]
        Nw : Normalized gamma intercept [m^-1 m^-3]
        Nt : Total concentration per time [1 / m^3]
        LWC : Liquid water content [g/m^3]
        '''
        
        # Convert the DSD to m (/1E6.) and make 2D for computations
        Diam2D,RhoAir2D = np.meshgrid(self.sizebins[:]*1E-6, self.rhoair)
        
        ts_params = get_timeseries_params(self.conc, Diam2D, RhoAir2D, \
                                       sumaxis = self.binDim)
                                       
        self.Ze = ts_params['Ze']
        self.dBZ = ts_params['dBZ']
        self.RR = ts_params['RR']
        self.Dm = ts_params['Dm']
        self.Nw = ts_params['Nw']
        self.Nt = ts_params['Nt']
        self.LWC = ts_params['LWC']
        
        del ts_params

    ##########################################
        
    def get_bin_parameters(self):
        '''
        Calculates the time series parameters by summing droplets at
        each time interval.  These are equivalent to disdrometer observations.
        
        Parameters returned:
        Nd : Total concentration per bin [1 / m^3]
        Nd_Avg : Average concentration per bin [1 / m^3]
        Nd_Med : Median concentration per bin [1 / m^3]
        '''
        
        bin_params = get_bins_params(self.conc, sumaxis = self.timeDim)
                                       
        self.Nd = bin_params['Nd']
        self.Nd_Avg = bin_params['Nd_Avg']
        self.Nd_Med = bin_params['Nd_Med']
        
        del bin_params

    ##########################################
        
    def get_model_parameters(self, modelName='gamma_ua98'):
        '''
        Calculates model fit parameters for particle distributions.
        The different models supported are:
        'gamma_ua98' : Gamma model as in Ulbrich and Atlas 1998
        'gamma_ts96' : Gamma model as in Tokay and Short 1996
        'exp'        : Exponential model as in Waldvogel 1974
        
        The default is 'gamma_ua98'. 
        
        Parameters returned:
        Shape : Shape parameter [m]
        D0 : Median volume diameter [m]
        Slope : Slope parameter [m^-1]
        N0 : Intercept parameter [m^(-1-shape) m^-3]
        '''
        
        # The variables are created before because if exponential
        # is chosen, the shape and D0 remain empty due to exponential
        # being a two parameter solution while the gamma is 3 parameter.
        self.shape = []
        self.slope = []
        self.N0 = []
        self.D0 = []
    
        # Convert the DSD to m (/1E6.) and make 2D for computations
        Diam2D,RhoAir2D = np.meshgrid(self.sizebins[:]*1E-6, self.rhoair)
        
        self.dsd_model = modelName
        
        mod_params = get_model_params(self.conc, Diam2D, RhoAir2D, \
                                modelName=modelName, sumaxis=self.binDim)
                                          
        self.shape = mod_params['shape']
        self.D0 = mod_params['D0']
        self.slope = mod_params['slope']
        self.N0 = mod_params['N0']
        
        del mod_params

    ##########################################
        
    def get_drop_paramaters(self):
        '''
        Calculates parameters for each droplet.
        This returns a dictionary of values for:
        
        Parameters returned:
        Ze : Equivalent radar reflectivity [linear units]
        dBZe : Equivalent radar reflectivity [log units]
        RR : Rainfall rate [mm /h]
        water_mass : Estimated water mass of droplet [g]
        '''
    
        # Convert the DSD to m (/1E6.) and make 2D for computations
        Diam2D,RhoAir2D = np.meshgrid(self.sizebins[:]*1E-6, self.rhoair)
        
        
        
        drop_params = get_drops_params(self.conc, Diam2D, RhoAir2D)
        
        self.Ze_drop = drop_params['Ze']
        self.dBZe_drop = drop_params['dBZe']
        self.RR_drop = drop_params['rain_rate']
        self.water_mass = drop_params['water_mass']
        
        del drop_params

    ##########################################
        
    def calc_radar_parameters(self, wavelength=tmatrix_aux.wl_X):
        '''
        Calculates the radar parameters and stores them in the object.
        Defaults to X-Band wavelength and Thurai et al. 2007 axis ratio setup.
        Sets object radar parameters:
        Zh, Zdr, Kdp, Ai

        Parameter:
        wavelength = tmatrix supported wavelength.
        '''
        self._setup_scattering(wavelength)
        print self.bin_edges/1000.
        print self.Nt.min(),self.Nt.max()

        for t in range(0, len(self.time)):
            BinnedDSD = pytmatrix.psd.BinnedPSD(self.bin_edges/1000., self.Nt[t]/1E9)
            self.scatterer.psd = BinnedDSD
            self.scatterer.set_geometry(tmatrix_aux.geom_horiz_back)
            self.Zdr[t] = 10 * np.log10(radar.Zdr(self.scatterer))
            self.Zh[t] = 10 * np.log10(radar.refl(self.scatterer))
            self.scatterer.set_geometry(tmatrix_aux.geom_horiz_forw)
            self.Kdp[t] = radar.Kdp(self.scatterer)
            self.Ai[t] = radar.Ai(self.scatterer)

    ##########################################

    def _setup_scattering(self, wavelength):
        self.scatterer = Scatterer(wavelength=wavelength,
                                   m=refractive.m_w_10C[wavelength])
        self.scatterer.psd_integrator = PSDIntegrator()
        self.scatterer.psd_integrator.axis_ratio_func = lambda D: 1.0 / \
            axis_ratio.axis_ratio_THBRS07(D)
        self.scatterer.psd_integrator.D_max = 8.0
        self.scatterer.psd_integrator.geometries = (
            tmatrix_aux.geom_horiz_back, tmatrix_aux.geom_horiz_forw)
        self.scatterer.or_pdf = orientation.gaussian_pdf(20.0)
        self.scatterer.orient = orientation.orient_averaged_fixed
        print "MADE IT CC!!!!!!!!"
        self.scatterer.psd_integrator.init_scatter_table(self.scatterer)
        print "MADE IT HERE!!!!!!!!"
    ##########################################
    
    def _value_range_check(self):
        '''
        Print some values to the terminal to check the processed ranges.
        '''
        
        print "=================================================="
        print "---- Water----"
        print "**** "+self.filename+" ****"
        print "LWC Min: %f, Max: %f, Avg: %f, std: %f" % \
               (self.LWC.min(),self.LWC.max(),self.LWC.mean(),self.LWC.std())
        print "Dm Min: %f, Max:  %f, Avg: %f, std: %f" % \
               (self.Dm.min(),self.Dm.max(),self.Dm.mean(),self.Dm.std())
        print "Shape Min: %f, Max:  %f, Avg: %f, std: %f" % \
               (self.shape.min(),self.shape.max(),self.shape.mean(),self.shape.std())
        print "Slope Min: %f, Max:  %f, Avg: %f, std: %f" % \
               (self.slope.min(),self.slope.max(),self.slope.mean(),self.slope.std())
        print "log10(N0) Min: %f, Max:  %f, Avg: %f, std: %f" % \
               (np.ma.log10(self.N0.min()),np.ma.log10(self.N0.max()),np.ma.log10(self.N0.mean()),np.ma.log10(self.N0.std()))
        print "RR Min: %f, Max:  %f, Avg: %f, std: %f" % \
               (self.RR.min(),self.RR.max(),self.RR.mean(),self.RR.std())
        print "  "
        print "D0 Min: %f, Max:  %f, Avg: %f, std: %f" % \
               (self.D0.min(),self.D0.max(),self.D0.mean(),self.D0.std())
        print "log10(Nw) Min: %f, Max:  %f, Avg: %f, std: %f" % \
               (np.ma.log10(self.Nw.min()),np.ma.log10(self.Nw.max()),np.ma.log10(self.Nw.mean()),np.ma.log10(self.Nw.std()))
        
        print "=================================================="
        print "Mean W = ",self.w_vel_air.mean(),", SD = ",self.w_vel_air.std()

    ##########################################
    
    def conv_strat_sep(self, method=None, window=None, \
                       ThreshR=None, ThreshSD=None, WThresh=None,\
                       NwD0_slope=None, NwD0_int=None, cThresh=None, sThresh=None, lOut=True):
        '''
        Apply a convective-stratiform-(transition) separation algorithm to the time
        series data.
        
        Parameters:
        method : Str
            Method to use for calculation.  The choices for method are:
                'b09' -  Bringi et al. (2009) uses a Nw-D0 separator line and adds 
                a 'transition' category according to Thurai et al. (2010).
                't01' - Testud et al. (2001) uses rainfall magnitude, with separation at
                a threshold value of 10 mm/h.
                'b03' - Bringi et al. (2003) uses the standard deviation of rainfall with
                separation at a threshold of 1.5 mm/h.
                'i12' - Islam et al. (2012) uses both rainfall and magnitude combining the 
                previous two methods.
                'a00' - Atlas et al. (2000) uses the vertical velocity magnitude of air to 
                find updrafts and downdraft areas, with a
        window : int
            The number of points in which the window calculation (e.g. averaging) is performed.
            Only used in some methods.
        ThreshR : float
            Rain rate threshold to use as separator.
            Only used in some methods.
        ThreshSD : float
            Rain rate standard deviation threshold to use as separator.
            Only used in some methods.
            
        Output:
        An array the same length as input is created.  The results are recorded as:
        1 = Stratiform
        2 = Convective
        3 = Transition
        
        '''
        
        if method is None:
            method = 'b09'
        else:
            method = method
        
        # Create arrays to hold C/S identifier (1=Sf, 2=Con, 3=Trans)
#B09_cst = np.zeros_like(Time)

        # Create a padded array and apply the mirroring technique for end points
        #  Caveat should be noted that if variability exists at end points 
        #  (convective cell) it may be  "watered down" or amplified.
        RR_pad = np.pad(self.RR,5,'reflect')

        if method == 't01':  #----Testud et al. 2001 method
            if ThreshR is None:
                ThreshR = 10.
            if window is None:
                window = 3 # This corresponds to a distance travelled by aircraft
            CST = CSpart.Rmag_cs(self.RR, Thresh=ThreshR, window=window)
        
        if method == 'b03':  #----Bringi et al. 2003 method
            if ThreshSD is None:
                ThreshSD = 1.5
            if window is None:
                window = 4 # This corresponds to a distance travelled by aircraft
            CST = CSpart.Rsd_cs(self.RR, Thresh=ThreshSD, window=window)
            
        if method == 'i12':  #----Islam et al. 2012 method 
            if ThreshR is None:
                ThreshR = 10.
            if ThreshSD is None:
                ThreshSD = 1.5
            if window is None:
                window = 4 # This corresponds to a distance travelled by aircraft
            CST = CSpart.RmagSD_cs(self.RR, ThreshR=ThreshR, \
                                    ThreshSD=ThreshSD, window=window)
        
        if method == 'a00':  #----Atlas et al. 2000 method
            if WThresh is None:
                WThresh = 1.3
            CST = CSpart.Wvel_cs(self.w_vel_air,Wthresh=WThresh)
            
        if method == 'b09':  #----Bringi et al 2009 method, Thurai et al. 2010 method
            if NwD0_slope is None:
                NwD0_slope = -1.65
            if NwD0_int is None:
                NwD0_int = 6.5
            if cThresh is None:
                cThresh = 0.2
            if sThresh is None:
                sThresh = -1.0
            CST = CSpart.Nw_D0_cst(self.Nw/1000.,self.D0*1000.,slope=NwD0_slope,intercept=NwD0_int,
                                Cthresh=cThresh,Sthresh=sThresh)
        return CST

    ##########################################
    
    def b09_sep_line(self, slope=None, intercept=None):
        '''
        Returns the separator line based upon slope and intercept used.
        If none is chosen, the values default to those in Bringi et al. (2009)
        '''
        if slope is None:
            slope = -1.65
        if intercept is None:
            intercept = 6.5

        # Create index line to use for separation
        B09_separator = slope * self.D0 * 1000. + intercept
        
        return B09_separator

    ##########################################
    
    def b09_get_index(self, slope=None, intercept=None):
        '''
        Returns the index for convective-stratiform-transition classification
        via the separator line based upon slope and intercept used.
        If none is chosen, the values default to those in Bringi et al. (2009)
        '''
        if slope is None:
            slope = -1.65
        if intercept is None:
            intercept = 6.5
            
        Nwtmp = (self.Nw.copy()) / 1000.

        # Create index line to use for separation
        ND0_Index = np.ma.log10(self.Nw/1000.) - (slope * self.D0 * 1000. + intercept)
                             
        # Mask values where D0 or Nw is equal zero.  Gives a false value of 1 due
        # to math of separator line
        np.ma.masked_where((self.D0 == 0.), ND0_Index)
        np.ma.masked_where((self.Nw == 0.), ND0_Index)
        
        return ND0_Index

    ##########################################
    
    def calc_R_Z_moment_relationship(self, filter=False, SDmult=1., \
                                     limLo=1E-2, limHi=1E6):
        '''
        Calculates the power law fit based upon moment derived parameters.
        Uses the linear regression technique, with optional filtering applied
        '''
        ZRrel = get_zr_linreg(self.Ze, self.RR, filter=filter, SDmult=SDmult, \
                             limLo=limLo, limHi=limHi)
                             
        return ZRrel

    ##########################################
    