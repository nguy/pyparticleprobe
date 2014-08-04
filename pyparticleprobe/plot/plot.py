"""
muphyspy.plot
===============================

A grouping of functions for plotting microphysical analyses.

Adapted by Nick Guy.

"""
# HISTORY::
#    9 Mar 2014 - Nick Guy. NRC, NOAA/NSSL (nick.guy@noaa.gov)
#                 Adapted from other programs to modularize. 
# FUNCTIONS::
# NwD0_scatter - Nw-D0 scatterplot
# Nw_D0_density_hex - Scatter density plot of Nw vs. D0 using hex method
# Nw_D0_density_pdf - 2D density plot of Nw vs. D0
# zr_scatter - Z-R scatterplot
# zr_density - Z-R density scatterplot
# zr_density_hex - Scatter density plot of Z vs. R using hex method
# zr_density_pdf - 2D density plot of Z vs. R
# plotDate_ts - Time series plot
# plotHov - Hovmoeller plot
#-------------------------------------------------------------------
# Load the needed packages
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
from matplotlib import ticker as mtic
from  matplotlib.dates import DateFormatter#, MinuteLocator
from matplotlib.colors import LogNorm
import general.library as gl
#-------------------------------------------------------------------
# Define various constants that may be used for calculations
#
#===============================================================
# BEGIN FUNCTIONS
#===============================================================
def NwD0_scatter(ax,Nw,D0,title=' ',col='k',msize=20,
                 xlabFontSize=16,xpad=7,xMajspace=0.5,xMinspace=0.1,
                 ylabFontSize=16,ypad=7,yMajspace=0.5,yMinspace=0.1,
                 xlims=(0.,5.)):
    """Create a scatterplot of Nw vs. D0
 INPUT::
  ax              = Axis to plot against (produced externally)
  Nw              = Normalized intercept parameter [mm^-1 m^-3]
  D0              = Median volume diameter [mm]
     OPTIONAL
  title           = Plot title
  col             = Color of markers
  xlabFontSize    = Font size for x-axis labels
  xpad            = Padding to apply to x-axis label
  ylabFontSize    = Font size for y-axis labels
  ypad            = Padding to apply to y-axis label
 OUTPUT::
  p               = Plot
 USAGE::
  p = NwD0_scatter(ax,Nw,D0,**args)
    """
# HISTORY::
#  10 Mar 2014 - Nick Guy NOAA/NSSL/WRDD, NRC
#--------------------------------------------------------
    # create the scatter plot
    plt.scatter(D0,np.log10(Nw),c=col,s=msize,edgecolors='none')
    
    # First remove italics from mathtext then label axes
    plt.rc('mathtext',fontset='stixsans',default='regular') # Removes italics mathtext
    plt.xlabel(r'D$_0$ (mm)',labelpad=xpad,fontsize=xlabFontSize)
    plt.ylabel(r'log$_{10}$[N$_w$] (mm$^{-1}$ m$^{-3}$)',labelpad=ypad,fontsize=ylabFontSize)

    plt.xlim(xlims) # Set X-axis limits
    
    # Set the plot title
    ax.set_title(title) 
    ax.xaxis.set_major_locator(mtic.MultipleLocator(xMajspace)) # Set the major tickmarks 
    ax.xaxis.set_minor_locator(mtic.MultipleLocator(xMinspace)) # Set the minor tickmarks
    ax.yaxis.set_major_locator(mtic.MultipleLocator(yMajspace)) # Set the major tickmarks 
    ax.yaxis.set_minor_locator(mtic.MultipleLocator(yMinspace)) # Set the minor tickmarks
    
    return
#**===============================================================
def Nw_D0_density_hex(ax,Nw,D0,title=' ',grid=(36,61),minct=0.01,
                   xlabFontSize=16,xpad=7,ylabFontSize=16,ypad=7,
                   xlims=(0.1,3.5),ylims=(1.0,7.0)):
    """Create a density scatterplot using the hex bin method (color marker based on density)
 INPUT::
  ax              = Axis to plot against (produced externally)
  Nw              = Normalized intercept parameter [mm^-1 m^-3]
  D0              = Median volume diameter [mm]
     OPTIONAL
  title           = Plot title
  grid            = The bin size in (x,y) to be used
  minct           = Minimum count to be displayed
 OUTPUT::
  p               = Plot
 USAGE::
  p = NW_D0_density_hex(ax,Nw,D0,**args)
    """
# HISTORY::
#  10 Mar 2014 - Nick Guy NOAA/NSSL/WRDD, NRC
#--------------------------------------------------------
    # This will plot a hex bin plot (density scatter plot)
    # Establish an array from 0-100 to show colors as density
    c = np.ones_like(D0) * 100 / len(D0)
    
    # Create the hex plot
    plt.hexbin(D0,np.log10(Nw),C=c,gridsize=grid,mincnt=minct,
                   reduce_C_function=np.sum)

    plt.ylim(xlims) # Set X-axis limits
    plt.xlim(ylims) # Set Y-axis limits
    
    # First remove italics from mathtext then label axes
    plt.rc('mathtext',fontset='stixsans',default='regular') # Removes italics mathtext
    plt.xlabel(r'D$_0$ (mm)',labelpad=xpad,fontsize=ylabFontSize)
    plt.ylabel(r'log$_{10}$[N$_w$] (mm$^{-1}$ m$^{-3}$)',labelpad=ypad,fontsize=xlabFontSize)
    
    ax.set_title(title) # Set the plot title
    
    # Add colorbar
    cb = plt.colorbar(shrink=0.85)
    cb.set_label('Normalized Counts')
    
    del c
    
    return
#**===============================================================
def Nw_D0_density_pdf(ax,Nw,D0,title=None,grid=(36,61),minct=0.01,
                   xlabFontSize=16,xpad=7,ylabFontSize=16,ypad=7,
                   xlims=(0.1,3.5),ylims=(1.0,7.0)):
    """Create a density scatterplot using 2D histogram
 INPUT::
  ax              = Axis to plot against (produced externally)
  Nw              = Normalized intercept parameter [mm^-1 m^-3]
  D0              = Median volume diameter [mm]
     OPTIONAL
  title           = Plot title
  grid            = The bin size in (x,y) to be used
  minct           = Minimum count to be displayed
 OUTPUT::
  p               = Plot
 USAGE::
  p = Nw_D0_density_pdf(ax,Nw,D0,**args)
    """
# HISTORY::
#  10 Mar 2014 - Nick Guy NOAA/NSSL/WRDD, NRC
#--------------------------------------------------------
    # First calculate the 2D PDF
    pdf2_Nw_D0, xedges, yedges = gl.hist2d_masked(D0,np.log10(Nw),bins=grid,
                                                rng=(ylims,xlims),norm=True)
    
    # Replace any zero values with missing data for nice plots
    pdf2_Nw_D0 = np.ma.masked_equal(pdf2_Nw_D0,0)
    
    # Flip and rotate the 2D Histogram
    pdf2_Nw_D0 = np.rot90(pdf2_Nw_D0)
    pdf2_Nw_D0 = np.flipud(pdf2_Nw_D0)
    
    # Create 2D arrays from xedges, yedges
    X, Y = np.meshgrid(xedges, yedges)
    
    # Plot the data using colormesh
    plt.pcolormesh(X,Y,pdf2_Nw_D0)#,cmap=

    plt.ylim(xlims) # Set X-axis limits
    plt.xlim(ylims) # Set Y-axis limits
    
    # First remove italics from mathtext then label axes
    plt.rc('mathtext',fontset='stixsans',default='regular') # Removes italics mathtext
    plt.xlabel(r'D$_0$ (mm)',labelpad=xpad,fontsize=ylabFontSize)
    plt.ylabel(r'log$_{10}$[N$_w$] (mm$^{-1}$ m$^{-3}$)',labelpad=ypad,fontsize=xlabFontSize)
    
    ax.set_title(title) # Set the plot title
    
    # Add colorbar
    cb = plt.colorbar(shrink=0.85)
    cb.set_label('Normalized Counts')
    
    del pdf2_Nw_D0,xedges,yedges,X,Y
    
    return
#**===============================================================
def zr_scatter(ax,Z,R,title=None,col='k',msize=20,
               xlabFontSize=16,xpad=5,ylabFontSize=16,ypad=5,
               xlims=(1E-2,1E6),ylims=(1E-3,1E3)):
    """Create a scatterplot of Z vs. R
 INPUT::
  ax              = Axis to plot against (produced externally)
  Z               = Reflectivity factor [mm^6/m^3]
  R               = Rainfall rate [mm/h]
     OPTIONAL
  title           = Plot title
  col             = Color of markers
  size            = Marker size (0-15)
  xlabFontSize    = Font size for x-axis labels
  xpad            = Padding to apply to x-axis label
  ylabFontSize    = Font size for y-axis labels
  ypad            = Padding to apply to y-axis label
  xlims           = Limits to apply to x-axis
  ylims           = Limits to apply to y-axis
 OUTPUT::
  p               = Plot
 USAGE::
  p = zr_scatter(ax,Z,R,**args)
    """
# HISTORY::
#  10 Mar 2014 - Nick Guy NOAA/NSSL/WRDD, NRC
#--------------------------------------------------------
    # Create the scatter plot
    plt.scatter(Z,R,c=col,s=msize,edgecolors='none')
    plt.xscale('log') # Make X-axis logarithmic
    plt.yscale('log') # Make Y-axis logarithmic
    plt.xlim(xlims) # Set X-axis limits
    plt.ylim(ylims) # Set Y-axis limits
    
    # First remove italics from mathtext then label axes
    plt.rc('mathtext',fontset='stixsans',default='regular') # Removes italics mathtext
    plt.xlabel(r'Reflectivity (mm$^{6}$ m$^{-3}$)',labelpad=xpad,fontsize=xlabFontSize)
    plt.ylabel(r'Rainfall Rate (mm h$^{-1}$)',labelpad=ypad,fontsize=ylabFontSize)
    
    ax.set_title(title) # Set the plot title
    
    return
#**===============================================================
def zr_density_hex(ax,Z,R,title=None,grid=(80,60),minct=0.01,
                   xlabFontSize=16,xpad=7,ylabFontSize=16,ypad=7,
                   xlims=(-3,3),ylims=(-20.,60.)):
    """Create a density scatterplot using the hex bin method (color marker based on density)
 INPUT::
  ax              = Axis to plot against (produced externally)
  Z               = Reflectivity factor [mm^6/m^3]
  R               = Rainfall rate [mm/h]
     OPTIONAL
  title           = Plot title
  grid            = The bin size in (x,y) to be used
  minct           = Minimum count to be displayed
 OUTPUT::
  p               = Plot
 USAGE::
  p = zr_density_hex(ax,Z,R,**args)
    """
# HISTORY::
#  10 Mar 2014 - Nick Guy NOAA/NSSL/WRDD, NRC
# NOTES::
#--------------------------------------------------------
    # This will plot a hex bin plot (density scatter plot)
    # Establish an array from 0-100 to show colors as density
    c = np.ones_like(Z) * 100 / len(Z)
    
    # Create the hex plot
    plt.hexbin(10.*np.log10(Z),np.log10(R),C=c,gridsize=grid,mincnt=minct,
                   reduce_C_function=np.sum)

    plt.ylim(xlims) # Set X-axis limits
    plt.xlim(ylims) # Set Y-axis limits
    
    # First remove italics from mathtext then label axes
    plt.rc('mathtext',fontset='stixsans',default='regular') # Removes italics mathtext
    plt.xlabel(r'Reflectivity (dBZ)',labelpad=xpad,fontsize=ylabFontSize)
    plt.ylabel(r'log$_{10}$Rainfall Rate (mm h$^{-1}$)',labelpad=ypad,fontsize=xlabFontSize)
    
    ax.set_title(title) # Set the plot title
    
    # Add colorbar
    cb = plt.colorbar(shrink=0.85)
    cb.set_label('Normalized Counts')
    
    del c
    
    return
#**===============================================================
def zr_density_pdf(ax,Z,R,title=None,grid=(80,60),minct=0.01,
                   xlabFontSize=16,xpad=7,ylabFontSize=16,ypad=7,
                   xlims=(-3,3),ylims=(-20.,60.)):
    """Create a density scatterplot using the hex bin method (color marker based on density)
 INPUT::
  ax              = Axis to plot against (produced externally)
  Z               = Reflectivity factor [mm^6/m^3]
  R               = Rainfall rate [mm/h]
     OPTIONAL
  title           = Plot title
  grid            = The bin size in (x,y) to be used
  minct           = Minimum count to be displayed
 OUTPUT::
  p               = Plot
 USAGE::
  p = zr_density_pdf(ax,Z,R,**args)
    """
# HISTORY::
#  10 Mar 2014 - Nick Guy NOAA/NSSL/WRDD, NRC
#--------------------------------------------------------
    # First calculate the 2D PDF
    pdf2_Z_R, xedges, yedges = gl.hist2d_masked(10.*np.log10(Z),np.log10(R),bins=grid,
                                                rng=(ylims,xlims),norm=True)
    
    # Replace any zero values with missing data for nice plots
    pdf2_Z_R = np.ma.masked_equal(pdf2_Z_R,0)
    
    # Flip and rotate the 2D Histogram
    pdf2_Z_R = np.rot90(pdf2_Z_R)
    pdf2_Z_R = np.flipud(pdf2_Z_R)
    
    # Create 2D arrays from xedges, yedges
    X, Y = np.meshgrid(xedges, yedges)
    
    # Plot the data using colormesh
    plt.pcolormesh(X,Y,pdf2_Z_R)#,cmap=

    plt.ylim(xlims) # Set X-axis limits
    plt.xlim(ylims) # Set Y-axis limits
    
    # First remove italics from mathtext then label axes
    plt.rc('mathtext',fontset='stixsans',default='regular') # Removes italics mathtext
    plt.xlabel(r'Reflectivity (dBZ)',labelpad=xpad,fontsize=ylabFontSize)
    plt.ylabel(r'log$_{10}$Rainfall Rate (mm h$^{-1}$)',labelpad=ypad,fontsize=xlabFontSize)
    
    ax.set_title(title) # Set the plot title
    
    # Add colorbar
    cb = plt.colorbar(shrink=0.85)
    cb.set_label('Normalized Counts')
    
    del pdf2_Z_R,xedges,yedges,X,Y
    
    return
#**===============================================================
def plotDate_ts(ax,Time,Var,colF='ko',msize=1.5,lw=2,dForm='%H:%M',tz=None,
                title=None,xlab=None,xlabFontSize=16,xpad=7,ylab=' ',ylabFontSize=16,ypad=7,
                xMinTicker='minute',yMajTicks=None,yMinTicks=None):
    """Create a time series plot, with time on X-axis and variable on Y-axis.
 INPUT::
  ax              = Axis to plot against (produced externally)
  Time            = Time array to plot on x-axis
  Var             = Variable to plot as time series
     OPTIONAL
  colF            = Color and marker shortcut (see python documentation)
  msize           = Marker size
  dForm           = Format of the time string for x-axis labels
  tz              = Time zone info to use when creating axis labels
  title           = Plot title
  xlab          = X-axis label
  ylab          = Y-axis label
  xpad          = Padding for X-axis label
  ypad          = Padding for Y-axis label
  xMinTicker    = Sting to set minor ticks, 'second','minute','hour','day' supported
  yMajTicks     = Values for major tickmark spacing, y-axis
  yMinTicks     = Values for minor tickmark spacing, y-axis
 OUTPUT::
  p               = Plot
 USAGE::
  p = plotDate_ts(ax,Time,Var,**args)
    """
# HISTORY::
#  17 Mar 2014 - Nick Guy NOAA/NSSL/WRDD, NRC
#--------------------------------------------------------
    # Set the date format
    xFmt = DateFormatter(dForm,tz=tz)

    plt.plot_date(Time,Var,fmt=colF,markersize=msize,lw=lw) # Create plot
    
    # Set the x-axis ticks
    ax.xaxis.set_major_formatter(xFmt) # Set the date format
    if xMinTicker == 'second':
        from  matplotlib.dates import SecondLocator
        ax.xaxis.set_minor_locator(SecondLocator()) 
    elif xMinTicker == 'minute':
        from  matplotlib.dates import MinuteLocator
        ax.xaxis.set_minor_locator(MinuteLocator()) 
    elif xMinTicker == 'hour':
        from  matplotlib.dates import HourLocator
        ax.xaxis.set_minor_locator(HourLocator())
    elif xMinTicker == 'day':
        from  matplotlib.dates import DayLocator
        ax.xaxis.set_minor_locator(DayLocator())
        
    # Set the y-axis ticks
    try:
        yMajTicks
        ax.yaxis.set_major_locator(mtic.MultipleLocator(yMajTicks)) # Set the major tickmarks
    except:
        pass
    try:
        yMinTicks
        ax.yaxis.set_minor_locator(mtic.MultipleLocator(yMinTicks)) # Set the minor tickmarks
    except:
        pass
#    ax.yaxis.set_major_locator(mtic.MultipleLocator(yMajTicks)) # Set the y-axis major ticks
#    ax.yaxis.set_minor_locator(mtic.MultipleLocator(yMinTicks)) # set the y-axis minor ticks
    
    plt.tick_params(which='both',direction='out') # Turn the tick marks outward
    
    # Set the Y label
    plt.ylabel(ylab,labelpad=ypad,fontsize=ylabFontSize)
    if xlab == None:
        pass
    else:
        plt.xlabel(xlab,labelpad=xpad,fontsize=xlabFontSize)
#    except:
#        pass
        
    # Set the title
#    try:
#        title
#        ax.set_title(title)
#    except:
#        pass
    if title == None:
        pass
    else:
        ax.set_title(title)
    
    return
#**===============================================================
def plotHov(ax,X,Time,Var,dForm='%H:%M',tz=None,title=None,clevs=7,vmin=-2,vmax=4,
                ylab=None,ylabFontSize=16,ypad=7,xlab=' ',xlabFontSize=16,xpad=7,
                yMajTicker='minute',xMajTicks=None,xMinTicks=None,xmax=6.5,
                cbar=True,cbor='vertical',cblab=' ',cmap='jet'):
    """Create a Hovmoeller plot with time on Y-axis.
 INPUT::
  ax              = Axis to plot against (produced externally)
  X               = Variable array to plot along x-axis
  Time            = Time array to plot on y-axis
  Var             = Variable to contour as Hovmoeller
     OPTIONAL
  dForm           = Format of the time string for x-axis labels
  tz              = Time zone info to use when creating axis labels
  title           = Plot title
  clevs           = Number of contour levels
  vmin            = Minimum contour value to display
  vmax            = Maximum contour value to display
  xlab            = X-axis label
  xlabFontSize    = X-axis Font size
  ylabFontSize    = Y-axis Font size
  ylab            = Y-axis label
  xpad            = Padding for X-axis label
  ypad            = Padding for Y-axis label
  xMinTicker      = Sting to set minor ticks, 'second','minute','hour','day' supported
  yMajTicks       = Values for major tickmark spacing, y-axis
  yMinTicks       = Values for minor tickmark spacing, y-axis
  cbar            = Set True for colorbar
  cbor            = Orientation of colorbar (vertical or horizontal)
  cblab           = Colorbar label
  cmap            = Colortable to use for colormap
 OUTPUT::
  p               = Plot
 USAGE::
  p = plotHov(ax,X,Time,Var,**args)
    """
# HISTORY::
#  17 Mar 2014 - Nick Guy NOAA/NSSL/WRDD, NRC
#--------------------------------------------------------
    # Set the date format
    yFmt = DateFormatter(dForm,tz=tz)
    
    # Create contour level array
    clevels = np.logspace(vmin,vmax,clevs)

    cs = plt.contourf(X,Time,Var,clevels,norm=LogNorm(),cmap=cmap) # Create plot
#    cs = plt.contourf(X,Time,Var,vmin=10.**vmin,vmax=10.**vmax,norm=LogNorm(),cmap=cmap) # Create plot
    plt.xscale('log') # Make X-axis logarithmic
    ax.set_xlim(0.,xmax)
    
    # First remove italics from mathtext then label axes
    plt.rc('mathtext',fontset='stixsans',default='regular') # Removes italics mathtext
      
    # Set the y-axis ticks
    ax.yaxis.set_major_formatter(yFmt) # Set the date format
    if yMajTicker == 'second':
        from  matplotlib.dates import SecondLocator
        ax.yaxis.set_major_locator(SecondLocator()) 
    elif yMajTicker == 'minute':
        from  matplotlib.dates import MinuteLocator
        ax.yaxis.set_major_locator(MinuteLocator()) 
    elif yMajTicker == 'hour':
        from  matplotlib.dates import HourLocator,MinuteLocator
        ax.yaxis.set_major_locator(HourLocator())
        ax.yaxis.set_minor_locator(MinuteLocator(interval=10)) 
    elif yMajTicker == 'day':
        from  matplotlib.dates import DayLocator
        ax.yaxis.set_major_locator(DayLocator())
        
    # Set the x-axis ticks
    try:
        xMajTicks
        ax.xaxis.set_major_locator(mtic.MultipleLocator(xMajTicks)) # Set the major tickmarks
        ax.xaxis.set_major_formatter(mtic.LogFormatter(labelOnlyBase=False))
    except:
        pass
    try:
        xMinTicks
        ax.xaxis.set_minor_locator(mtic.MultipleLocator(xMinTicks)) # Set the minor tickmarks
    except:
        pass
    
    # Set the X and Y labels
    plt.xlabel(xlab,labelpad=xpad,fontsize=xlabFontSize)
    if ylab == None:
        pass
    else:
        plt.ylabel(ylab,labelpad=ypad,fontsize=ylabFontSize)
        
    # Set the title
    if title == None:
        pass
    else:
        ax.set_title(title)
        
    # Add Colorbar   
    l_f = mtic.LogFormatter(10, labelOnlyBase=False) 
    if cbor == 'vertical':
        frac = 0.05
        shrink = 0.6
    cb = plt.colorbar(cs,orientation=cbor,fraction=frac,pad=.05)#,ticks=clevels)#format="%2G",
    cb.set_label(cblab) # Set the colorbar label
    # Set the number of ticks in the colorbar based upon number of contours
#    tick_locator = mtic.MaxNLocator(nbins=int(clevs/2))
#    cb.locator = tick_locator
#    cb.update_ticks()
    
    return
#**===============================================================