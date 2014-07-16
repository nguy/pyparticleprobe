pyparticleprobe
===============

A python package to process cloud and precipitation particle probe  data.

Author: Nick Guy nick.guy@noaa.gov

History:
Created:  21 Jan 2014	Nick Guy (NRC; NOAA/NSSL)
Updated:  14 May 2014 – Restructured functions and pathways
Updated:  15 Jul 2014 - Refactored to class structure, heavily influenced by the 
                        PyDisdrometer package by Joseph Hardin.

## Installation
As of now, I have used these as a standalone package that was not installed as a package.
To access, I just added the path where the folder was unpacked to my 
PYTHONPATH environmental variable path (in my .bashrc file)

e.g. export PYTHONPATH=/Users/nickguy/programs/python/pythonlib

There is a setup.py file provided, HOWEVER this may not work properly and I would suggest
 using the above method to access the package until the software becomes more mature.
 
 
## Usage
```python
import pyparticleprobe as probe

dsd = probe.read_file(filename)
```

Note that there are a number of optional arguments that help to define the reader.  
Currently the code largely depends on a particular format, but this could be modified 
as more examples of data are found.

To add parameters to the dsd class:
```python
dsd.get_bin_parameters()  # Bin the data over diameters bin

dsd.get_ts_parameters() # Bin the data over each time step

dsd.get_model_parameters(modelName='gamma_ua98') # Get the model fit parameters
```

Note that the default for the model fit is the gamma model as in Ulbrich and Atlas (1998).

## Requirements
This library draws heavily from the [PyDisdrometer] (https://github.com/josephhardinee/PyDisdrometer) package
and also requires the 
[PyTMatrix](https://github.com/jleinonen/pytmatrix) package

NOTES::

This is open-source software, with no warranties extended.

Development is ongoing, with features added on an as needed or suggested basis.

The plot folders currently require a library I use locally.  If you would like that 
library just let me know.   Eventually it will be folded into this project.
