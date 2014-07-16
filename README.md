pyparticleprobe
===============

A python package to process cloud and precipitation particle probe  data.

Author: Nick Guy (nick.guy@noaa.gov)

Created:  21 Jan 2014
Updated:  14 May 2014 – Restructured functions and pathways
Updated:  15 Jul 2014 - Refactored to take class structure.  Heavily influenced by the PyDisdrometer package by Joseph Hardin.

## Installation
As of now, I have used these as a standalone package that was not installed as a package.
To access, I just added the path where the folder was unpacked to my 
PYTHONPATH environmental variable path (in my .bashrc file)

e.g. export PYTHONPATH=/Users/nickguy/programs/python/pythonlib

There is a setup.py file provided, HOWEVER this may not work properly and I would suggest
 using the above method to access the package until the software becomes more mature.
 
## Usage
Initiate class by opening a file.  There are a number of optional arguments that help in processing the file.
As of now these attributes are tailored to files created after processing the CIP/PIP probes aboard the NOAA P-3 aircraft.

```python
dsd = pyparticleprobe(filename)
```

To add additional analyis information:
```python
dsd.get_bin_parameters() # Bin the data over diameters bin

dsd.get_ts_parameters() # Bin the data over each time step

dsd.get_model_parameters(modelName='gamma_ua98') # Get the model fit parameters, default uses gamma model as in Ulbrich and Atlas (1998)
```

## Requirements
This library borrow heavily from the [PyDisdrometer](https://github.com/josephhardinee/PyDisdrometer) package.
Therefore it uses the typical scientific python stack

It also requires the [PyTMatrix Package](https://github.com/jleinonen/pytmatrix).

## Notes
This code is in beta format and works with data from the NOAA P-3 aircraft that has been processed.  It should be extendable to other probe data.  

Future enhancements will hopefully include the ability to read other data files.  
Also, the plotting portion depends on a local library on my machine and therefore won't work.  This will eventually be changed to allow some plots to be produced.


