"""pyparticleprobe: Python cloud and precipitation particle probe data analysis.

pyparticleprobe is a toolkit that contains utilities that can be used
 to analyze processed (from raw) microphysical particle probe data.  


"""

from numpy.distutils.core import setup, Extension
import os

#- Pull the header into a variable 
doclines = __doc__.split("\n")

#- Set variables for setup
packages = ['pyparticleprobe']
package_dirs={'pyparticleprobe'}
datafiles = glob.glob(os.path.join(pathout,'*'))
datafiles = [os.path.join('data',os.path.basename(f)) for f in datafiles]
package_data = {'pyparticleprobe':datafiles}

#- Run setup
setup (name = 'pyparticleprobe',
       version = '0.1.0',
       author = 'Nick Guy',
       author_email = 'nick.guy@noaa.gov'
       packages = packages,
       package_dir = package_dirs,
       package_data = package_data,
       url = 'https://github.com/nguy/pyparticleprobe',
       license='LICENSE.txt',
       description = doclines[0],
       long_description = """A toolkit that contains utilities that can be used
 to analyze processed (from raw) microphysical particle probe data.
""",
       install_requires = ['Numpy >=1.7.2',
                           'SciPy >=0.13.3',
                           'pytmatrix>=0.2.0',
                           'os'],
       )
