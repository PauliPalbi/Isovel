#in case you dont know the shape of the gas of your disk, use this
#this calls eddy.fit_cube
# https://github.com/richteague/eddy/tree/master/docs/tutorials

import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from eddy.fit_cube import rotationmap
from multiprocessing import Pool

# use bettermoments to create v0.fits and dv0.fits
# on the terminal do the following
# pip install bettermoments
# bettermoments path/to/cube.fits

def input_variables(path, name_file, 
            d, x0, y0, mstar, vlsr, z0, psi, PA, inc, 
            r_min=1.5, r_max=2.3, 
            downsample=20, clip=3, 
            nwalkers=50, nburnin=500, nsteps=3000, 
            beam=False, just_results=True):
    """
    # if you dont know distance and know par: d = 1000/par #pc  par = xx #mas
    # p0 = [0, 0, 1.248, 4753., 0.27, 1.22,  180+Pa, Inc]

    """
    p0 = [x0, y0, mstar, vlsr, z0, psi, PA, inc] # 8 dim list in this order

    cube, samples, percentiles, dicti = fit_shape(path, name_file, 
                                                d, p0, r_min, r_max, 
                                                downsample=downsample, clip=clip, 
                                                nwalkers=nwalkers, nburnin=nburnin, nsteps=nsteps, 
                                                beam=beam)
    
    if just_results:
        return dicti
    else:
        return cube, samples, percentiles, dicti


# read file and make mcmc to fit parameters to get the shape
def fit_shape(path, name_file, d, p0, r_min, r_max, downsample=20,
            clip=3, nwalkers=50, nburnin=500, nsteps=3000, beam=False):
    """

    """
    cube = rotationmap(path=path+name_file+'v0.fits',
                    uncertainty=path+name_file+'dv0.fits',
                    downsample=downsample,
                    clip=clip)

    # you can fix inc and PA if you know them or let them as free parameters
    # here they are free parameters
    # read eddy documentation
    # for this purpose we will use a simple flared surface as z = z0*r**psi

    params = {}

    # Dictate where the mask is.
    r_min = r_min * cube.bmaj
    r_max = r_max
    returns= ['samples', 'percentiles', 'dict']

    # Start with the positions of the free variables in p0.
    # Note that as we have non-zero z values we must keep Mstar
    # a free parameter to account for this change.

    params['x0'] = 0      #
    params['y0'] = 1      #
    params['mstar'] = 2   #
    params['vlsr'] = 3    #
    params['z0'] = 4      #
    params['psi'] = 5     #
    params['PA'] = 6      #
    params['inc'] = 7     # degrees

    # Fix the other parameters.
    params['dist'] = d  # parsec
    params['beam'] = beam # should we convolve the model? **MUST BE A BOOLEAN**

    # Run the MCMC.
    with Pool() as pool:
        samples, percentiles, dicti = cube.fit_map(p0=p0, 
                                                params=params, 
                                                r_min=r_min, 
                                                r_max=r_max, 
                                                nwalkers=nwalkers, 
                                                nburnin=nburnin, 
                                                nsteps=nsteps, 
                                                pool=pool, 
                                                optimize=False,
                                                returns=returns)
    # Return cube, samples, percentiles and results of best fit in dicti format
    return cube, samples, percentiles, dicti