#in case you dont know the shape of the gas of your disk, use this
#this calls eddy.fit_cube
# https://github.com/richteague/eddy/tree/master/docs/tutorials

# use bettermoments to create v0.fits and dv0.fits
# on the terminal do the following
# pip install bettermoments
# bettermoments path/to/cube.fits

import numpy as np
import matplotlib.pyplot as plt
from eddy.fit_cube import rotationmap
from multiprocessing import Pool


def fit_Gas(path, name_file, 
            d, vlsr, z0, psi, PA, inc, mstar=1,
            x0=0, y0=0,
            r_min=1.5, r_max=2.3, 
            downsample=20, clip=3, 
            nwalkers=50, nburnin=500, nsteps=3000, 
            beam=False, just_results=True):
    """
    This perform a fit to the gas data that has to be saved as the results after applying bettermoments
    It can return just the fit best values or more features

    Args:
        path (str)          :The path were you saved the bettermoments results (name_v0.fits and name_dv0.fits)
        name_file (str)     :name of your Cube fits file
        d (float)           :distance in parsecs
        vlsr (float)        :velosity of the system in m/s
        z0 (float)          : amplitude value of power law in arcsec
        psi (float)         : flared value of power law
        PA (float)          : position angle of the source in deg
        inc (float)         : inclination of the disk in deg
        mstar=1 (float)     :(optional) mass of the central star in units of Msun
        x0=0 (float)        :(optional) offset of the cube in horizontal direction in units of arcsec
        y0=0 (float)        :(optional) offset of the cube in vertical direction in units of arcsec
        r_min=1.5 (float)   :(optional) min radius for the mask #can be check by using just_results=False in units of arcsec
        r_max=2.3 (float)   :(optional) max radius for the mask #can be check by using just_results=False in units of arcsec
        downsample=20 (int) :(optional) Downsample value
        clip=3 (int)        :(optional) Clip value
        nwalkers=50 (int)   :(optional) Number of walkers
        nburnin=500 (int)   :(optional) Number of stepts that will be burn
        nsteps=3000 (int)   :(optional) Tot number of steps
        beam=False (bool)   :(optional) Use this as a default
        just_results(bool)  :(optional) To return just best fit. Use =False to return more features

    Returns:
        dicti (dictionary)  : Disctionary with all best fit values and some more important or useful values
        cube (class)        :(optional)  Cube to work with in eddys code
        samples (np.array)  :(optional) array that has the shape (iterations, nwalkers) with all steps of mcmc
        percentiles         :(optional) Return bestfit values and percentiles
    """
    # if you dont know distance and know par: d = 1000/par #pc  par = xx #mas
    # example of p0 for J1615 : p0 = [0, 0, 1.248, 4753., 0.27, 1.22,  180+Pa, Inc]

    p0 = [x0, y0, mstar, vlsr, z0, psi, PA, inc] # 8 dim list in this order

    cube, samples, percentiles, dicti = fit_eddy(path, name_file, 
                                                d, p0, r_min, r_max, 
                                                downsample=downsample, clip=clip, 
                                                nwalkers=nwalkers, nburnin=nburnin, nsteps=nsteps, 
                                                beam=beam)
    
    if just_results:
        return dicti
    else:
        return cube, samples, percentiles, dicti


# read file and make mcmc to fit parameters to get the shape
def fit_eddy(path, name_file, d, p0, r_min, r_max, downsample=20,
            clip=3, nwalkers=50, nburnin=500, nsteps=3000, beam=False):
    """
    This perform the fit using eddy.fit_cube

    Args:
        path (str)              : The path were you saved the bettermoments results (name_v0.fits and name_dv0.fits)
        name_file (str)         : name of your Cube fits file
        d (float)               : distance in parsecs
        p0 (list)               : python 8 dim list with first guest of bestfit
        r_min=1.5 (float)       :(optional) min radius for the mask #can be check by using just_results=False in units of arcsec
        r_max=2.3 (float)       :(optional) max radius for the mask #can be check by using just_results=False in units of arcsec
        downsample=20 (int)     :(optional) Downsample value
        clip=3 (int)            :(optional) Clip value
        nwalkers=50 (int)       :(optional) Number of walkers
        nburnin=500 (int)       :(optional) Number of stepts that will be burn
        nsteps=3000 (int)       :(optional) Tot number of steps
        beam=False (bool)       :(optional) Use this as a default
        just_results(bool)      :(optional) To return just best fit. Use =False to return more features

    Returns:
        dicti (dictionary)      :Disctionary with all best fit values and some more important or useful values
        cube (class)            :Cube to work with in eddys code
        samples (np.array)      :Array that has the shape (iterations, nwalkers) with all steps of mcmc
        percentiles             :Return bestfit values and percentiles
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

    