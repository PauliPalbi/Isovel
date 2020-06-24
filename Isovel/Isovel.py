# This to start documentation
import numpy as np
import astropy.constants as const
#'PA': 325.03970248937446,
# 'mstar': 1.1835185573440012,
# 'inc': 46.88375837648006,
# 
# 'dist': 157.68621074959017,
# 'vlsr': 4744.619070925755,
#
# 'z0': 0.24717562173400856,
# 'psi': 1.1814702194561688

class disk:
    def __init__(self,mstar,PA,inc,r0,z0, psi,vsou,vels_lev,Rmax = 300,vmax = 5,southcloser=True)
        self.mstar = mstar*const.M_sun.value
        if southcloser:
            self.inc= np.radians(180+inc) # south is closer
        else:
            self.inc = np.radians(360-34.88677) # north is closer
        self.PA = PA
        #pa = 90-85.4526 #check later for correct formula
        r0=1
        z=z0*r0**psi
        self.phi= np.arctan(z)
        self.vsou=vsou
        self.Rmax = Rmax # au
        self.vmax = vmax
        self.vels_lev = vels_lev
        #vels_lev = np.array([3,3.7,4.4,5.1,5.8,6.5]) - vsou
        print vels_lev



def Isovel(PA, mstar, inc, dist, vlsr, z0, psi):
    """
    multiplication

    Args: 
        a (float): first number
        b (float): second number
    Returns:
        float: multiplication
    """


    return a*b
