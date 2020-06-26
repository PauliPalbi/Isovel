# Class that calculate the velocity of the gas 
import numpy as np
import astropy.constants as const
import astropy.units as u
from scipy.ndimage.interpolation import rotate
from scipy import optimize
import matplotlib.pyplot as plt
from astropy.io import fits

class Isovel:
    """
    Import the parameters of a disk 
    Return the velocities of the mid plane and of the surface of the gas. (Front and back side)

    Args:
        mstar (float): Msun units
        pa (float): in degrees
        inc (float): in degrees
        z0 (float): in arcsec
        psi (float): 
        vlsr (float): in km/s
        pix_to_au (float): change of units from pixels to au
        R_max (float): Max radius in au
        v_max (float):
        nx (int): number of points in x
        ny (int): number of points in y
        border (float): Border of the velocity fit in au 
    """
    def __init__(self, mstar, pa, inc, z0, psi, vlsr, pix_to_au,
                R_max=300, v_max=5, nx=300, ny=300, border=550):

        M_sun = const.M_sun.to(u.g)

        self.mstar = mstar * M_sun.value

        self.inc= np.radians(inc) # south is closer

        self.pa = pa
        #pa = 90-85.4526 #check later for correct formula

        z = z0*1**psi
        self.z0 = z0
        self.psi = psi
        self.phi = np.arctan(z)

        self.vlsr=vlsr #km/s
        self.R_max = R_max # au
        self.v_max = v_max
        
        self.nx = nx
        self.ny = ny
        self.pix_to_au = pix_to_au 
        self.border = border


    def velocities(self):
        """
        This perform the fit using eddy.fit_cube

        Args:
            #No inputs directly

        Returns:
            v_near (array): array with keplerian vel in front side 
            v_far (array): array with keplerian vel in back side 
            v_thin (array): array with keplerian velocity in the mid plane
        """
            # projected coordinates (should be xprime, yprime)
        G = const.G.value*1000
        au = 1.496e13

        xp = np.linspace(-self.R_max, self.R_max, self.nx)
        yp = np.linspace(-self.R_max, self.R_max, self.ny)
        xxp, yyp = np.meshgrid(xp, yp)

        xx = xxp 
        yy = yyp/np.cos(self.inc)
        rr_cm = np.sqrt(xx**2+yy**2) * self.pix_to_au *au # cm
        theta = np.arctan2(yy,xx)

            # projected line of sight velocity at x',y'
        vel_thin = np.sqrt(G*self.mstar/rr_cm)*np.sin(self.inc)*np.cos(theta) # cm/s
        vel_thin /= 1e5 # km/s
        
            # rotate
        vel_thin = rotate(vel_thin,-self.pa, reshape=False)

        seci = 1/np.cos(self.inc)
        def f(t,ixp,iyp):
            return (np.cos(2*self.inc) + np.cos(2*self.phi))*t**2 - (2*np.sin(self.phi)**2*2*iyp*np.tan(self.inc))*t - (2*np.sin(self.phi)**2*(ixp**2 + iyp**2*seci**2 ))

        tt_near = np.zeros(yyp.shape) # near side of the disk
        tt_far = np.zeros(yyp.shape) # far side of the disk

        for ix in range(self.nx):
            for iy in range(self.ny):
                sol = optimize.root(f, [0.], args=(xxp[ix,iy],yyp[ix,iy]), method='lm')
                if np.degrees(self.inc) < 180:
                    tt_near[ix,iy] = np.abs(sol.x) # au
                    tt_far[ix,iy] = -np.abs(sol.x)
                else:
                    tt_far[ix,iy] = np.abs(sol.x) # au
                    tt_near[ix,iy] = -np.abs(sol.x)

        xx = xxp
        yy_near = yyp/np.cos(self.inc) + tt_near*np.sin(self.inc)
        zz_near = tt_near*np.cos(self.inc)
        rr_near_cm = np.sqrt(xx**2+yy_near**2+zz_near**2) * self.pix_to_au *au # cm
        theta_near = np.arctan2(yy_near,xx)

        yy_far = yyp/np.cos(self.inc) + tt_far*np.sin(self.inc)
        zz_far = tt_far*np.cos(self.inc)
        rr_far_cm = np.sqrt(xx**2+yy_far**2+zz_far**2) * self.pix_to_au*au # cm
        theta_far = np.arctan2(yy_far,xx)

        vel_near = np.sqrt(G* self.mstar/rr_near_cm)*np.sin(self.inc)*np.cos(theta_near)/1e5 # km/s
        vel_far = np.sqrt(G* self.mstar/rr_far_cm)*np.sin(self.inc)*np.cos(theta_far)/1e5 # km/s

        # rotate
        vel_near = rotate(vel_near,self.pa, reshape=False)
        vel_far = rotate(vel_far,self.pa, reshape=False)

        '''
        mask_contour_near = rotate(rr_near_cm/au, self.pa, reshape=False)[:][:]>self.border
        mask_contour_far= rotate(rr_far_cm/au, self.pa, reshape=False)[:][:]>self.border
        mask_contour_thin = rotate(rr_cm/au, self.pa, reshape=False)[:][:]>self.border
        
        vel_near = np.ma.array(vel_near, mask=mask_contour_near)
        vel_far = np.ma.array(vel_far, mask=mask_contour_far)
        vel_thin = np.ma.array(vel_thin, mask=mask_contour_thin)
        '''
        return vel_near, vel_far, vel_thin
