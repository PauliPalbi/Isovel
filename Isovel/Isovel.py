# This to start documentation
import numpy as np
import astropy.constants as const
from scipy.ndimage.interpolation import rotate
from scipy import optimize
#'PA': 325.03970248937446,
# 'self.mstar': 1.1835185573440012,
# 'inc': 46.88375837648006,
# 
# 'dist': 157.68621074959017,
# 'vlsr': 4744.619070925755,
#
# 'z0': 0.24717562173400856,
# 'psi': 1.1814702194561688

class Isovel:
    def __init__(self, mstar, pa, inc, z0, psi, vlsr, vels_lev, 
                R_max=300, v_max=5, nx=300, ny=300, southcloser=True, pix_to_au =1):

        M_sun = const.M_sun

        self.mstar = mstar * M_sun.value

        if southcloser:
            self.inc= np.radians(180+inc) # south is closer
        else:
            self.inc = np.radians(360-inc) # north is closer

        self.pa = pa
        #pa = 90-85.4526 #check later for correct formula

        z = z0*1**psi
        self.z0 = z0
        self.psi = psi
        self.phi = np.arctan(z)

        self.vlsr=vlsr #km/s
        self.R_max = R_max # au
        self.v_max = v_max
        self.vels_lev = vels_lev - self.vlsr
        #vels_lev = np.array([3,3.7,4.4,5.1,5.8,6.5]) - vlsr

        self.nx = nx
        self.ny = ny
        self.pix_to_au = pix_to_au 


    def Keplerian_rotation(self, rr_cm, theta):
        """
        multiplication

        Args: 
            
        Returns:
        
            float: multiplication
        """
        G = const.G
        vel = np.sqrt(G.value * self.mstar / rr_cm) 
        vel *= np.cos(theta) * np.sin(np.radians(self.inc)) # km/s
        
        # rotate
        pa = np.radians(self.pa)
        vel = rotate(vel, pa)
        return vel

    def create_shape(self):
        """
        """
        xp = np.linspace(-self.R_max, self.R_max, self.nx)
        yp = np.linspace(-self.R_max, self.R_max, self.ny)
        xxp, yyp = np.meshgrid(xp, yp)
        return xxp, yyp


    def deproject_coordinates(self, x, y):
        """
        multiplication

        Args: 
            
        Returns:
        
            float: multiplication
        """
        inc = np.radians(self.inc)
        xx = x
        yy = y/np.cos(inc)
        return xx, yy


    def cilindrical_coords(self, xx, yy, zz, Dim_3 = False):
        if not Dim_3:
            rr = np.sqrt(xx**2 + yy**2)*self.pix_to_au # au
        else:
            rr = np.sqrt(xx**2 + yy**2 + zz**2)*self.pix_to_au  # au
        au = 1.496e8 # from au to km
        rr_km = rr*au # km
        theta = np.arctan2(yy,xx)
        return rr_km, theta


    def f(self, t, ixp, iyp):
        """
        we assume flared is a cone

        intersection of line with cone
        # Flared disk -> cone

        # intersection of line with cone:
        # yp = (y - t*sin(i))*cosi
        # zp = t*cos(i)

        # where t are the solution of a line intersecting a cone
        # t**2 * (np.cos(2*i) + np.cos(2*phi)) - 2*np.sin(phi)**2 * (xxp**2 + yyp**2*np.sec(i)**2 + 2*t*yyp*np.tan(i)) = 0
        # The positive and negative roots of this equation correspond tothe front and back halves
        # of the double cone, respectively

        # solve for t

        Args: 
            t (array):are the solution of a line intersecting a cone
            ixp ():
            iyp ():
        Returns:
            array:

        """
        inc = np.radians(self.inc)
        seci = 1/np.cos(inc)

        a = (np.cos( 2 * inc) + np.cos( 2 * self.phi)) * t**2 
        b = - (2 * np.sin(self.phi)**2 * 2 * iyp * np.tan(inc)) * t 
        c = - (2 * np.sin(self.phi)**2 * (ixp**2 + iyp**2 * seci**2 ))
        return a+b+c

    def surface_projection(self, xxp, yyp):
        """
        """
        #vel_front
        tt_front = np.zeros(yyp.shape) # front side of the disk
        #vel_back
        tt_back = np.zeros(yyp.shape) # back side of the disk
        for ix in range(self.nx):
            for iy in range(self.ny):
                sol = optimize.root(f, [0.], args=(xxp[ix,iy],yyp[ix,iy]), method='lm')
                if self.inc<180:
                    tt_front[ix,iy] = np.abs(sol.x) # au
                    tt_back[ix,iy] = -np.abs(sol.x)
                else:
                    tt_back[ix,iy] = np.abs(sol.x) # au
                    tt_front[ix,iy] = -np.abs(sol.x)
        return tt_front, tt_back


    def calculate_vel(self, border=550):
        """
        """
        
        xxp, yyp = self.create_shape()
        x_dep, y_dep = self.deproject_coordinates(xxp, yyp)

        z=np.array([])
        rr, theta = self.cilindrical_coords(x_dep, y_dep, z) #km and rad

        vel_thin = self.Keplerian_rotation(rr, theta)

        tt_front, tt_back = self.surface_projection(x_dep, y_dep)

        inc = np.radians(self.inc)

        yy_front = y_dep + tt_front*np.sin(inc)
        zz_front = tt_front*np.cos(inc)
        print(zz_front)
        rr_front, theta_front = self.cilindrical_coords(x_dep, yy_front, zz_front, Dim_3=True)
        
        yy_back = y_dep + tt_back*np.sin(inc)
        zz_back = tt_back*np.cos(inc)
        rr_back, theta_back = self.cilindrical_coords(x_dep, yy_back, zz=zz_back, Dim_3=True)
        

        vel_front = self.Keplerian_rotation(rr_front, theta_front)
        vel_back = self.Keplerian_rotation(rr_back, theta_back)

        pa = np.radians(self.pa)
        mask_contour_front = rotate(rr_front, pa, reshape=False)[:][:]>border
        mask_contour_back= rotate(rr_back, pa, reshape=False)[:][:]>border
        mask_contour_thin = rotate(rr, pa, reshape=False)[:][:]>border

        vel_front = np.ma.array(vel_front, mask=mask_contour_front)
        vel_back = np.ma.array(vel_back, mask=mask_contour_back)
        vel_thin = np.ma.array(vel_thin, mask=mask_contour_thin)

        return vel_thin, vel_front, vel_back




    def power_z_function(self, r):
        """
        power law of surface shape, not used for now


        """
        z = self.z0 * np.power(r, self.psi)
        return z



    '''
    def plot_velocities(self):
            # Plot
        fig, (ax1,ax2,ax3,ax4) = pl.subplots(1,4,figsize=(10,10))
        ext=[-Rmax,Rmax,-Rmax,Rmax]

        im1 = ax1.imshow(vel_thin,extent=ext,origin='lower',vmin=-vmax,vmax=vmax)
        ax1.contour(vel_thin,vels_lev,extent=ext,origin='lower',colors='white',linestyles='solid',linewidths=2)

        im2 = ax2.imshow(vel_front,extent=ext,origin='lower',vmin=-vmax,vmax=vmax)
        ax2.contour(vel_front,vels_lev,extent=ext,origin='lower',colors='white',linestyles='dashed',linewidths=2)

        im3 = ax3.imshow(vel_back,extent=ext,origin='lower',vmin=-vmax,vmax=vmax)
        ax3.contour(vel_back,vels_lev,extent=ext,origin='lower',colors='white',linestyles='dotted',linewidths=2)

        cols = ['red','red','red','blue','blue','blue']
        im4 = ax4.imshow(vel_thin*0,extent=ext,origin='lower',cmap='gray_r')
        #ax4.contour(vel_thin,[-1.05,1.05],extent=ext,origin='lower',colors=cols,linestyles='solid',linewidths=2)
        #ax4.contour(vel_front,[-1.05,1.05],extent=ext,origin='lower',colors=cols,linestyles='dashed',linewidths=2)
        #ax4.contour(vel_back,[-1.05,1.05],extent=ext,origin='lower',colors=cols,linestyles='dotted',linewidths=2)

        ax4.contour(vel_front,vels_lev,extent=ext,origin='lower',colors=cols,linestyles='dashed',linewidths=2)
        ax4.contour(vel_back,vels_lev,extent=ext,origin='lower',colors=cols,linestyles='dotted',linewidths=2)

        for ax in ax1,ax2,ax3,ax4:
            ax.set_xlim([-Rmax,Rmax])
            ax.set_ylim([-Rmax,Rmax])


        pl.show()
    '''


mstar, pa, inc, z0, psi, vlsr =(1.2047920150052118, 325.15077530944075, 46.41703391829195, 0.2494024079038776, 1.211576566441788, 4.741979863714648)
vels_lev = np.array([3,3.7,4.4,5.1,5.8,6.5])

Isovel = Isovel(mstar, pa, inc, z0, psi, vlsr, vels_lev)
vel_thin, vel_front, vel_back = Isovel.calculate_vel()