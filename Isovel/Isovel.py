# This to start documentation
import numpy as np
import astropy.constants as const
from scipy.ndimage.interpolation import rotate
from scipy import optimize
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
    def __init__(self, mstar, pa, inc, z0, psi, vsou, vels_lev, 
                    R_max=300, v_max=5, southcloser=True):
        M_sun = const.M_sun
        self.mstar = mstar * M_sun.value

        if southcloser:
            self.inc= np.radians(180+inc) # south is closer
        else:
            self.inc = np.radians(360-inc) # north is closer

        self.pa = pa
        #pa = 90-85.4526 #check later for correct formula

        z = z0*1**psi
        self.phi= np.arctan(z)

        self.vsou=vsou
        self.R_max = R_max # au
        self.v_max = v_max
        self.vels_lev = vels_lev
        #vels_lev = np.array([3,3.7,4.4,5.1,5.8,6.5]) - vsou
        print(vels_lev)

        nx, ny = (100,100)


def Keplerian_rotation(mstar, inc, pa, rr_cm, theta):
    """
    multiplication

    Args: 
        
    Returns:
    
        float: multiplication
    """
    G = const.G
    vel = np.sqrt(G.value * mstar / rr_cm) 
    vel *= np.cos(theta) * np.sin(inc) # km/s
    
    # rotate
    vel = rotate(vel, pa)
    return vel

########################################################


def create_shape(R_max, nx=100, ny=100):
    """
    """
    nx, ny = (nx,ny)
    xp = np.linspace(-R_max, R_max, nx)
    yp = np.linspace(-R_max, R_max, ny)
    xxp, yyp = np.meshgrid(xp, yp)
    return xxp, yyp


def deproject_coordinates(x, y, inc):
    """
    multiplication

    Args: 
        
    Returns:
    
        float: multiplication
    """
    inc = np.radians(inc)
    xx = x
    yy = y/np.cos(inc)
    return xx, yy

def cilindrical_coords(xx, yy, zz=0):
    rr = np.sqrt(xx**2 + yy**2 + zz**2) # au
    au = 1.496e8 # from au to km
    rr_cm = rr*au # km
    theta = np.arctan2(yy,xx)
    return rr_cm, theta


def flared_disk(t, ixp, iyp, inc, phi):
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
    seci = 1/np.cos(inc)

    a = (np.cos( 2 * inc) + np.cos( 2 * phi)) * t**2 
    b = - (2 * np.sin(phi)**2 * 2 * iyp * np.tan(inc)) * t 
    c = - (2 * np.sin(phi)**2 * (ixp**2 + iyp**2 * seci**2 ))
    return a+b+c

def surface_projection(xxp, yyp, inc, phi, nx, ny):
    """
    """
    #vel_front
    tt_front = np.zeros(yyp.shape) # front side of the disk
    #vel_back
    tt_back = np.zeros(yyp.shape) # back side of the disk
    for ix in range(nx):
        for iy in range(ny):
            sol = optimize.root(flared_disk, [0.], args=(xxp[ix,iy],yyp[ix,iy], inc, phi), method='lm')
            if np.degrees(inc)<180:
                tt_front[ix,iy] = np.abs(sol.x) # au
                tt_back[ix,iy] = -np.abs(sol.x)
            else:
                tt_back[ix,iy] = np.abs(sol.x) # au
                tt_front[ix,iy] = -np.abs(sol.x)
    return tt_front, tt_back


def do_stuff(R_max, mstar, inc, pa, phi, nx, ny):
    """
    """
    
    xxp, yyp = create_shape(R_max)
    x_dep, y_dep = deproject_coordinates(xxp, yyp, inc)

    rr, theta = cilindrical_coords(x_dep, y_dep) #km and rad

    vel_thin = Keplerian_rotation(mstar, inc, pa, rr, theta)

    tt_front, tt_back = surface_projection(x_dep, y_dep, inc, phi, nx, ny)


    yy_front = y_dep + tt_front*np.sin(inc)
    zz_front = tt_front*np.cos(inc)
    rr_front, theta_front = cilindrical_coords(x_dep, yy_front, zz=zz_front)
    
    yy_back = y_dep + tt_back*np.sin(inc)
    zz_back = tt_back*np.cos(inc)
    rr_back, theta_back = cilindrical_coords(x_dep, yy_back, zz=zz_back)
    

    vel_front = Keplerian_rotation(mstar, inc, pa, rr_front, theta_front)
    vel_back = Keplerian_rotation(mstar, inc, pa, rr_back, theta_back)

    return vel_thin, vel_front, vel_back




def power_z_function(r, z0, psi):
    """

    """
    z = z0 * np.power(r, psi)
    return z



'''
def plot_velocities():
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


