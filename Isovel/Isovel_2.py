# This to start documentation
import numpy as np
import astropy.constants as const
from scipy.ndimage.interpolation import rotate
from scipy import optimize



def Keplerian_rotation(mstar,theta,inc,pa, rr_cm):
    """
    multiplication

    Args: 
        
    Returns:
    
        float: multiplication
    """
    G = const.G
    vel = np.sqrt(G.value * mstar / rr_cm) 
    vel *= np.cos(theta) * np.sin(np.radians(inc)) # km/s
    
    # rotate
    pa = np.radians(pa)
    vel = rotate(vel, pa)
    return vel

def create_shape(R_max, nx, ny):
    """
    """
    xp = np.linspace(-R_max, R_max, nx)
    yp = np.linspace(-R_max, R_max, ny)
    xxp, yyp = np.meshgrid(xp, yp)
    return xxp, yyp


def deproject_coordinates(inc, x, y):
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


def cilindrical_coords(pix_to_au, xx, yy, zz, Dim_3 = False):
    if not Dim_3:
        rr = np.sqrt(xx**2 + yy**2)*pix_to_au # au
    else:
        rr = np.sqrt(xx**2 + yy**2 + zz**2)*pix_to_au  # au
    au = 1.496e8 # from au to km
    rr_km = rr*au # km
    theta = np.arctan2(yy,xx)
    return rr_km, theta


def f(t, ixp, iyp):
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
    inc = np.radians(inc)
    seci = 1/np.cos(inc)

    a = (np.cos( 2 * inc) + np.cos( 2 * phi)) * t**2 
    b = - (2 * np.sin(phi)**2 * 2 * iyp * np.tan(inc)) * t 
    c = - (2 * np.sin(phi)**2 * (ixp**2 + iyp**2 * seci**2 ))
    return a+b+c

def surface_projection(inc, phi, xxp, yyp, nx, ny):
    """
    """
    #vel_front
    tt_front = np.zeros(yyp.shape) # front side of the disk
    #vel_back
    tt_back = np.zeros(yyp.shape) # back side of the disk
    inc= inc
    phi=phi
    for ix in range(nx):
        for iy in range(ny):
            sol = optimize.root(f, [0.], args=(xxp[ix,iy],yyp[ix,iy]), method='lm')
            if inc<180:
                tt_front[ix,iy] = np.abs(sol.x) # au
                tt_back[ix,iy] = -np.abs(sol.x)
            else:
                tt_back[ix,iy] = np.abs(sol.x) # au
                tt_front[ix,iy] = -np.abs(sol.x)
    return tt_front, tt_back


def calculate_vel(mstar, pa, inc, z0, psi, vlsr, vels_lev, R_max=300, v_max=5, nx=300, ny=300, southcloser=True, pix_to_au =1, border=550):
    """
    """
    phi = np.arctan(z0*1**psi)
    
    xxp, yyp = create_shape(R_max, nx, ny)
    x_dep, y_dep = deproject_coordinates(inc, xxp, yyp)

    z=np.array([])
    rr, theta = cilindrical_coords(pix_to_au, x_dep, y_dep, z, Dim_3 = False) #km and rad

    vel_thin = Keplerian_rotation(mstar,theta,inc,pa, rr)

    tt_front, tt_back = surface_projection(inc, phi, xxp, yyp, nx, ny)

    Inc = np.radians(inc)

    yy_front = y_dep + tt_front*np.sin(Inc)
    zz_front = tt_front*np.cos(Inc)
    print(zz_front)
    rr_front, theta_front = cilindrical_coords(pix_to_au, x_dep, yy_front, zz_front, Dim_3=True)
    
    yy_back = y_dep + tt_back*np.sin(Inc)
    zz_back = tt_back*np.cos(Inc)
    rr_back, theta_back = cilindrical_coords(pix_to_au, x_dep, yy_back, zz_back, Dim_3=True)
    

    vel_front = Keplerian_rotation(mstar,theta_front, inc, pa, rr_front)
    vel_back = Keplerian_rotation(mstar,theta_back, inc, pa, rr_back)

    pa = np.radians(pa)
    mask_contour_front = rotate(rr_front, pa, reshape=False)[:][:]>border
    mask_contour_back= rotate(rr_back, pa, reshape=False)[:][:]>border
    mask_contour_thin = rotate(rr, pa, reshape=False)[:][:]>border

    vel_front = np.ma.array(vel_front, mask=mask_contour_front)
    vel_back = np.ma.array(vel_back, mask=mask_contour_back)
    vel_thin = np.ma.array(vel_thin, mask=mask_contour_thin)

    return vel_thin, vel_front, vel_back




def power_z_function(z0, psi, r):
    """
    power law of surface shape, not used for now


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

'''
R_max=300
v_max=5
nx=300
ny=300
southcloser=True
pix_to_au =1

mstar, pa, inc, z0, psi, vlsr =(1.2047920150052118, 325.15077530944075, 46.41703391829195, 0.2494024079038776, 1.211576566441788, 4.741979863714648)
vels_lev = np.array([3,3.7,4.4,5.1,5.8,6.5])

Isovel = Isovel(mstar, pa, inc, z0, psi, vlsr, vels_lev)
vel_thin, vel_front, vel_back = Isovel.calculate_vel()
'''