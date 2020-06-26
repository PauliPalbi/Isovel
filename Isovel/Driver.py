import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
import astropy.constants as const

import fit_shape
from isovel import Isovel 


###########################
#The way of using this and to do a plot
#the values were obtain with fit_shape
###########################

path = '/Volumes/Transcend/J1615_2020/DataWork/Results/CO/'
name_12CO = path+'J1615_12CO_im_uv006'
fits_image = name_12CO +'.fits'

channels= np.arange(-7.99991, 16.1501, 0.35)

image_1 = fits.open(fits_image)
data_1 = image_1[0].data[0][0]
imsize_1 = len(data_1[0])
Rmax = imsize_1/2. # pix
header = image_1[0].header
pixscales=np.abs(3600*header['CDELT2'])
arcs = Rmax * pixscales
arcs = [arcs,-arcs,-arcs,arcs]


mstar, pa, inc, z0, psi, vlsr =(1.2047920150052118, 325.15077530944075, 46.41703391829195, 0.2494024079038776, 1.211576566441788, 4.741979863714648)
pa = 90-pa-180
border = 550.
par=6.341708607533389 #mas
d = 1000/par #pc
pix_to_au = (pixscales*d)

Isovel_disk = Isovel(mstar, pa, inc, z0, psi, vlsr, pix_to_au, border = border)
vel_thin, vel_front, vel_back  = Isovel_disk.velocities()

print(vel_thin, vel_front, vel_back)

fig, axs = plt.subplots(3, 3, sharey=True, figsize=(16, 16))
for i in range(32,41):
    b=i-32
    a=0
    if b>2:
        a=1
        b=i-32-3
    if b>2:
        a=2
        b=i-32-6
    
    image_1 = fits.open(fits_image)
    data_1 = image_1[0].data[0][i]

    image_12mm = np.zeros([imsize_1, imsize_1])
    image_12mm[:][:] = data_1[:][:]
    
    vels_lev = (np.array([channels[i]])-vlsr)
    axs[a,b].imshow(image_12mm,origin='lower', cmap='gist_heat', extent=arcs)
    axs[a,b].contour(vel_front,vels_lev,origin='lower',colors='cyan',linestyles='dashed',linewidths=2.1, extent=arcs, corner_mask=True)
    axs[a,b].contour(vel_back,vels_lev,origin='lower',colors='cyan',linestyles='dotted',linewidths=1.2, extent=arcs, corner_mask=True)

fig.savefig('J1615_isovelocities.png',dpi=600)
