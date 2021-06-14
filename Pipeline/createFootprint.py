################################################################################
### Authors: Tiago Castro, Pierluigi Monaco                                  ###
### This is an example script to create a footprint                          ###
###                                                                          ###
################################################################################
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astropy.io import fits
import healpy as hp
from healpy.rotator import Rotator
import sys
from os import path

outdir = '../Products/'

# these are the parameters of the footprint, please copy them to the input file

footprint_fname = outdir + 'Footprints/100sqdeg.fits'
footprint_res = 2048
footprint_tag = "100sqdeg"
footprint_zrange = [0.8, 2.0]

# linear size of the field in radians
size = 10. * np.pi/180.

# reddening map for checking where the footprint has been placed
#reddening = hp.ud_grade(hp.read_map(outdir+'ExtinctionMaps/HFI_CompMap_ThermalDustModel_2048_R1.20.fits', field=2), footprint_res)
reddening = hp.ud_grade(hp.read_map(outdir+'ExtinctionMaps/COM_CompMap_ThermalDust-commander_2048_R2.00.fits', field=2), footprint_res)


# creating boolean healpy map
(theta,phi)=hp.pix2ang(footprint_res,np.arange(hp.nside2npix(footprint_res)))
footprint=np.zeros(hp.nside2npix(footprint_res),dtype=bool)

# square of 10 deg of size on the equator
footprint[(theta>=np.pi/2.-size/2.) & (theta<np.pi/2.+size/2.) & (phi<size)]=True

# rotation to a convenient place in galactic coordinates
rot=Rotator(rot=[-40.0, -40.0, 0.0])
footprint=rot.rotate_map_pixel(footprint)>0.5

sky_fraction = footprint.sum()/footprint.size

# plot of the location
foot2 = reddening.copy()
foot2[footprint] *= 2
hp.mollview(foot2,max=1000)
plt.savefig(outdir+'Plots/100sqdeg.png')

# writes footprint on fits file
print("## writing footprint on file {}".format(footprint_fname))

fg = fits.Column( name='FOOTPRINT_G', array=footprint, format='L')
tb = fits.BinTableHDU.from_columns([fg])
tb.header.append('RES')
tb.header['RES']=footprint_res
tb.header.append('MINZ')
tb.header['MINZ']=footprint_zrange[0]
tb.header.append('MAXZ')
tb.header['MAXZ']=footprint_zrange[1]
tb.header.append('TAG')
tb.header['TAG']=footprint_tag
tb.header.append('SKYFRAC')
tb.header['SKYFRAC']=sky_fraction

tb.writeto(footprint_fname, overwrite=True)

