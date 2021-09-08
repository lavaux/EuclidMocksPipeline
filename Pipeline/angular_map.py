################################################################################
### Authors: Tiago Castro, Pierluigi Monaco                                  ###
###                                                                          ###
################################################################################
import numpy as np
import matplotlib.pyplot as plt
import healpy as hp
from astropy.io import fits
import sys

print("# Running angular_map.py")

if len(sys.argv)<5:
    print("Usage: pyton {} [galaxy catalog] [selection] [z1] [z2]".format(sys.argv[0]))
    sys.exit(0)

Nside=256

fname = sys.argv[1]
selection = sys.argv[2]
if selection == 'None':
    selection = None

z1 = np.float(sys.argv[3])
z2 = np.float(sys.argv[4])

print("Plotting map of catalog {} in z=[{},{}]".format(fname,z1,z2))

print("Reading catalog...")
cat = fits.getdata(fname)

if selection is not None:
    print("Reading selection...")
    mysel = fits.getdata(selection)['SELECTION']
else:
    mysel = np.ones(len(cat),dtype=bool)

print("Filtering catalog...")
filter = (cat['true_redshift_gal']>=z1) & (cat['true_redshift_gal']<z2) & mysel

print("Finding pixels for {} galaxies...".format(filter.sum()))
conv = np.pi/180.
pix  = hp.ang2pix(Nside, np.pi/2. - cat['dec_gal'][filter]*conv, cat['ra_gal'][filter]*conv)

print("Constructing map...")
map = np.zeros((hp.nside2npix(Nside)), dtype=int)
for p in pix:
    map[p]+=1

hp.mollview(map,rot=[0,0,0],title=fname)

plt.show()

print("Done!")

