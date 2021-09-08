import numpy as np
import matplotlib.pyplot as plt
import healpy as hp
from healpy.rotator import Rotator
from astropy.io  import fits

rotE2G = Rotator(coord=['C','G'])

map_E_fname = 'sc8_varmap_E.fits'
map_G_fname = 'sc8_varmap_G.fits'

print("reading {}...".format(map_E_fname))
map_E = hp.read_map(map_E_fname,partial=True)
print("rotating...")
map_G = rotE2G.rotate_map_pixel(map_E)
print("fixing unseen...")
map_G[map_G<0]=hp.UNSEEN
print("writing {}...".format(map_G_fname))
hp.write_map(map_G_fname, map_G, partial=True)

map_E_fname = 'sc8_expmap_E.fits'
map_G_fname = 'sc8_expmap_G.fits'

print("reading {}...".format(map_E_fname))
map_E = hp.read_map(map_E_fname,partial=True)
print("rotating...")
map_G = rotE2G.rotate_map_pixel(map_E)
print("fixing unseen...")
map_G[map_G<0]=hp.UNSEEN
print("writing {}...".format(map_G_fname))
hp.write_map(map_G_fname, map_G, partial=True)

print("### writing footprint on file ../Footprints/BenSC8.fits")

footprint_G = map_G != hp.UNSEEN
sky_fraction = footprint_G.sum()/footprint_G.size

fg = fits.Column( name='FOOTPRINT_G', array=footprint_G, format='L')
tb = fits.BinTableHDU.from_columns([fg])
tb.header.append('RES')
tb.header['RES']=4096
tb.header.append('MINZ')
tb.header['MINZ']=0.85
tb.header.append('MAXZ')
tb.header['MAXZ']=1.85
tb.header.append('TAG')
tb.header['TAG']='BenSC8'
tb.header.append('SKYFRAC')
tb.header['SKYFRAC']=sky_fraction

tb.writeto('../Footprints/BenSC8.fits', overwrite=True)

print("done")
