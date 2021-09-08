################################################################################
### Author: Tiago Castro, Pierluigi Monaco                                   ###
###                                                                          ###
################################################################################
import numpy as np
import zarr
from astropy.io import fits
import healpy as hp
import sys
from os import path
from euclid_obssys import readConfig

#
# TO DO: INSERT ROTATION AND REPORT IT IN THE FOOTPRINT
#

if len(sys.argv)<2:
    print("Usage: python {} [my input file]".format(sys.argv[0]))
    sys.exit(0)
try: 
    input = readConfig(sys.argv[1])
except Exception as e:
    print(e)
    print("input file not found")
    print("Usage: python {} [my input file]".format(sys.argv[0]))
    sys.exit(0)

print("# Running applyFootprintToMasterCatalog.py with {}".format(sys.argv[1]))
print("# loading catalog...")

# input raw catalog
cat = fits.getdata(input['build_fname']('RawCatalogs',[input['query'],None]))

# loads the survey footprint in equatorial coordinates
footprint_res, footprint_zrange, sky_fraction, footprint = input['read_footprint']()
print("# this footprint covers {}% of the sky".format(sky_fraction*100.))

print("# selecting galaxies...")

zused = cat['true_redshift_gal']
redshift_sel = (zused>=footprint_zrange[0]) & (zused<=footprint_zrange[1])
ra_gal=cat['ra_gal'][redshift_sel]
dec_gal=cat['dec_gal'][redshift_sel]

conv=np.pi/180.
# if input.rgal is not None:
#     print("# rotating catalog...")
#     theta_eq, phi_eq = input.rgal( np.pi/2 - dec_gal * conv, ra_gal * conv )
# else:
theta_eq = np.pi/2. - dec_gal * conv
phi_eq = ra_gal * conv

## Get galaxy pixels in the sky
print("# finding sky pixels...")

pix      = hp.ang2pix( footprint_res, theta_eq, phi_eq )
foot_sel = np.zeros_like(zused,dtype=bool)
fp_small = footprint[pix]
foot_sel[redshift_sel] = fp_small

Nextract = foot_sel.sum()

print(f"# Nextract={Nextract}")

store = zarr.open_group(input['master_fname'](), mode="w")
extract = store.empty(shape=(Nextract,), dtype=cat.dtype, chunks=(10000000,), name="catalog")

#extract = np.empty(Nextract, dtype=cat.dtype)

for field in cat.dtype.names:
   print(f"# Doing field {field}")
   extract[field] = cat[field][foot_sel]
#foot_sel] = zarr.copy(cat[foot_sel], extract, log=sys.stdout)
#=cat[field][foot_sel]

del cat

#print('# writing catalog to file {}/RawCatalogs/{}_{}.fits'.format(input['outdir'],input['query'],input['footprint_tag']))

#fits.writeto('{}/RawCatalogs/{}_{}.fits'.format(input['outdir'],input['query'],input['footprint_tag']), extract, overwrite=True)

print("# done!")

