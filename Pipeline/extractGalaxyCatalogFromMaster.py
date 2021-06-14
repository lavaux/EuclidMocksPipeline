################################################################################
### Author: Tiago Castro, Pierluigi Monaco                                   ###
###                                                                          ###
################################################################################
import numpy as np
from astropy.io import fits
import sys
from os import path
from euclid_obssys.config import readConfig
import zarr


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

print("# Running extractGalaxyCatalogFromMaster.py with {}".format(sys.argv[1]))

footprint_res, footprint_zrange, sky_fraction, footprint = input['read_footprint']()
del footprint

print("# loading catalog {}...".format(input['master_fname']()))

# input raw catalog
cat = zarr.open_group(input['master_fname'](), mode="r")['catalog']

print("# selecting galaxies...")

selection = cat[input['flux_key']] > input['logflux_limit']

Nextract = selection.sum()

print("# extracting {} galaxies".format(Nextract))

extract = np.empty(Nextract, 
                   dtype=[('x_gal', float), ('y_gal', float), ('z_gal', float),
                          ('ra_gal', float), ('dec_gal', float), ('kind', int), 
                          ('true_redshift_gal', float), ('observed_redshift_gal', float), 
                          ('halo_lm', float), ('galaxy_id', int), ('halo_id', int), 
                          (input['flux_key'], float), ('sh_'+input['flux_key'], float)])

for field in extract.dtype.names:
    print("    processing {}".format(field))

    if field=='sh_'+input['flux_key']:
        print("# shuffling galaxies...")
        
        shuffled=np.copy(cat[input['flux_key']][selection])

        zbins = np.linspace(footprint_zrange[0],footprint_zrange[1],round( (footprint_zrange[1]-footprint_zrange[0] )/input['deltazbin'] + 1 ) )
        zindex = np.digitize( extract['true_redshift_gal'], zbins )
        for iz in np.arange(zbins.size-1):
            ff = zindex==iz
            this = shuffled[ff]
            np.random.shuffle(this)
            shuffled[ff] = this
        extract[field]=shuffled
    else:

        extract[field]=cat[field][selection]

del cat

fname=input['flagcat_fname']()

print("# writing file {}".format(fname))
out = zarr.open_group(fname, mode="w")
out['catalog'] = extract
#fits.writeto(fname, extract, overwrite=True)

print("# done!")
