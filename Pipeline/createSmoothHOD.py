################################################################################
### Authors: Tiago Castro, Pierluigi Monaco                                  ###
###                                                                          ###
################################################################################
from astropy.io import fits
import numpy as np
from colossus.cosmology import cosmology
from scipy.ndimage import gaussian_filter
import sys
import healpy as hp

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

print("# Running createSmoothHOD.py with {}".format(sys.argv[1]))

footprint_res, footprint_zrange, sky_fraction, footprint = input.read_footprint()
print("# The footprint covers {}% of the sky".format(100*sky_fraction))

# Import the Raw Catalog
print("# Reading indices from file {}...".format(input.indices_fname()))
rawcat_indices = zarr.open_group(input.indices_fname())['indices']

print("# Reading Catalog {}...".format(input.master_fname()))
root = zarr.open_group(input.master_fname())
rawcat = root['catalog']

print("# Extracting colums from the catalog...")
good = np.array(rawcat_indices) > 0
zs   = rawcat['true_redshift_gal'][rawcat_indices][good]

del rawcat_indices

flux = rawcat[input.flux_key][good]
lm   = rawcat['halo_lm'][good]
kind = rawcat['kind'][good]

del rawcat
del good

print("# Binning galaxies in flux...")

# Flux cuts from the minimum flux defined in input.py to
# two magnitude higher
Ncut = 400
DeltaFlux = 1.6
fcut = np.linspace(input.logflux_limit, input.logflux_limit+DeltaFlux, Ncut)

# Selecting the galaxies inside the fcut bins 
gal_fluxbin = np.digitize( flux, fcut )

# Redshift bins
z1=footprint_zrange[0]
z2=footprint_zrange[1]
zbins = np.linspace(z1,z2,round( ( z2-z1)/input.deltazbin + 1 ) )

# Volume bins
Vbins = 4*np.pi/3 * (input.cosmo.comovingDistance(0.0, zbins[1:])**3 - input.cosmo.comovingDistance(0.0, zbins[:-1])**3) * sky_fraction

# Halo mass Bins
Mbins = np.linspace( lm.min(), lm.max() )

# Number of Centrals and Satelites inside each bin of redshift and flux
Ncen    = np.zeros( (fcut.size, zbins.size-1, Mbins.size-1) )
Nsat    = np.zeros( (fcut.size, zbins.size-1, Mbins.size-1) )
Nhalos  = np.zeros( (zbins.size-1, Mbins.size-1) )
Ncen_g  = np.zeros( (fcut.size, zbins.size-1, Mbins.size-1) )
Nsat_g  = np.zeros( (fcut.size, zbins.size-1, Mbins.size-1) )
Nhal_g  = np.zeros( (zbins.size-1, Mbins.size-1) )

print("# Selecting centrals and satellites...")

cen  = kind==0

print("# Binning galaxies in mass...")
Nhalos, bins1, bins2 = np.histogram2d(zs[cen], lm[cen], bins=[zbins, Mbins])
Nhalos = (Nhalos.T/Vbins).T

print("# Constructing the SDHOD...")

# Working on the flux limited sample
for i in range(1, fcut.size + 1):

    print("    Working on flux cut {} of {}...".format(i,Ncut))

    # Reading kind array
    idsi = (gal_fluxbin == i)
    this_kind = kind[ idsi ]
    this_zs   = zs  [ idsi ]
    this_lm   = lm  [ idsi ]
    this_cen  = this_kind==0
    this_sat  = ~this_cen

    Ncen[i-1], bins1, bins2 = np.histogram2d(this_zs[this_cen], this_lm[this_cen], bins=[zbins, Mbins])
    Nsat[i-1], bins1, bins2 = np.histogram2d(this_zs[this_sat], this_lm[this_sat], bins=[zbins, Mbins])

    Ncen[i-1] = ( (Ncen[i-1].T) /Vbins ).T
    Nsat[i-1] = ( (Nsat[i-1].T) /Vbins ).T

    # Smoothing the dn/dz
    for j in range(Mbins.size-1):

        Ncen_g[i-1,:,j] = gaussian_filter( Ncen[i-1,:,j], input.smoothing_length )
        Nsat_g[i-1,:,j] = gaussian_filter( Nsat[i-1,:,j], input.smoothing_length ) 


print("# Smoothing the halo SDHOD...")

# Smoothing the halo dn/dz
for j in range(Mbins.size-1):
    Nhal_g[:,j] = gaussian_filter( Nhalos[:,j], input.smoothing_length )


print("# Cumulating the number of galaxies...")

# Cumulative Summation on the fluxes 
Ncen   = Ncen.sum(axis=0) - Ncen.cumsum(axis=0) + Ncen
Nsat   = Nsat.sum(axis=0) - Nsat.cumsum(axis=0) + Nsat
Ncen_g = Ncen_g.sum(axis=0) - Ncen_g.cumsum(axis=0) + Ncen_g
Nsat_g = Nsat_g.sum(axis=0) - Nsat_g.cumsum(axis=0) +  Nsat_g


print("# Writing file {}...".format(input.SDHOD_fname()))

# Write SDHOD in fits file
cat = np.empty(1,  dtype=[('n_halos', float, (zbins.size-1, Mbins.size-1)),
                          ('n_cen', float, (fcut.size, zbins.size-1, Mbins.size-1)), 
                          ('n_sat', float, (fcut.size, zbins.size-1, Mbins.size-1)),
                          ('n_hal_gaus', float, (zbins.size-1, Mbins.size-1)),
                          ('n_cen_gaus', float, (fcut.size, zbins.size-1, Mbins.size-1)), 
                          ('n_sat_gaus', float, (fcut.size, zbins.size-1, Mbins.size-1)),
                          ('z_bins', float, zbins.size),
                          ('V_bin', float, Vbins.size),
                          ('M_bins', float, Mbins.size),
                          ('f_bins', float, fcut.size)])

cat['n_halos'] = Nhalos
cat['n_cen'] = Ncen;
cat['n_sat'] = Nsat;
cat['n_hal_gaus'] = Nhal_g;
cat['n_cen_gaus'] = Ncen_g;
cat['n_sat_gaus'] = Nsat_g;
cat['z_bins'] = zbins; 
cat['V_bin'] = Vbins;
cat['M_bins'] = Mbins;
cat['f_bins'] = fcut;

fits.writeto(input.SDHOD_fname(), cat, overwrite=True)

print("# done!")
