################################################################################
### Authors: Tiago Castro, Pierluigi Monaco                                  ###
###                                                                          ###
################################################################################
import numpy as np
import healpy as hp
import sys
from astropy.io import fits


if len(sys.argv)<2:
    print("Usage: pyton {} [my input file]".format(sys.argv[0]))
    sys.exit(0)
try: 
    input = __import__(sys.argv[1],  globals(), locals(), [], 0)
except ModuleNotFoundError:
    print("input file not found")
    print("Usage: pyton {} [my input file]".format(sys.argv[0]))
    sys.exit(0)


print("# Running createRandom.py with {}".format(sys.argv[1]))

np.random.seed(seed=input.SEED_random)

# TO BE IMPLEMENTED: read smoothed dndz and force the random to follow it

# loads the footprint
footprint_res, footprint_zrange, sky_fraction, footprint = input.read_footprint()

print("# Computing angular extent of footprint...")
theta,phi = hp.pix2ang(footprint_res,np.arange(hp.nside2npix(footprint_res)))
meanspacing = np.sqrt(4.*np.pi / footprint.size)

thmin = np.min(theta[footprint]) - 2.*meanspacing
thmax = np.max(theta[footprint]) + 2.*meanspacing
phmin = np.min(phi[footprint])   - 2.*meanspacing
phmax = np.max(phi[footprint])   + 2.*meanspacing

print("# I will populate an area with theta=[%f,%f], phi=[%f,%f]"%(thmin,thmax,phmin,phmax))

print ("# loading data catalog {}...".format(input.galcat_fname()))
datacat       = fits.getdata(input.galcat_fname())
data_redshift = datacat[input.redshift_key]
data_flux     = datacat[input.my_flux_key]
del datacat

Ngal = data_redshift.size
Nrandom = np.int(input.alpha * Ngal)

print("# Starting to create {} random galaxies...".format(Nrandom))
Nbunch = Nrandom//10

Nstored = 0
ra_gal  = np.empty(Nrandom, dtype=float)
dec_gal = np.empty(Nrandom, dtype=float)

while (Nstored<Nrandom):
    randra   = np.random.uniform(phmin,phmax,Nbunch)
    randdec  = np.arccos( np.random.uniform( np.cos(thmin), np.cos(thmax), Nbunch) )

    pix = hp.ang2pix(footprint_res, randdec, randra)

    select  = footprint[pix]
    Nselect = select.sum()
    Nup2now = Nstored + Nselect
    if Nup2now>Nrandom:
        Nup2now = Nrandom
        Nselect = Nrandom - Nstored
    print("    selected %d random galaxies out of %d, total: %d"%(Nselect,Nbunch,Nup2now))
    ra_gal[Nstored:Nup2now]  = randra[select][:Nselect]
    dec_gal[Nstored:Nup2now] = np.pi/2. - randdec[select][:Nselect]

    Nstored = Nup2now

ra_gal *= 180./np.pi
dec_gal *= 180./np.pi

print("# Assigning redshifts and fluxes...")
Nbunch  = Nrandom // 10
Nstored = 0
redshift = np.empty(Nrandom, dtype=float)
flux     = np.empty(Nrandom, dtype=float)
while (Nstored<Nrandom):

    Nadd = Nbunch
    if Nbunch + Nstored > Nrandom:
        Nadd = Nrandom - Nstored

    thesegals = np.random.uniform(0,Ngal,Nadd).astype(int)
    # This is experimental...
    if input.smooth_dndz_in_random:
        redshift[Nstored:Nstored+Nadd] = data_redshift[thesegals] * np.random.normal(1., input.deltazbin * input.smoothing_length, Nadd)
    else:
        redshift[Nstored:Nstored+Nadd] = data_redshift[thesegals]
    flux[Nstored:Nstored+Nadd]     = data_flux[thesegals]
    Nstored += Nadd
    print("    added %d random galaxies, total: %d"%(Nadd,Nstored))


print("# Saving random to file {}...".format(input.random_fname()))
################# Saving random galaxies and random redshift################
red  = fits.Column( name=input.redshift_key,  array=redshift, format='E' )
ra   = fits.Column( name='ra_gal',            array=ra_gal,   format='E' )
dec  = fits.Column( name='dec_gal',           array=dec_gal,  format='E' )
flu  = fits.Column( name=input.my_flux_key,   array=flux,     format='E' )
tw   = fits.BinTableHDU.from_columns([red, ra, dec, flu])
tw.writeto(input.random_fname(), overwrite=True)

print("# DONE!")
