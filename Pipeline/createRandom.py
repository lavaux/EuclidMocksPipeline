################################################################################
### Authors: Tiago Castro, Pierluigi Monaco                                  ###
###                                                                          ###
################################################################################
import numpy as np
import healpy as hp
import sys
from astropy.io import fits
import filenames

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

print("# Assigning redshifts and fluxes...")

if (input.cat_type is not 'pinocchio') | (input.pinocchio_last_run is None):

    fname = filenames.galcat(input,input.pinocchio_first_run)
    print ("# loading data catalog {}...".format(fname))
    datacat = fits.getdata(fname)

    # apply selection to data catalog if required
    if (input.apply_dataselection_to_random) & (input.selection_data_tag is not None):
        print("# Applying selection {} to data...".format(input.selection_data_tag))
        sel = fits.getdata(filenames.selection_data(input))['SELECTION']
        data_redshift = datacat[input.redshift_key][sel]
        data_flux     = datacat[input.flux_key][sel]
    else:
        data_redshift = datacat[input.redshift_key]
        data_flux     = datacat[input.flux_key]

    del datacat

    Ngal = data_redshift.size
    Nrandom = np.int(input.alpha * Ngal)

    # for a single data catalog, it adds galaxies randomly from the catalog
    # until the wanted number is reached, so the sequence of galaxies is not 
    # simply alpha replications of the data catalog
    Nbunch  = Nrandom // 10
    Nstored = 0
    redshift = np.empty(Nrandom, dtype=float)
    flux     = np.empty(Nrandom, dtype=float)
    while (Nstored<Nrandom):

        Nadd = Nbunch
        if Nbunch + Nstored > Nrandom:
            Nadd = Nrandom - Nstored

        thesegals = np.random.uniform(0,Ngal,Nadd).astype(int)
        redshift[Nstored:Nstored+Nadd] = data_redshift[thesegals]
        flux[Nstored:Nstored+Nadd]     = data_flux[thesegals]
        Nstored += Nadd
        print("    added %d random galaxies, total: %d"%(Nadd,Nstored))

else:

    # for a set of data catalogs, it adds them randomly to the random vector.
    # Here each data mock is used as a whole, and it can be replicated several times

    # we extract here alpha+1 mocks to avoid that the number of galaxies is insufficient
    tobeused = np.sort(np.random.uniform(input.pinocchio_first_run,input.pinocchio_last_run+1,input.alpha+1).astype(int))
    howmanytimes = np.array([np.in1d(tobeused,i).sum() for i in range(0,input.pinocchio_last_run+1)])

    Nrandom=0
    redshift = np.array([])
    flux     = np.array([])

    for myrun in np.arange(input.pinocchio_first_run, input.pinocchio_last_run + 1):
        if howmanytimes[myrun]>0:

            fname = filenames.galcat(input,myrun)
            print ("# loading data catalog {} to be used {} times...".format(fname,howmanytimes[myrun]))
            datacat = fits.getdata(fname)

            # applies selection if required
            if (input.apply_dataselection_to_random) & (input.selection_data_tag is not None):
                print("# Applying selection {} to data...".format(input.selection_data_tag))
                sel = fits.getdata(filenames.selection_data(input,myrun))['SELECTION']
                data_redshift = datacat[input.redshift_key][sel]
                data_flux     = datacat[input.flux_key][sel]
            else:
                data_redshift = datacat[input.redshift_key]
                data_flux     = datacat[input.flux_key]
            del datacat

            Ndata = data_redshift.size

            for i in range(howmanytimes[myrun]):
                redshift = np.concatenate((redshift,data_redshift))
                flux     = np.concatenate((flux,data_flux))

                Nrandom+=Ndata
                print("    added %d random galaxies, total: %d"%(Ndata,Nrandom))

print("# Starting to create {} random positions...".format(Nrandom))
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

print("# Saving random to file {}...".format(filenames.random(input)))
################# Saving random galaxies and random redshift################
red  = fits.Column( name=input.redshift_key,  array=redshift, format='E' )
ra   = fits.Column( name='ra_gal',            array=ra_gal,   format='E' )
dec  = fits.Column( name='dec_gal',           array=dec_gal,  format='E' )
flu  = fits.Column( name=input.flux_key,   array=flux,     format='E' )
tw   = fits.BinTableHDU.from_columns([red, ra, dec, flu])
tw.writeto(filenames.random(input), overwrite=True)

print("# DONE!")
