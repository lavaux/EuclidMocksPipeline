################################################################################
### Authors: Tiago Castro, Pierluigi Monaco, Guilhem Lavaux                  ###
###                                                                          ###
################################################################################
from . import register_tool
from ..config import readConfig
import sys
from typing import Tuple
from numpy.typing import *


def nside2norder(nside):
    """
    Give the HEALpix order for the given HEALpix nside parameter.

    :param nside: nside of the healpy pixelization
    :return: norder, norder of the healpy pixelization
    """
    import numpy as np

    norder = np.log2(nside)
    if not (norder.is_integer()):
        raise ValueError("Wrong nside number (it is not 2**norder)")
    return int(norder)


def rand_vec_in_pix(nside, ipix, nest=False):
    """
    Draw vectors from a uniform distribution within a HEALpixel.

    (Taken from astrotools: https://git.rwth-aachen.de/astro/astrotools/-/blob/master/astrotools/healpytools.py

    :param nside: nside of the healpy pixelization
    :param ipix: pixel number(s)
    :param nest: set True in case you work with healpy's nested scheme
    :return: vectors containing events from the pixel(s) specified in ipix
    """
    import healpy as hp
    import numpy as np

    if not nest:
        ipix = hp.ring2nest(nside, ipix=ipix)

    n_order = nside2norder(nside)
    n_up = 29 - n_order
    i_up = ipix * 4 ** n_up
    i_up += np.random.randint(0, 4 ** n_up, size=np.size(ipix))

    return np.asarray(hp.pix2vec(nside=2 ** 29, ipix=i_up, nest=True)).transpose()


def generate_random_simple(input: dict) -> Tuple[ArrayLike,ArrayLike]:
    from ..disk import DefaultCatalogRead, DefaultCatalogWrite
    import numpy as np
    from ..filenames import galcat, selection_data

    fname = galcat(input, input.pinocchio_first_run)
    print ("# loading data catalog {}...".format(fname))
    with DefaultCatalogRead(fname) as gal_store:
        datacat = gal_store['catalog']

        # apply selection to data catalog if required
        if (input.apply_dataselection_to_random) & (input.selection_data_tag is not None):
            print("# Applying selection {input.selection_data_tag} to data...")
            with DefaultCatalogRead(selection_data(input)) as sel_store:
                sel = sel_store['SELECTION']['SELECTION']
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

    return redshift, flux

def generate_random_pinocchio(input: dict) -> Tuple[ArrayLike, ArrayLike]:
    import numpy as np
    from ..disk import DefaultCatalogRead
    from ..filenames import galcat, dndz

    # for a set of data catalogs, it adds them randomly to the random vector.
    # Here each data mock is used as a whole, and it can be replicated several times
    with DefaultCatalogRead(dndz(r1=input.pinocchio_first_run,r2=input.pinocchio_last_run)) as dn_store:
        dndz  = dn_store['dn_dz']
    Ngal  = np.int(dndz['N_gal'].sum())
    Nrandom = np.int(input.alpha * Ngal)

    # we extract here alpha+1 mocks to avoid that the number of galaxies is insufficient
    tobeused = np.sort(np.random.uniform(input.pinocchio_first_run,input.pinocchio_last_run+1,input.alpha+1).astype(int))
    howmanytimes = np.array([np.in1d(tobeused,i).sum() for i in range(0,input.pinocchio_last_run+1)])

    Nstored  = 0
    redshift = np.empty(Nrandom, dtype=float)
    flux     = np.empty(Nrandom, dtype=float)
    for myrun in np.arange(input.pinocchio_first_run, input.pinocchio_last_run + 1):
        if howmanytimes[myrun]>0:

            fname = galcat(myrun)
            print ("# loading data catalog {} to be used {} times...".format(fname,howmanytimes[myrun]))
            with DefaultCatalogRead(fname) as galstore:
                datacat = galstore['catalog']
            
            # applies selection if required
            if (input.apply_dataselection_to_random) & (input.selection_data_tag is not None):
                print("# Applying selection {} to data...".format(input.selection_data_tag))
                with DefaultCatalogRead(selection_data(input, myrun)) as selstore:
                    sel = sel['SELECTION']['SELECTION']
                data_redshift = datacat[input.redshift_key][sel]
                data_flux     = datacat[input.flux_key][sel]
            else:
                data_redshift = datacat[input.redshift_key]
                data_flux     = datacat[input.flux_key]
            del datacat

            Ndata = data_redshift.size
            # no selection is applied to pinocchio mocks

            for i in range(howmanytimes[myrun]):
                redshift = np.concatenate((redshift,data_redshift))
                flux     = np.concatenate((flux,data_flux))

                Nrandom+=Ndata
                print("    added %d random galaxies, total: %d"%(Ndata,Nrandom))
    
    return redshift, flux

@register_tool
def createRandom(config: str, legacy_algorithm: bool = False) -> None:
    """Creates a random for the galaxy catalog

    The random catalog is created  with abundance alpha times the data catalog, alpha as
    specified in the input catalog. TO BE IMPLEMENTED: on request it can
    follow a specified galaxy number density.

    Args:
        config (str): Pipeline config file
    """
    from ..disk import DefaultCatalogRead, DefaultCatalogWrite
    from .. import filenames
    import healpy as hp
    import numpy as np

    input = readConfig(config)

    print(f"# Running createRandom.py with {config}")

    np.random.seed(seed=input.SEED_random)

    # TO BE IMPLEMENTED: read smoothed dndz and force the random to follow it

    # loads the footprint
    footprint_res, footprint_zrange, sky_fraction, footprint = input.read_footprint()

    print("# Computing angular extent of footprint...")
    theta, phi = hp.pix2ang(footprint_res, np.arange(hp.nside2npix(footprint_res)))
    meanspacing = np.sqrt(4.0 * np.pi / footprint.size)

    thmin = np.min(theta[footprint]) - 2.0 * meanspacing
    thmax = np.max(theta[footprint]) + 2.0 * meanspacing
    phmin = np.min(phi[footprint]) - 2.0 * meanspacing
    phmax = np.max(phi[footprint]) + 2.0 * meanspacing

    print(
        "# I will populate an area with theta=[%f,%f], phi=[%f,%f]"
        % (thmin, thmax, phmin, phmax)
    )

    print("# Assigning redshifts and fluxes...")
    if (input.cat_type is not 'pinocchio') | (input.pinocchio_last_run is None):
        redshift,flux=generate_random_simple(input)
    else:
        redshift,flux=generate_random_pinocchio(input)

    Nrandom = redshift.size
    print("# Saving random to file {}...".format(filenames.random(input)))

    with DefaultCatalogWrite(filenames.random(input)) as store:
        catalog = store.new_array(
            "catalog",
            shape=(Nrandom,),
            dtype=[
                (input.redshift_key, float),
                ("ra_gal", float),
                ("dec_gal", float),
                (input.my_flux_key, float),
            ],
            tags={"legacy": legacy_algorithm},
        )

        ra_gal = catalog["ra_gal"]
        dec_gal = catalog["dec_gal"]
        print(f"footprint sum = {footprint.sum()}")

        catalog[input.redshift_key] = redshift
        catalog[input.flux_key] = flux

        print("# Starting to create {} random galaxies...".format(Nrandom))
        Nbunch = Nrandom // 10

        Nstored = 0

        selected_pixels = np.where(footprint == 1)[0]
        while Nstored < Nrandom:
            # pick a random pixel amont the footprint pixels that are 1
            pix_list = selected_pixels[
                np.random.randint(low=0, high=selected_pixels.size, size=Nbunch)
            ]

            # Pick a random vector inside each pixel
            randvec = rand_vec_in_pix(
                nside=footprint_res, ipix=pix_list, nest=False
            )
            randdec, randra = hp.vec2ang(randvec)
            Nselect = Nbunch

            Nup2now = Nstored + Nselect
            if Nup2now > Nrandom:
                Nup2now = Nrandom
                Nselect = Nrandom - Nstored
            print(
                "    selected %d random galaxies out of %d, total: %d (target %d)"
                % (Nselect, Nbunch, Nup2now, Nrandom)
            )
            ra_gal[Nstored:Nup2now] = randra[:Nselect]
            dec_gal[Nstored:Nup2now] = np.pi / 2.0 - randdec[:Nselect]

            Nstored = Nup2now

        ra_gal *= 180.0 / np.pi
        dec_gal *= 180.0 / np.pi

    
    print("# DONE!")
