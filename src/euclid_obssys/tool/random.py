################################################################################
### Authors: Tiago Castro, Pierluigi Monaco                                  ###
###                                                                          ###
################################################################################
from . import register_tool
from ..config import readConfig
import sys



def nside2norder(nside):
    """
    Give the HEALpix order for the given HEALpix nside parameter.

    :param nside: nside of the healpy pixelization
    :return: norder, norder of the healpy pixelization
    """
    import numpy as np

    norder = np.log2(nside)
    if not (norder.is_integer()):
        raise ValueError('Wrong nside number (it is not 2**norder)')
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


@register_tool
def createRandom(config: str) -> None:
    """Creates a random for the galaxy catalog

    The random catalog is created  with abundance alpha times the data catalog, alpha as
    specified in the input catalog. TO BE IMPLEMENTED: on request it can
    follow a specified galaxy number density.

    Args:
        config (str): Pipeline config file
    """
    from euclid_obssys.disk import DefaultCatalogRead, DefaultCatalogWrite
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

    print(f"# loading data catalog {input.galcat_fname()}...")
    with DefaultCatalogRead(input.galcat_fname()) as store:
        datacat = store["catalog"]
        data_redshift = datacat[input.redshift_key]
        data_flux = datacat[input.my_flux_key]
        del datacat

    Ngal = data_redshift.size
    Nrandom = np.int(input.alpha * Ngal)

    print("# Starting to create {} random galaxies...".format(Nrandom))
    Nbunch = Nrandom // 10

    Nstored = 0

    with DefaultCatalogWrite(input.random_fname()) as store:
        catalog = store.new_array(
            "catalog",
            shape=(Nrandom,),
            dtype=[
                (input.redshift_key, float),
                ("ra_gal", float),
                ("dec_gal", float),
                (input.my_flux_key, float),
            ],
        )

        ra_gal = catalog["ra_gal"]
        dec_gal = catalog["dec_gal"]
        redshift = catalog[input.redshift_key]
        flux = catalog[input.flux_key]
        print(f"footprint sum = {footprint.sum()}")

        selected_pixels = np.where(footprint==1)[0]
        while Nstored < Nrandom:
            #pick a random pixel amont the footprint pixels that are 1
            pix_list = selected_pixels[np.random.randint(low=0, high=selected_pixels.size, size=Nbunch)]

            # Pick a random vector inside each pixel
            randvec = rand_vec_in_pix(nside=footprint_res, ipix=pix_list, nest=False)
            randdec,randra = hp.vec2ang(randvec)

            #randra = np.random.uniform(phmin, phmax, Nbunch)
            #randdec = np.arccos(np.random.uniform(np.cos(thmin), np.cos(thmax), Nbunch))

            #pix = hp.ang2pix(footprint_res, randdec, randra)

            #select = footprint[pix]
            #Nselect = select.sum()
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

        print("# Assigning redshifts and fluxes...")
        Nbunch = Nrandom // 10
        Nstored = 0
        while Nstored < Nrandom:

            Nadd = Nbunch
            if Nbunch + Nstored > Nrandom:
                Nadd = Nrandom - Nstored

            thesegals = np.random.uniform(0, Ngal, Nadd).astype(int)
            # This is experimental...
            if input.smooth_dndz_in_random:
                redshift[Nstored : Nstored + Nadd] = data_redshift[
                    thesegals
                ] * np.random.normal(
                    1.0, input.deltazbin * input.smoothing_length, Nadd
                )
            else:
                redshift[Nstored : Nstored + Nadd] = data_redshift[thesegals]
            flux[Nstored : Nstored + Nadd] = data_flux[thesegals]
            Nstored += Nadd
            print("    added %d random galaxies, total: %d" % (Nadd, Nstored))

        print("# Saving random to file {}...".format(input.random_fname()))

    print("# DONE!")
