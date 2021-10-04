################################################################################
### Author: Tiago Castro, Pierluigi Monaco                                   ###
###                                                                          ###
################################################################################
from . import register_tool
from ..config import readConfig


@register_tool
def numbercounts(config: str):
    """Build number counts per magnitude bin

    Args:
      config (str): Configuration
    """
    import numpy as np
    from astropy.io import fits
    import healpy as hp
    import sys
    from os import path
    from ..disk import DefaultCatalogRead, DefaultCatalogWrite

    print("# Running numbercounts.py with {}".format(config))
    input = readConfig(config)

    myrun=None
    if len(sys.argv)>=3:
        if sys.argv[2].isdigit and input.cat_type is 'pinocchio':
            myrun=int(sys.argv[2])
            print("# I will process run number {}".format(myrun))
        else:
            print("# WARNING: unrecognised command-line option {}".format(sys.argv[2]))
    else:
        if input.cat_type is 'pinocchio':
            myrun=input.pinocchio_first_run
            print("# I will process run number {}".format(myrun))
    
    # data catalog
    fname = input.galcat_fname(myrun)

    print("# loading catalog {}...".format(fname))

    if not path.exists(fname):
        print("ERROR: galaxy catalog {} does not exist".format(fname))
        sys.exit(-1)

    with DefaultCatalogRead(fname) as store:
        cat = store["catalog"]

    # loads the survey footprint in equatorial coordinates
    footprint_res, footprint_zrange, sky_fraction, footprint = input.read_footprint()
    sky_coverage = input.sqdegonthesky * sky_fraction

    # selection
    if input.selection_data_tag is not None:
        print("# loading selection {}...".format(input.selection_data_fname()))
        with DefaultCatalogRead(input.selection_data_fname(run=myrun)) as store:
            selection = store["SELECTION"]["SELECTION"]
    else:
        selection = np.ones(len(cat), dtype=bool)

    print("# setting binnings...")

    nzbins = len(input.finalCatZShell)
    zbins = np.asarray(input.finalCatZShell)

    deltaf = 0.1
    diff = -14.0 - input.logflux_limit
    nflux = int(diff / deltaf) + 1
    fluxmax = input.logflux_limit + deltaf * nflux
    flux_bins = np.linspace(input.logflux_limit, fluxmax, num=nflux + 1)
    xflux = 10.0 ** (0.5 * (flux_bins[1:] + flux_bins[:-1]))
    DeltaF = 2.5 * (flux_bins[1:] - flux_bins[:-1])

    for i in range(nzbins):

        z1 = zbins[i, 0]
        z2 = zbins[i, 1]

        print("# Processing redshift interval [{},{}]".format(z1, z2))

        zsel = (cat[input.redshift_key] >= z1) & (cat[input.redshift_key] < z2)
        Ngal = np.histogram(cat[input.flux_key][selection & zsel], bins=flux_bins)[0]
        isCen = cat["kind"][selection & zsel] == 0
        Ncen = np.histogram(
            cat[input.flux_key][selection & zsel][isCen], bins=flux_bins
        )[0]

        LF = Ngal / sky_coverage / (z2 - z1) / DeltaF
        LF_cen = Ncen / sky_coverage / (z2 - z1) / DeltaF
        LF_sat = (Ngal - Ncen) / sky_coverage / (z2 - z1) / DeltaF

        ## Writes on file
        fname = input.numbercounts_fname(z1, z2, run=myrun)
        print("# Writing results in file {}".format(fname))

        with DefaultCatalogWrite(fname) as store:
            flux_counts = store.new_array(
                "NUMBER_COUNT",
                shape=(Ngal.size,),
                dtype=[
                    ("LF", np.float),
                    ("LF_cen", np.float),
                    ("LF_sat", np.float),
                    ("f_center", np.float),
                    ("f_lower", np.float),
                    ("f_upper", np.float),
                ],
            )

            flux_counts["LF"] = LF
            flux_counts["LF_cen"] = LF_cen
            flux_counts["LF_sat"] = LF_sat
            flux_counts["f_center"] = xflux
            flux_counts["f_lower"] = 10.0 ** flux_bins[:-1]
            flux_counts["f_upper"] = 10.0 ** flux_bins[1:]

    print("# done!")
