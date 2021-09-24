################################################################################
### Authors: Tiago Castro, Pierluigi Monaco, Guilhem Lavaux                  ###
###                                                                          ###
################################################################################
from . import register_tool
from ..config import readConfig
import sys
from typing import Optional


@register_tool
def createSelection(
    config: str, run_number: Optional[int] = None, use_data: bool = True
) -> None:
    import numpy as np
    from astropy.io import fits
    import healpy as hp
    import sys
    import os
    from scipy.interpolate import RegularGridInterpolator
    from euclid_obssys.disk import DefaultCatalogRead, DefaultCatalogWrite

    def extinction_curve(z):

        if sel_input.extinction_curve is "standard":
            wavel = 6562.8 * (1.0 + z)
            mwl = np.array(
                [
                    3600.0,
                    4500.0,
                    5500.0,
                    6600.0,
                    8000.0,
                    12500.0,
                    16500.0,
                    22000.0,
                    35000.0,
                    48000.0,
                ]
            )
            mwA = np.array(
                [1.531, 1.324, 1.000, 0.748, 0.482, 0.282, 0.175, 0.112, 0.058, 0.023]
            )
            return 3.086 * np.interp(wavel, mwl, mwA)
        else:
            return None

    input = readConfig(config)

    print(f"# Running createSelection.py with {config}")

    # check if data or random catalog is to be selected
    myrun = None
    if use_data:
        print("# Selection will be applied to the random catalog")

    if run_number is not None and input.cat_type != "pinocchio":
        print("# WARNING: run_number used for non Pinocchio catalog.")

    if input.cat_type == "pinocchio":
        if run_number is not None:
            myrun = run_number
        else:
            myrun = input.pinocchio_first_run

        print("# I will process run number {}".format(myrun))

    if use_data:
        print("# Selection will be applied to the data catalog")

    # check selection tags
    if use_data:
        if input.selection_data_tag is None:
            print("No selection specified for galaxy catalog, exiting")
            sys.exit(0)
        else:
            sel_input_fname = "sel_input_{}.py".format(input.selection_data_tag)
    else:
        if input.selection_random_tag is None:
            print("No selection specified for random catalog, exiting")
            sys.exit(0)
        else:
            sel_input_fname = "sel_input_{}.py".format(input.selection_random_tag)

        if input.apply_dataselection_to_random:
            print(
                "apply_dataselection_to_random has been specified, nothing to do here, exiting"
            )
            sys.exit(0)

    # load selection input
    sel_input = readConfig(sel_input_fname)

    # list of allowed selection types
    allowed_selections = [
        "extinction",
        "visibilitymask",
        "fluxcut",
        "central",
        "satellite",
        "lookup",
    ]
    for key in sel_input.selection_keys:
        if key not in allowed_selections:
            print("Error: selection key {} not recognised".format(key))
            sys.exit(1)

    # load galaxy catalog
    if myrun is None:
        myrun = input.cat_type

    if use_data:
        cat_fname = input.galcat_fname(myrun)
        print(f"Opening galaxy catalog {cat_fname}...")
    else:
        cat_fname = input.random_fname()
        print(f"Opening random catalog {cat_fname}...")

    with DefaultCatalogRead(cat_fname) as store:
        cat = store["catalog"]

    Ngal = len(cat)

    selection = np.ones((Ngal), dtype=bool)

    print("The catalog contains {} galaxies".format(Ngal))

    if "lookup" in sel_input.selection_keys:

        # this option is alternative to extinction, visibilitymask, fluxcut

        print(
            "Loading lookup table {}".format(
                input.outdir + sel_input.lookup_table_fname
            )
        )
        lut = fits.getdata(os.path.join(input.outdir, sel_input.lookup_table_fname))

        # THESE SHOULD BE IN THE HEADER...
        nred = 14
        nHal = 14
        nexp = 4
        nbkg = 13

        ndata = nred * nHal * nexp * nbkg
        if ndata != len(lut):
            print("Error: table length is %d while I expected %d" % (ndata, len(lut)))
            sys.exit(0)

        print(
            "Look-up table dimension: {} * {} * {} * {} = {}".format(
                nred, nHal, nexp, nbkg, ndata
            )
        )

        xred = lut["z"].reshape(nbkg, nexp, nHal, nred)[0, 0, 0, :]
        xHal = np.log10(lut["flux_Ha"]).reshape(nbkg, nexp, nHal, nred)[0, 0, :, 0]
        xexp = lut["exp_time"].reshape(nbkg, nexp, nHal, nred)[0, :, 0, 0]
        xbkg = lut["sky_bg"].reshape(nbkg, nexp, nHal, nred)[:, 0, 0, 0]

        one_exposure = xexp[0]

        prob = lut["galaxy_frac"].reshape(nbkg, nexp, nHal, nred)

        interpolator = RegularGridInterpolator(
            (xbkg, xexp, xHal, xred),
            prob,
            method="linear",
            bounds_error=False,
            fill_value=None,
        )

        rand_numbers = np.random.rand(Ngal)

        pixels = None

        # applies extinction
        my_flux = cat[input.flux_key]
        if "extinction" in sel_input.selection_keys:

            print("# applying extinction in lookup table...")
            conv = np.pi / 180.0
            fname = os.path.join(input.outdir, sel_input.extinctionmap_fname)
            print("# loading reddening map {}...".format(fname))
            reddening = hp.read_map(fname, field=sel_input.extinctionmap_field)
            if sel_input.extinctionmap_res != hp.npix2nside(reddening.size):
                print("# resampling extinction map...")
                reddening = hp.ud_grade(reddening, sel_input.extinctionmap_res)
            print("# finding pixels for galaxies...")
            # NB pixels in principle can be computed only once by resampling one of the maps
            pix = hp.ang2pix(
                sel_input.extinctionmap_res,
                np.pi / 2.0 - cat["dec_gal"] * conv,
                cat["ra_gal"] * conv,
            )
            print("# computing extinction...")
            my_flux -= 0.4 * reddening[pix] * extinction_curve(cat[input.redshift_key])

        # background noise
        if sel_input.lookup_noise_fname is None:
            print("# Background set to minimal")
            background = xbkg[0] * np.ones(Ngal, dtype=np.float)  # minimal background
        else:
            if os.path.isfile(input.outdir + sel_input.lookup_noise_fname):
                print(
                    "# Reading noise map {}".format(
                        input.outdir + sel_input.lookup_noise_fname
                    )
                )
                noise_map = hp.read_map(
                    input.outdir + sel_input.lookup_noise_fname, partial=True
                )
                if pixels is None:
                    conv = np.pi / 180.0
                    pixels = hp.ang2pix(
                        sel_input.lookup_Nside,
                        np.pi / 2.0 - cat["dec_gal"][selection] * conv,
                        cat["ra_gal"][selection] * conv,
                    )
                background = noise_map[pixels]
            else:
                print(
                    "ERROR: noise map file {} not found".format(
                        input.outdir + sel_input.lookup_noise_fname
                    )
                )
                sys.exit(0)

        # exposure time
        if sel_input.lookup_exptime_fname is None:
            print("# Exposure map set to 4 exposures for all")
            exposure = (
                4.0 * one_exposure * np.ones(Ngal, dtype=np.float)
            )  # four exposures
        else:
            if os.path.isfile(input.outdir + sel_input.lookup_exptime_fname):
                print(
                    "# Reading exposure time map {}".format(
                        input.outdir + sel_input.lookup_exptime_fname
                    )
                )
                exptime_map = hp.read_map(
                    input.outdir + sel_input.lookup_exptime_fname, partial=True
                )
                if pixels is None:
                    conv = np.pi / 180.0
                    pixels = hp.ang2pix(
                        sel_input.lookup_Nside,
                        np.pi / 2.0 - cat["dec_gal"][selection] * conv,
                        cat["ra_gal"][selection] * conv,
                    )
                exposure = exptime_map[pixels] * one_exposure
            else:
                print(
                    "ERROR: exposure time map file {} not found".format(
                        input.outdir + sel_input.lookup_exptime_fname
                    )
                )
                sys.exit(0)

        points = np.array(
            [background, exposure, my_flux, cat[input.redshift_key]]
        ).transpose()
        selection &= interpolator(points) > rand_numbers

    if "flux_limit" in sel_input.selection_keys:
        my_flux_limit = sel_input.selection_logflux_limit
    else:
        my_flux_limit = input.logflux_limit

    if ("extinction" in sel_input.selection_keys) & (
        "lookup" not in sel_input.selection_keys
    ):

        print("# applying extinction...")

        conv = np.pi / 180.0

        # THIS PROCESS MAY BE HEAVY FOR THE RANDOM CATALOG, WE MAY SPLIT IT INTO SEVERAL SECTIONS

        fname = input.outdir + sel_input.extinctionmap_fname
        print("# loading reddening map {}...".format(fname))
        reddening = hp.read_map(fname, field=sel_input.extinctionmap_field)
        if sel_input.extinctionmap_res != hp.npix2nside(reddening.size):
            print("# resampling extinction map...")
            reddening = hp.ud_grade(reddening, sel_input.extinctionmap_res)
        print("# finding pixels for galaxies...")
        pix = hp.ang2pix(
            sel_input.extinctionmap_res,
            np.pi / 2.0 - cat["dec_gal"] * conv,
            cat["ra_gal"] * conv,
        )
        print("# computing extinction...")
        ext_mag = reddening[pix] * extinction_curve(cat[input.redshift_key])

        print("# constructing selection...")
        selection &= cat[input.flux_key] >= my_flux_limit + ext_mag / 2.5

    elif "flux_limit" in sel_input.selection_keys:

        print("# applying flux limit {}...".format(my_flux))

        selection &= cat[input.flux_key] >= sel_input.selection_logflux_limit

    elif "visibilitymask" in sel_input.selection_keys:

        print(
            "# applying visibility mask {}...".format(
                input.outdir + sel_input.selection_VM_fname
            )
        )

        VM = fits.getdata(input.outdir + sel_input.selection_VM_fname)["VM"]
        VM = hp.ud_grade(VM, sel_input.selection_VM_res)

        conv = np.pi / 180.0
        pix = hp.ang2pix(
            sel_input.selection_VM_res,
            np.pi / 2.0 - cat["dec_gal"] * conv,
            cat["ra_gal"] * conv,
        )

        selection &= cat[input.flux_key] >= VM[pix]

    # selection of centrals and satellites applies only to data catalogs
    if use_data:
        if "central" in sel_input.selection_keys:

            print("# applying selection of centrals...")

            selection &= cat["kind"] == 0

        elif "satellite" in sel_input.selection_keys:

            print("# applying selection of satellites...")

            selection &= cat["kind"] == 1

    if use_data:
        fname = input.selection_data_fname()
    else:
        fname = input.selection_random_fname()

    print(f"# Writing file {fname}...")

    with DefaultCatalogWrite(fname) as store:

        tofits = store.new_array("SELECTION", (Ngal,), dtype=[("SELECTION", bool)])
        tofits["SELECTION"] = selection

    print("# Done!")
