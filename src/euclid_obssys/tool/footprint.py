################################################################################
### Author: Tiago Castro, Pierluigi Monaco, Guilhem Lavaux                   ###
###                                                                          ###
################################################################################
from . import register_tool
import sys
from ..config import readConfig


@register_tool
def applyFootprintToMaster(config: str):
    """Apply an indicated footprint in the configuration to a large master catalog.

    It applies a footprint to the master catalog, extracting another (smaller) master
    catalog.

    Args:
        config (str): Pipeline config file
    """
    import numpy as np
    import healpy as hp
    from euclid_obssys.disk import DefaultCatalogRead, DefaultCatalogWrite

    input = readConfig(config)

    print(f"# Running applyFootprintToMasterCatalog.py with {config}")
    print("# loading catalog...")

    # input raw catalog
    cat = DefaultCatalogRead(input.build_fname("RawCatalogs", [input.query, None]))[
        "catalog"
    ]

    # loads the survey footprint in equatorial coordinates
    footprint_res, footprint_zrange, sky_fraction, footprint = input.read_footprint()
    print(f"# this footprint covers {sky_fraction*100.}% of the sky")
    print("# selecting galaxies...")

    zused = cat["true_redshift_gal"]
    redshift_sel = (zused >= footprint_zrange[0]) & (zused <= footprint_zrange[1])
    ra_gal = cat["ra_gal"][redshift_sel]
    dec_gal = cat["dec_gal"][redshift_sel]

    conv = np.pi / 180.0
    # if input.rgal is not None:
    #     print("# rotating catalog...")
    #     theta_eq, phi_eq = input.rgal( np.pi/2 - dec_gal * conv, ra_gal * conv )
    # else:
    theta_eq = np.pi / 2.0 - dec_gal * conv
    phi_eq = ra_gal * conv

    ## Get galaxy pixels in the sky
    print("# finding sky pixels...")

    pix = hp.ang2pix(footprint_res, theta_eq, phi_eq)
    foot_sel = np.zeros_like(zused, dtype=bool)
    fp_small = footprint[pix]
    foot_sel[redshift_sel] = fp_small

    Nextract = foot_sel.sum()

    print(f"# Nextract={Nextract} to {input.master_fname()}")

    with DefaultCatalogWrite(input.master_fname()) as store:
        extract = store.add_array("catalog", shape=(Nextract,), dtype=cat.dtype)

        for field in cat.dtype.names:
            print(f"# Doing field {field}")
            extract[field] = cat[field][foot_sel]

    print("# done!")
