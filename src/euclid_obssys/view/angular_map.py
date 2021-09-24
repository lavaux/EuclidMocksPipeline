################################################################################
### Authors: Tiago Castro, Pierluigi Monaco, Guilhem Lavaux                  ###
###                                                                          ###
################################################################################
from . import register_view_tool


@register_view_tool
def angular_map(
    catalog: str,
    selection: str,
    z1: float,
    z2: float,
    Nside: int = 256,
    output: str = None,
) -> None:
    """[summary]

    Args:
        catalog (str): [description]
        selection (str): [description]
        z1 (float): [description]
        z2 (float): [description]
        Nside (int, optional): [description]. Defaults to 256.
    """

    import numpy as np
    import matplotlib.pyplot as plt
    import healpy as hp
    from astropy.io import fits
    import sys
    from euclid_obssys.disk import DefaultCatalogRead

    print("# Running angular_map.py")

    Nside = 256

    fname = catalog
    if selection == "None":
        selection = None

    print("Plotting map of catalog {} in z=[{},{}]".format(fname, z1, z2))

    print("Reading catalog...")
    with DefaultCatalogRead(fname) as store:
        cat = store["catalog"]

    if selection is not None:
        print(f"Reading selection {selection}...")
        with DefaultCatalogRead(selection) as store:
            mysel = store['SELECTION']['SELECTION']
    else:
        mysel = np.ones(len(cat), dtype=bool)

    print("Filtering catalog...")
    filter = (cat["true_redshift_gal"] >= z1) & (cat["true_redshift_gal"] < z2) & mysel

    print("Finding pixels for {} galaxies...".format(filter.sum()))
    conv = np.pi / 180.0
    pix = hp.ang2pix(
        Nside, np.pi / 2.0 - cat["dec_gal"][filter] * conv, cat["ra_gal"][filter] * conv
    )

    print("Constructing map...")
    map = np.zeros((hp.nside2npix(Nside)), dtype=int)
    for p in pix:
        map[p] += 1

    fig = plt.figure(num=1)
    hp.mollview(map,rot=[0,0,0],title=fname,fig=1)
    if output is not None:
        fig.savefig(output)
    else:
        plt.show()

    print("Done!")
