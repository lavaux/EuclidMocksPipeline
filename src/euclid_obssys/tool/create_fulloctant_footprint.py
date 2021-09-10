################################################################################
### Authors: Tiago Castro, Pierluigi Monaco                                  ###
###                                                                          ###
################################################################################
from . import register_tool
from ..config import readConfig


@register_tool
def createFullOctantFootprint(outdir: str = "Products") -> None:
    import numpy as np
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from astropy.io import fits
    import healpy as hp
    from healpy.rotator import Rotator
    import sys
    from os import path
    from euclid_obssys.disk import DefaultCatalogWrite

    # these are the parameters of the footprint

    footprint_fname = path.join(outdir, "Footprints", "FullOctant.fits")
    footprint_res = 2048
    footprint_tag = None
    footprint_zrange = [0.8, 2.0]

    # reddening map for checking where the footprint has been placed
    # reddening = hp.ud_grade(hp.read_map(outdir+'ExtinctionMaps/HFI_CompMap_ThermalDustModel_2048_R1.20.fits', field=2), footprint_res)
    reddening = hp.ud_grade(
        hp.read_map(
            path.join(
                outdir,
                "ExtinctionMaps",
                "COM_CompMap_ThermalDust-commander_2048_R2.00.fits",
            ),
            field=2,
        ),
        footprint_res,
    )

    # creating boolean healpy map
    (theta, phi) = hp.pix2ang(footprint_res, np.arange(hp.nside2npix(footprint_res)))
    footprint = np.zeros(hp.nside2npix(footprint_res), dtype=bool)

    # square of 10 deg of size on the equator
    footprint[(theta < np.pi / 2.0) & (phi < np.pi / 2.0)] = True
    sky_fraction = footprint.sum() / footprint.size

    # plot of the location
    foot2 = reddening.copy()
    foot2[footprint] *= 2
    # hp.mollview(footprint.astype(int)+reddening,max=100)
    hp.mollview(foot2, max=1000)
    plt.savefig(path.join(outdir, "Plots", "FullOctant.png"))

    # writes footprint on fits file
    print("## writing footprint on file {}".format(footprint_fname))

    with DefaultCatalogWrite(footprint_fname) as out_file:

        footprint = footprint.astype([("FOOTPRINT_G", long)])
        out_file.set_array("footprint", footprint)
        out_file.add_tag("RES", footprint_res)
        out_file.add_tag("MINZ", footprint_zrange[0])
        out_file.add_tag("MAXZ", footprint_zrange[1])
        out_file.add_tag("TAG", footprint_tag)
        out_file.add_tag("SKYFRAC", sky_fraction)
