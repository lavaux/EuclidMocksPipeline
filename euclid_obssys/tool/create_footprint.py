################################################################################
### Authors: Tiago Castro, Pierluigi Monaco                                  ###
###                                                                          ###
################################################################################
from . import register_tool
from ..config import readConfig


@register_tool
def createFootprint(outdir: str = "Products") -> None:
    import numpy as np
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from euclid_obssys.disk import DefaultCatalogWrite
    import healpy as hp
    from healpy.rotator import Rotator
    import sys
    from os import path

    # these are the parameters of the footprint, please copy them to the input file

    footprint_fname = path.join(outdir, "Footprints", "100sqdeg.fits")
    footprint_res = 2048
    footprint_tag = "100sqdeg"
    footprint_zrange = [0.8, 2.0]

    # linear size of the field in radians
    size = 10.0 * np.pi / 180.0

    # reddening map for checking where the footprint has been placed
    # reddening = hp.ud_grade(hp.read_map(outdir+'ExtinctionMaps/HFI_CompMap_ThermalDustModel_2048_R1.20.fits', field=2), footprint_res)
    # NOTE: This is the new name for this map
    reddening = hp.ud_grade(
        hp.read_map(
            path.join(
                outdir,
                "ExtinctionMaps/COM_CompMap_ThermalDust-commander_2048_R2.00.fits",
            ),
            field=2,
        ),
        footprint_res,
    )

    # creating boolean healpy map
    (theta, phi) = hp.pix2ang(footprint_res, np.arange(hp.nside2npix(footprint_res)))
    footprint = np.zeros(hp.nside2npix(footprint_res), dtype=bool)

    # square of 10 deg of size on the equator
    footprint[
        (theta >= np.pi / 2.0 - size / 2.0)
        & (theta < np.pi / 2.0 + size / 2.0)
        & (phi < size)
    ] = True

    # rotation to a convenient place in galactic coordinates
    rot = Rotator(rot=[-40.0, -40.0, 0.0])
    footprint = rot.rotate_map_pixel(footprint) > 0.5

    sky_fraction = footprint.sum() / footprint.size

    # plot of the location
    foot2 = reddening.copy()
    foot2[footprint] *= 2
    hp.mollview(foot2, max=1000)
    plt.savefig(path.join(outdir, "Plots", "100sqdeg.png"))

    # writes footprint on fits file
    print("## writing footprint on file {}".format(footprint_fname))

    with DefaultCatalogWrite(footprint_fname) as out_file:

        footprint = footprint.astype([("FOOTPRINT_G", int)])
        out_file.set_array("footprint", footprint)
        out_file.add_tag("RES", footprint_res)
        out_file.add_tag("MINZ", footprint_zrange[0])
        out_file.add_tag("MAXZ", footprint_zrange[1])
        out_file.add_tag("TAG", footprint_tag)
        out_file.add_tag("SKYFRAC", sky_fraction)

    print("DONE!!!")
