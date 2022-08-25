################################################################################
### Authors: Tiago Castro, Pierluigi Monaco                                  ###
###                                                                          ###
################################################################################
from . import register_tool
from ..config import readConfig


@register_tool
def createFootprint(outdir: str = "Repo", tag: str = "100sqdeg", size: float = 10.0) -> None:
    """Create a survey footprint on the sky, starting from the north pole.
    
    Args:
        * outdir (str): output repository directory for the footprint
        * tag (str): Tag for the footprint
        * size (float): Size in degrees of the footprint
    """
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

    footprint_fname = path.join(outdir, "Footprints", f"{tag}.fits")
    footprint_res = 2048
    footprint_tag = tag
    footprint_zrange = [0.8, 2.0]

    # linear size of the field in radians
    size *= np.pi / 180.0

    # reddening map for checking where the footprint has been placed
    # reddening = hp.ud_grade(hp.read_map(outdir+'ExtinctionMaps/HFI_CompMap_ThermalDustModel_2048_R1.20.fits', field=2), footprint_res)
    # NOTE: This is the new name for this map
    reddening = hp.ud_grade(
        hp.read_map(
            path.join(
                outdir,
                "ExtinctionMaps/HFI_CompMap_ThermalDustModel_2048_R1.20.fits",
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
    maxval = foot2[footprint].max()
    foot2[footprint] =1
    foot2[footprint==False] = 0#foot2[footprint].min()
    hp.mollview(foot2, max=1)
    plt.savefig(path.join(outdir, "Footprints", f"{tag}.png"))

    # writes footprint on fits file
    print("## writing footprint on file {}".format(footprint_fname))

    with DefaultCatalogWrite(footprint_fname) as out_file:

        footprint = footprint.astype([("FOOTPRINT_G", int)])
        out_file.set_array(
            "footprint",
            footprint,
            tags={
                "RES": footprint_res,
                "MINZ": footprint_zrange[0],
                "MAXZ": footprint_zrange[1],
                "TAG": footprint_tag,
                "SKYFRAC": sky_fraction,
            },
        )

    print("DONE!!!")
