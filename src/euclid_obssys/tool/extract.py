################################################################################
### Author: Tiago Castro, Pierluigi Monaco, Guilhem Lavaux                   ###
###                                                                          ###
################################################################################
from . import register_tool
from euclid_obssys.config import readConfig


@register_tool
def extractGalaxyCatalogFromMaster(config: str):
    """Extract a galaxy catalog from a master catalog based on pipeline config.

    extracts a galaxy
    catalog from a master catalog by applying the standard flux limit,
    according to the model. Run it separately for each model.

    Args:
        config (str): The pipeline configuration file
    """
    from euclid_obssys.disk import DefaultCatalogWrite, DefaultCatalogRead
    import numpy as np

    input = readConfig(config)

    print(f"# Running extractGalaxyCatalogFromMaster.py with {config}")

    footprint_res, footprint_zrange, sky_fraction, footprint = input.read_footprint()
    del footprint

    print(f"# loading catalog {input.master_fname()}...")

    # input raw catalog
    with DefaultCatalogRead(input.master_fname()) as cat_file:
      cat = cat_file["catalog"]

    print("# selecting galaxies...")

    selection = cat[input.flux_key] > input.logflux_limit

    Nextract = selection.sum()

    print(f"# extracting {Nextract} galaxies")

    extract = np.empty(
        Nextract,
        dtype=[
            ("x_gal", float),
            ("y_gal", float),
            ("z_gal", float),
            ("ra_gal", float),
            ("dec_gal", float),
            ("kind", int),
            ("true_redshift_gal", float),
            ("observed_redshift_gal", float),
            ("halo_lm", float),
            ("galaxy_id", int),
            ("halo_id", int),
            (input["flux_key"], float),
            ("sh_" + input["flux_key"], float),
        ],
    )

    for field in extract.dtype.names:
        print("    processing {}".format(field))

        if field == "sh_" + input["flux_key"]:
            print("# shuffling galaxies...")

            shuffled = np.copy(cat[input["flux_key"]][selection])

            zbins = np.linspace(
                footprint_zrange[0],
                footprint_zrange[1],
                round(
                    (footprint_zrange[1] - footprint_zrange[0]) / input["deltazbin"] + 1
                ),
            )
            zindex = np.digitize(extract["true_redshift_gal"], zbins)
            for iz in np.arange(zbins.size - 1):
                ff = zindex == iz
                this = shuffled[ff]
                np.random.shuffle(this)
                shuffled[ff] = this
            extract[field] = shuffled
        else:

            extract[field] = cat[field][selection]

    del cat

    fname = input.flagcat_fname()

    print(f"# writing file {fname}")
    with DefaultCatalogWrite(fname) as out:
        out.set_array("catalog", extract)

    print("# done!")
