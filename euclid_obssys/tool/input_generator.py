from . import register_tool


@register_tool
def generateConfig(
    o: str,
    repodir: str,
    projectdir: str,
    lf_model: int = 1,
    rsd: bool = False,
    type: str = "sdhod",
    footprintTag: str = None,
) -> None:
    """This generate a configuration file for the Euclid Observational systematics pipeline.

    Args:
        repodir (str): Repository directory of Euclid data. 
        projectdir (str): Name of the directory holding Project data.
        lf_model (int, optional): Luminosity function model to be adopted.
        rsd (bool, optional): Whether redshift space distortions must be applied. Defaults to False.
        o (str): [description]. Output configuration file Defaults to None.
        footprintTag (str): footprint to use , e.g. "100sqdeg"
        type (str, optional): Type of mock catalog to generate. Must be one of
            "sdhod", "flagship", "pinocchio", "box".

    Raises:
        ValueError: [description]
        ValueError: [description]
    """
    import pkg_resources
    import re
    import os
    import sys

    input_template = pkg_resources.resource_string(
        "euclid_obssys", "templates/input.py"
    ).decode("utf-8")

    if not type in ["sdhod", "flagship", "pinocchio", "box"]:
        raise ValueError("Invalid catalog type")

    if repodir is None:
        repodir = os.getcwd()

    if o is None:
        raise ValueError("Missing output file")

    replacements = {
        "REPODIR": repodir,
        "PROJECTDIR": projectdir,
        "FOOTTAG": f'"{footprintTag}"' if not footprintTag is None else "None",
        "LFMODEL": repr(lf_model),
        "CATTYPE": type,
        "SHUFFLE": "False",
        "SELDATA": "None",
        "SELRAND": "None",
        "RSDFLAG": repr(rsd),
    }

    print("Output directory is {}".format(outdir))

    for k, v in replacements.items():
        input_template = re.sub(k, v, input_template)

    with open(o, mode="wt") as f:
        f.write(input_template)

@register_tool
def generateSelConfig(
    o: str,
) -> None:
    """This generate a configuration file for the Euclid Observational systematics pipeline (selection function).

    Args:
        o (str): [description]. Output configuration file Defaults to None.

    Raises:
        ValueError: [description]
    """
    import pkg_resources
    import re
    import os
    import sys

    input_template = pkg_resources.resource_string(
        "euclid_obssys", "templates/sel_input.py"
    ).decode("utf-8")

    if o is None:
        raise ValueError("Missing output file")

    #replacements = {
    #    "OUTDIR": outdir,
    #    "FOOTTAG": f'"{footprintTag}"' if not footprintTag is None else "None",
    #    "LFMODEL": repr(lf_model),
    #    "CATTYPE": type,
    #    "SHUFFLE": "False",
    #    "SELDATA": "None",
    #    "SELRAND": "None",
    #    "RSDFLAG": repr(rsd),
    #}

    #for k, v in replacements.items():
    #    input_template = re.sub(k, v, input_template)

    with open(o, mode="wt") as f:
        f.write(input_template)

__all__ = ["generateConfig", "generateSelConfig"]
