from . import register_tool

@register_tool
def generateConfig(
    o: str,
    lf_model: int = 1,
    rsd: bool = False,
    outdir: str = None,
    type: str = "sdhod",
) -> None:
    """This generate a configuration file for the Euclid Observational systematics pipeline.

    Args:
        lf_model (int, optional): Luminosity function model to be adopted.
        rsd (bool, optional): Whether redshift space distortions must be applied. Defaults to False.
        outdir (str, optional): Output directory. Defaults to the current working directory.
        o (str): [description]. Defaults to None.
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

    if outdir is None:
        outdir = os.getcwd()

    if o is None:
        raise ValueError("Missing output file")

    replacements = {
        "OUTDIR": outdir,
        "FOOTTAG": "None",
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


__all__ = ["generateConfig"]
