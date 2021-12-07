from typing import Tuple
import numpy as np
from .disk import ensurePath

"""
This constructs the file names for all files accessible by the pipeline
"""

# critical point: there is an inconsistency in the use of smoothing length


def tags(pars: object, runstart: int, runstop: int, redshift: Tuple[float, float]):
    """this returns a dictionary with all the possible tags to be used

    Args:
        pars ([type]): [description]
        runstart (int): [description]
        runstop (int): [description]
        redshift (Tuple[float,float]): [description]

    Returns:
        dict: all possible tags t be used
    """
    taglist = {
        "random": "random",
        "sel_data": pars.selection_data_tag,  # selection tag for the data catalog
        "sel_rand": pars.selection_random_tag,  # selection tag for the random catalog
        "sel_rand2": None,  # here random is constructed from selected data
        "sel_both": "ds{}_rs{}".format(
            pars.selection_data_tag, pars.selection_random_tag
        ),  # combined selection tag
        "cat_type": pars.cat_type,
        "query": pars.query,
        "footp": pars.footprint_tag,
        "other": pars.otherlab,
        "lf_model": "m{}".format(pars.lf_model),
        "sm": "smL{}".format(pars.smoothing_length),
        "cm": "cM{}".format(pars.cmrelation),
        "random_sm": "nosm",
        "shuffled": None,
        "rsd": "truez",
        "alpha": "{}x".format(pars.alpha),
        "ngrid": "gr{}".format(pars.ngrid),
        "Lbox": "box{}".format(pars.Lbox),
        "redshift": None,
    }

    # fixes specific tags
    if pars.random_realization_number is not None:
        taglist["random"] = "random{}".format(pars.random_realization_number)

    if pars.smooth_dndz_in_random:
        taglist["random_sm"] = taglist["sm"]

    if pars.RSD:
        taglist["rsd"] = "obsz"

    if pars.shuffled_fluxes:
        taglist["shuffled"] = "shuffled"

    if pars.cat_type is "pinocchio":
        if runstart is None:
            taglist["cat_type"] = "{}".format(pars.pinocchio_kernel)
        else:
            if runstop is None:
                taglist["cat_type"] = "{}{:04d}".format(pars.pinocchio_kernel, runstart)
            else:
                taglist["cat_type"] = "{}{:04d}-{:04d}".format(
                    pars.pinocchio_kernel, runstart, runstop
                )

    # the random can be constructed on the selected data, this tag highlights this choice
    if pars.apply_dataselection_to_random:
        taglist["sel_rand2"] = pars.selection_random_tag

    if redshift is not None:
        taglist["redshift"] = "z{}-{}".format(redshift[0], redshift[1])

    return taglist


def build_fname(
    pars: str,
    thisdir: str,
    tag_filter,
    ext: str = ".fits",
    head: str = None,
    tail: str = None,
    RepoDirectory: bool = False,
    runstart: int = None,
    runstop: int = None,
    redshift: Tuple[float, float] = None,
    skip_tags: bool = True,
    skip_project_dir: bool = False,
):
    """this constructs a file name by adding together a number of tags
    separated by a _ if they are not None, otherwise they are skipped

    Args:
        pars ([type]): the input parameter structure
        thisdir ([type]): the directory where to put the file
        tag_filter ([type]): selects the tags to be written in the file
        ext (str, optional): file extension. Defaults to ".fits".
        head ([type], optional): a tag to be added at the beginning of the filename. Defaults to None.
        tail ([type], optional): a tag to be added before the file extension. Defaults to None.
        RepoDirectory (bool, optional): if True the file is in the repository otherwise in the project. Defaults to False.
        runstart ([type], optional): gives the pinocchio run or its range. Defaults to None.
        runstop ([type], optional): gives the pinocchio run or its range. Defaults to None.
        redshift ([type], optional):  gives the redshift bin of the catalog slice. Defaults to None.
        skip_tags (bool, optional): skip tags contained in pars.skip_tags. Defaults to True.
        skip_project_dir (bool, optional): skip the project or repo directory. Defaults to False.

    Returns:
        str: the filename
    """

    # if needed, adds a / to the requested directory file
    if thisdir[-1] is not "/":
        thisdir += "/"
    # if needed, puts '.' at the beginning of ext
    if ext is not None and ext[0] is not ".":
        ext = "." + ext

    # the root directory is either the repository or the project directory
    if not skip_project_dir:
        if RepoDirectory:
            fname = pars.repo + thisdir
        else:
            fname = pars.project + thisdir
    else:
        fname = thisdir

    ensurePath(fname)

    # list of tags to skip
    if skip_tags and pars.skip_tags is not None:
        exclude = pars.skip_tags
    else:
        exclude = []

    mytags = tags(pars, runstart, runstop, redshift)

    first = True
    if head is not None:
        fname += head
        first = False
    for element in mytags:
        tag = mytags[element]
        if element in tag_filter and element not in exclude and tag is not None:
            if not first:
                if tag[0] is not "_":
                    tag = "_" + tag
            first = False
            fname += tag
    if tail is not None:
        fname += "_" + tail
    if ext is not None:
        fname += ext
    return fname


def exclude_dir(fname):
    slashes = [pos for pos, char in enumerate(fname) if char == "/"]
    return fname[slashes[-1] + 1 :]


# file names

# these are files found in the general repo
# for files in the repo tag skipping is suppressed
def master(pars):
    """
    file name of the master catalog
    """
    tag_filter = ["query", "footp"]
    return build_fname(
        pars, "RawCatalogs", tag_filter, RepoDirectory=True, skip_tags=False
    )


def indices(pars):
    """
    file name for the indices of a master catalog
    (NB: it accepts a footprint, but it is very likely that they are produced only for the octant)
    """
    tag_filter = ["query", "footp"]
    return build_fname(
        pars,
        "RawCatalogs",
        tag_filter,
        RepoDirectory=True,
        tail="indices",
        skip_tags=False,
    )


def SDHOD(pars):
    """
    name of SDHOD tables file
    """
    tag_filter = ["query", "lf_model", "sm", "cm"]
    return build_fname(
        pars, "SDHOD", tag_filter, RepoDirectory=True, head="SDHOD", skip_tags=False
    )


# the filenames of pinocchio runs is outside the tagging system
def pinplc(pars, myrun):
    run = "{}{:04d}".format(pars.pinocchio_kernel, myrun)
    return "{}/{}/pinocchio.{}.plc.out".format(pars.pinocchio_repo, run, run)


# these are files produced by the pipeline and written in the project directory
def flagcat(pars):
    """
    name of flagship galaxy catalog, bypasses cat_type
    """
    save = pars.cat_type
    pars.cat_type = "flagship"
    tag_filter = ["cat_type", "query", "footp", "lf_model"]
    fname = build_fname(pars, "GalaxyCatalogs", tag_filter)
    pars.cat_type = save
    return fname


def sel_input(pars, use_data: bool) -> str:
    if use_data:
       if pars.selection_data_tag is None:
            print("No selection specified for galaxy catalog, exiting")
            sys.exit(0)
       else:
            sel_input_tag = "sel_data"
    else:
        if pars.selection_random_tag is None:
            print("No selection specified for random catalog, exiting")
            sys.exit(0)
        else:
            sel_input_tag = "sel_rand"


    tag_filter = [sel_input_tag]
    sel_input_fname = build_fname(
        pars, "SelectionInputs", tag_filter, head="sel_input", RepoDirectory=True, ext=".py"
    )
    return sel_input_fname


def hodcat(pars):
    """
    name of smooth hod galaxy catalog, bypasses cat_type
    """
    save = pars.cat_type
    pars.cat_type = "sdhod"  # hodcat"
    tag_filter = [
        "cat_type",
        "query",
        "footp",
        "other",
        "lf_model",
        "cm",
    ]  # we should add here sm!
    fname = build_fname(pars, "GalaxyCatalogs", tag_filter)
    pars.cat_type = save
    return fname


def pincat(pars, myrun):
    """
    name of pinocchio galaxy catalog, bypasses cat_type
    """
    save = pars.cat_type
    pars.cat_type = "pinocchio"
    tag_filter = ["cat_type", "footp", "other", "lf_model", "cm"]
    fname = build_fname(pars, "GalaxyCatalogs", tag_filter, runstart=myrun)
    pars.cat_type = save
    return fname


def galcat_kernel(pars):
    """
    kernel for galaxy catalog name
    """
    if pars.cat_type is "flagship":
        tag_filter = ["cat_type", "query", "footp", "lf_model"]
    elif pars.cat_type is "sdhod":
        tag_filter = ["cat_type", "query", "footp", "other", "lf_model", "cm"]
    elif pars.cat_type is "pinocchio":
        tag_filter = ["cat_type", "footp", "other", "lf_model", "cm"]
    else:
        print("ERROR: unrecognised cat_type in configuration file")
        return None
    return tag_filter


def galcat(pars, myrun=None):
    """
    name of galaxy catalog
    """
    tag_filter = galcat_kernel(pars)
    return build_fname(pars, "GalaxyCatalogs", tag_filter, runstart=myrun)


def random(pars):
    """
    name of random catalog
    """
    tag_filter = galcat_kernel(pars)
    tag_filter.append("random")
    tag_filter.append("sel_rand2")
    tag_filter.append("random_sm")
    tag_filter.append("rsd")
    tag_filter.append("alpha")
    return build_fname(
        pars,
        "RandomCatalogs",
        tag_filter,
        runstart=pars.pinocchio_first_run,
        runstop=pars.pinocchio_last_run,
    )


def selection_data(pars, myrun=None):
    """
    name of selection file for data
    """
    if pars.selection_data_tag is None:
        return ""
    else:
        tag_filter = galcat_kernel(pars)
        tag_filter.append("sel_data")
        tag_filter.append("shuffled")
        return build_fname(pars, "Selections", tag_filter, runstart=myrun, head="data")


def selection_random(pars):
    """
    name of selection file for data
    """
    if pars.selection_random_tag is None:
        return ""
    else:
        tag_filter = galcat_kernel(pars)
        tag_filter.append("sel_rand")
        tag_filter.append("shuffled")
        return build_fname(
            pars,
            "Selections",
            tag_filter,
            runstart=pars.pinocchio_first_run,
            runstop=pars.pinocchio_last_run,
            head="random",
        )


def dndz(pars):
    """
    name for dndz
    """
    tag_filter = galcat_kernel(pars)
    tag_filter.append("sel_data")
    tag_filter.append("shuffled")
    tag_filter.append("sm")
    tag_filter.append("rsd")
    return build_fname(
        pars,
        "NumberCounts",
        tag_filter,
        runstart=pars.pinocchio_first_run,
        runstop=pars.pinocchio_last_run,
        head="dndz",
    )


def LE3_data(pars, z1, z2, myrun=None):
    """
    name of data catalog for LE3
    """
    tag_filter = galcat_kernel(pars)
    tag_filter.append("sel_data")
    tag_filter.append("shuffled")
    tag_filter.append("rsd")
    tag_filter.append("redshift")

    return build_fname(
        pars, "Catalogs4LE3", tag_filter, runstart=myrun, head="data", redshift=[z1, z2]
    )


def LE3_random(pars, z1, z2):
    """
    name of random catalog for LE3
    """
    tag_filter = galcat_kernel(pars)
    tag_filter.append("random")
    tag_filter.append("sel_rand")
    tag_filter.append("random_sm")
    tag_filter.append("rsd")
    tag_filter.append("alpha")
    tag_filter.append("redshift")

    return build_fname(
        pars,
        "Catalogs4LE3",
        tag_filter,
        runstart=pars.pinocchio_first_run,
        runstop=pars.pinocchio_last_run,
        redshift=[z1, z2],
    )


def numbercounts(pars, z1, z2, myrun=None):
    """
    name for number counts
    """
    tag_filter = galcat_kernel(pars)
    tag_filter.append("sel_data")
    tag_filter.append("shuffled")
    tag_filter.append("rsd")
    tag_filter.append("redshift")
    return build_fname(
        pars,
        "NumberCounts",
        tag_filter,
        runstart=myrun,
        head="numbercounts",
        redshift=[z1, z2],
    )


def delta_data(pars, z1, z2, myrun=None):
    """
    name of data catalog angular density map
    """
    tag_filter = galcat_kernel(pars)
    tag_filter.append("sel_data")
    tag_filter.append("shuffled")
    tag_filter.append("rsd")
    tag_filter.append("redshift")

    return build_fname(
        pars, "Cls", tag_filter, runstart=myrun, head="delta_data", redshift=[z1, z2]
    )


def delta_random(pars, z1, z2):
    """
    name of random catalog angular density map
    """
    tag_filter = galcat_kernel(pars)
    tag_filter.append("random")
    tag_filter.append("sel_rand")
    tag_filter.append("random_sm")
    tag_filter.append("rsd")
    tag_filter.append("alpha")
    tag_filter.append("redshift")

    return build_fname(
        pars,
        "Cls",
        tag_filter,
        head="delta",
        runstart=pars.pinocchio_first_run,
        runstop=pars.pinocchio_last_run,
        redshift=[z1, z2],
    )


def cls(pars, z1, z2, myrun=None):
    """
    name of Cls file
    """
    tag_filter = galcat_kernel(pars)
    tag_filter.append("sel_data")
    tag_filter.append("shuffled")
    tag_filter.append("rsd")
    tag_filter.append("redshift")

    return build_fname(
        pars, "Cls", tag_filter, runstart=myrun, head="Cls", redshift=[z1, z2]
    )


def estimator_params(pars, estimator, z1, z2, myrun=None):
    """
    name of parameter file for an estimator
    """
    tag_filter = galcat_kernel(pars)
    tag_filter.append("sel_both")
    tag_filter.append("sh")
    if pars.random_realization_number is not None:
        tag_filter.append("random")
    tag_filter.append("random_sm")
    tag_filter.append("rsd")
    tag_filter.append("alpha")
    tag_filter.append("ngrid")
    tag_filter.append("Lbox")
    tag_filter.append("redshift")
    tag_filter.append("shuffled")

    return build_fname(
        pars,
        estimator + "/Params",
        tag_filter,
        head=estimator + "params",
        runstart=myrun,
        redshift=[z1, z2],
        ext=".ini",
    )


def estimator_list(pars, estimator, z1, z2):
    """
    name of json list for CM-*
    """
    tag_filter = galcat_kernel(pars)
    tag_filter.append("sel_both")
    tag_filter.append("sh")
    if pars.random_realization_number is not None:
        tag_filter.append("random")
    tag_filter.append("random_sm")
    tag_filter.append("rsd")
    tag_filter.append("alpha")
    tag_filter.append("ngrid")
    tag_filter.append("Lbox")
    tag_filter.append("redshift")
    tag_filter.append("shuffled")

    return build_fname(
        pars,
        estimator + "/Params",
        tag_filter,
        head=estimator + "list",
        redshift=[z1, z2],
        ext=".json",
    )


def estimator_measure(pars, estimator, z1, z2, myrun=None, skip_project_dir=False):
    """
    name of parameter file for an estimator
    """
    tag_filter = galcat_kernel(pars)
    tag_filter.append("sel_both")
    tag_filter.append("sh")
    if pars.random_realization_number is not None:
        tag_filter.append("random")
    tag_filter.append("random_sm")
    tag_filter.append("rsd")
    tag_filter.append("alpha")
    tag_filter.append("ngrid")
    tag_filter.append("Lbox")
    tag_filter.append("redshift")
    tag_filter.append("shuffled")

    return build_fname(
        pars,
        estimator + "/Measures",
        tag_filter,
        head=estimator,
        runstart=myrun,
        redshift=[z1, z2],
        ext=None,
        skip_project_dir=skip_project_dir,
    )


def estimator_script(pars, estimator, number):
    """
    name of parameter file for an estimator
    """
    tag_filter = galcat_kernel(pars)
    tag_filter.append("sel_both")
    tag_filter.append("sh")
    if pars.random_realization_number is not None:
        tag_filter.append("random")
    tag_filter.append("random_sm")
    tag_filter.append("rsd")
    tag_filter.append("alpha")
    tag_filter.append("ngrid")
    tag_filter.append("Lbox")
    tag_filter.append("shuffled")

    return build_fname(
        pars,
        estimator + "/Scripts",
        tag_filter,
        head=estimator + "script",
        runstart=pars.pinocchio_first_run,
        runstop=pars.pinocchio_last_run,
        tail="{}".format(number),
        ext=".sh",
    )


def plot_estimator(pars, estimator, order, z1, z2, myrun=None):
    """
    name of parameter file for an estimator
    """
    tag_filter = galcat_kernel(pars)
    tag_filter.append("sel_both")
    tag_filter.append("sh")
    if pars.random_realization_number is not None:
        tag_filter.append("random")
    tag_filter.append("random_sm")
    tag_filter.append("rsd")
    tag_filter.append("alpha")
    tag_filter.append("ngrid")
    tag_filter.append("Lbox")
    tag_filter.append("redshift")
    tag_filter.append("shuffled")

    return build_fname(
        pars,
        "Plots",
        tag_filter,
        head=estimator + "{}".format(order),
        runstart=myrun,
        redshift=[z1, z2],
        ext=".png",
    )


def plot_dndz(pars):
    """
    name of Cls file
    """
    tag_filter = galcat_kernel(pars)
    tag_filter.append("sel_both")
    tag_filter.append("sm")
    tag_filter.append("shuffled")
    tag_filter.append("rsd")
    tag_filter.append("redshift")

    return build_fname(
        pars,
        "Plots",
        tag_filter,
        head="dndz",
        ext=".png",
        runstart=pars.pinocchio_first_run,
        runstop=pars.pinocchio_last_run,
    )


def plot_cls(pars, z1, z2, myrun=None):
    """
    name of Cls file
    """
    tag_filter = galcat_kernel(pars)
    tag_filter.append("sel_data")
    tag_filter.append("shuffled")
    tag_filter.append("rsd")
    tag_filter.append("redshift")

    return build_fname(
        pars,
        "Plots",
        tag_filter,
        runstart=myrun,
        head="Cls",
        redshift=[z1, z2],
        ext=".png",
    )


def plot_numbercounts(pars, z1, z2, myrun=None):
    """
    name for number counts
    """
    tag_filter = galcat_kernel(pars)
    tag_filter.append("sel_data")
    tag_filter.append("shuffled")
    tag_filter.append("rsd")
    tag_filter.append("redshift")

    return build_fname(
        pars,
        "Plots",
        tag_filter,
        runstart=myrun,
        head="numbercounts",
        redshift=[z1, z2],
        ext=".png",
    )

def plot_hod(pars, type="centrals", myrun=None):
    tag_filter = galcat_kernel(pars)
    tag_filter.append("sel_data")
    tag_filter.append("shuffled")
    tag_filter.append("rsd")
    tag_filter.append("redshift")

    return build_fname(pars, "Plots", tag_filter, runstart=myrun, head=f"SDHOD_{type}", ext=".png")
