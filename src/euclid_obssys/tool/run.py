from . import get_tools
from .input_generator import generateConfig
from .footprint import applyFootprintToMaster
from .extract import extractGalaxyCatalogFromMaster
from .hod import createSmoothHOD, createSDHOD_Catalog


class EuclidTool(object):
    """Euclid Observational systematics toolbox

    This is the euclid toolbox for generating mock catalogs from master catalogs, general
    galaxy high level models, and Euclid footprints.

    Commands should be executed in that order:
        * config
        *

    """

    def __init__(self):
        for tool in get_tools():
            setattr(self, tool.__name__, tool)


def run_tool():
    import fire

    fire.Fire(EuclidTool, name="euclid_obssys.tool")
