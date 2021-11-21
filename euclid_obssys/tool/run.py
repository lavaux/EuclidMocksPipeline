from . import get_tools
from . import (
    input_generator,
    footprint,
    extract,
    hod,
    random,
    dn_dz,
    numbercounts,
    create_footprint,
    create_fulloctant_footprint,
    indices_for_sats,
    apply_footprint,
    create_selection,
    writeCatalogs4LE3,
    setup_repo
)


class EuclidTool(object):
    """Euclid Observational systematics toolbox

    This is the euclid toolbox for generating mock catalogs from master catalogs, general
    galaxy high level models, and Euclid footprints.

    Commands should be executed in that order:
        * generateConfig
        * applyFootprintToMaster
        * extractGalaxyCatalogFromMaster
        * createSmoothHOD
        * createSDHOD_Catalog
        * createRandom

    """

    def __init__(self):
        for tool in get_tools():
            setattr(self, tool.__name__, tool)


def run_tool():
    import fire

    fire.Fire(EuclidTool, name="euclid_obssys.tool")
