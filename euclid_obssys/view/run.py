from . import get_tools
from . import angular_map, visualizeHOD, dn_dz


class EuclidViewTool(object):
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

    result = fire.Fire(EuclidViewTool, name="euclid_obssys.view")
    print(f"# Got result: {'OK' if result is None else result}")
    raise SystemExit(0 if result is None else 1)

