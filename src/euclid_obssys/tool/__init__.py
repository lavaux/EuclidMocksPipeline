from .input_generator import generateConfig
from .footprint import applyFootprintToMaster
from .extract import extractGalaxyCatalogFromMaster

class EuclidTool(object):
    """Euclid Observational systematics toolbox

    This is the euclid toolbox for generating mock catalogs from master catalogs, general
    galaxy high level models, and Euclid footprints.

    Commands should be executed in that order:
        * config
        * 

    """

    def __init__(self):
        self.config = generateConfig
        self.applyFootprintToMaster = applyFootprintToMaster
        self.extractGalaxyCatalogFromMaster = extractGalaxyCatalogFromMaster