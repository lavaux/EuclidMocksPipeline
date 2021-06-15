from .input_generator import generateConfig
from .footprint import applyFootprintToMaster


class EuclidTool(object):
    def __init__(self):
        self.config = generateConfig
        self.applyFootprintToMaster = applyFootprintToMaster
