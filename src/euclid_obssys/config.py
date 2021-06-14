from typing import Dict

def readConfig(fname: str) -> Dict:
    """Read a pipeline configuration file.

    That configuration file is written in python

    Args:
        fname (str): File path of the configuration file

    Returns:
        dict: A dictionnary with a number of definitions and functions
    """
    with open(fname, mode="rt") as f:
        code='\n'.join(f.read().splitlines())
    global_dict = {'__builtins__':__builtins__}
    exec(code, global_dict)
    return global_dict
