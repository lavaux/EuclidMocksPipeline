from .. import tool
from ..tool import run as tool_run


_module_symbols = globals()
__all__ = []


# https://stackoverflow.com/questions/15506971/recursive-version-of-reload
def _rreload(
    module,
    paths=None,
    mdict=None,
    base_module=None,
    blacklist=None,
    reloaded_modules=None,
):
    """Recursively reload modules."""
    from importlib import reload
    from types import ModuleType
    import sys
    import os

    if paths is None:
        paths = [""]
    if mdict is None:
        mdict = {}
    if module not in mdict:
        # modules reloaded from this module
        mdict[module] = []
    if base_module is None:
        base_module = module
    if blacklist is None:
        blacklist = ["importlib", "typing", "_frozen_importlib"]
    if reloaded_modules is None:
        reloaded_modules = []
    if not module.__name__.startswith("euclid_obssys"):
       return reloaded_modules
    try:
        reload(module)
    except NotImplementedError:
        pass
    reloaded_modules.append(module.__name__)
    for attribute_name in dir(module):
        attribute = getattr(module, attribute_name)
        if type(attribute) is ModuleType and attribute.__name__ not in blacklist:
            if attribute not in mdict[module]:
                if attribute.__name__ not in sys.builtin_module_names:
                    # if os.path.dirname(attribute.__file__) in paths:
                    mdict[module].append(attribute)
                    reloaded_modules = _rreload(
                        attribute,
                        paths,
                        mdict,
                        base_module,
                        blacklist,
                        reloaded_modules,
                    )
    try:
        reload(module)
    except NotImplementedError:
        pass
    return reloaded_modules


def _fillup_symbols():
    global __all__, _module_symbols
    _module_symbols.update({f.__name__: f for f in tool.get_tools()})
    __all__.clear()
    __all__.append("refresh")
    __all__.extend([f.__name__ for f in tool.get_tools()])


def refresh():
    from importlib import reload

    global tool, tool_run

    tool.clear_tools()
    old_dask = tool._dask_started

    tool = reload(tool)
    tool._dask_started = old_dask
    print(_rreload(tool_run))
    _fillup_symbols()


_fillup_symbols()
