import inspect
import traceback
from functools import wraps
import logging

_toolbox = []
_dask_started = False
_debug = False
_log = logging.getLogger("tools")

_log.setLevel(logging.INFO)

def register_tool(function):
    global _debug
    if _debug:
        print(f"Register {function.__name__}")
    def new_f(*args, **kwargs):
        try:
            print(f"# Starting up '{function.__name__}'")
            return function(*args, **kwargs)
        except Exception as e:
            print(f"Error while running {function.__name__}.")
            print(f"Exception was {e}")
            traceback.print_exc()
            return e

    new_f = wraps(function)(new_f)

    _toolbox.append(new_f)
    return new_f

def need_dask(function):
    def new_f(*args, **kwargs):
        from distributed import Client
        global _dask_started
        if not _dask_started:
           _log.info("Starting dask")
           client = Client(n_workers=8)
           _dask_started = True

        return function(*args, **kwargs)

    return wraps(function)(new_f)

def clear_tools():
    _toolbox.clear()

def get_tools():
    return _toolbox
