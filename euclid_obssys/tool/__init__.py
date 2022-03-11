import inspect
import traceback

_toolbox = []


def register_tool(function):
    def new_f(*args, **kwargs):
        try:
            function(*args, **kwargs)
        except Exception as e:
            print(f"Error while running {function.__name__}.")
            print(f"Exception was {e}")
            traceback.print_exc()
            return e

    new_f.__doc__ = function.__doc__
    new_f.__signature__ = inspect.signature(function)
    new_f.__name__ = function.__name__

    _toolbox.append(new_f)
    return new_f


def get_tools():
    return _toolbox
