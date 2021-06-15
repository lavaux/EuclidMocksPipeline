
_toolbox = []

def register_tool(function):
    _toolbox.append(function)
    return function

def get_tools():
    return _toolbox