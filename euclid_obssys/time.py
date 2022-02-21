import contextlib
import functools
import time


current_contexts = []


@contextlib.contextmanager
def check_time(blockname):
    current_contexts.append(blockname)
    _t0 = time.perf_counter()
    yield
    _t1 = time.perf_counter()
    print(f"TIME: {'/'.join(current_contexts)} has needed {(_t1-_t0):.4f} seconds")
    current_contexts.pop()


def check_time_func(f):
    def new_f(*args, **kwargs):
        with check_time(f.__qualname__):
            return f(*args, **kwargs)

    functools.update_wrapper(new_f, f)
    return new_f


__all__ = ["check_time", "check_time_func"]
