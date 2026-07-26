# perf_utils.py
import time
from functools import wraps

def timed(func):
    """Decorator to print execution time of a function."""
    @wraps(func)
    def wrapper(*args, **kwargs):
        start = time.perf_counter()
        result = func(*args, **kwargs)
        elapsed = time.perf_counter() - start
        print(f"⏱️ {func.__name__} took {elapsed:.5f} seconds")
        return result
    return wrapper