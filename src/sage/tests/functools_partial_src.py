r"""
Ensure that ``functools.partial`` is correctly handled by
:func:`~sage.misc.sageinspect.sage_getsourcelines`.
"""
from functools import partial

def base(x):
    """
    Test function to make sure
    :func:`~sage.misc.sageinspect.sage_getsourcelines` can get
    the code of a function created by ``functools.partial``.

    EXAMPLES::

        sage: from sage.tests.functools_partial_src import base, _test_func
        sage: base(3)
        21
        sage: _test_func()
        42
    """
    x = x * 7
    return x

_test_func = partial(base, 6)

