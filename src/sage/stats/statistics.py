r"""
A sage compatible version of Python statistics module

This module currently implement the following subsets of functions
also present in the Python module statistics

- ``mean()``
- ``stdev()``, ``pstdev()``
- ``variance()``, ``pvariance()``

For convenience, it also imports the following functions directly
from the Python statistics module

- ``mode()``
- ``multimode()``
- ``median()``

TESTS:

Some doctests from the former ``basic_stats`` module::

        sage: import sage.stats.statistics as statistics

        sage: v = [1,2,4,1,6,2,6,7,1]
        sage: statistics.multimode(v)
        [1]
        sage: v.count(1)
        3
        sage: statistics.multimode([])
        []

        sage: statistics.multimode([1,2,3,4,5])
        [1, 2, 3, 4, 5]
        sage: statistics.multimode([3,1,2,1,2,3])
        [3, 1, 2]
        sage: statistics.multimode([0, 2, 7, 7, 13, 20, 2, 13])
        [2, 7, 13]

        sage: statistics.multimode(['sage', 'four', 'I', 'three', 'sage', 'pi'])
        ['sage']

        sage: statistics.median([1,2,3,4,5])
        3
        sage: statistics.median([e, pi])
        1/2*pi + 1/2*e
        sage: statistics.median(['sage', 'linux', 'python'])
        'python'
        sage: statistics.median([])
        Traceback (most recent call last):
        ...
        StatisticsError: no median for empty data
"""

from statistics import mode, multimode, median
from warnings import warn

from sage.structure.coerce import is_numpy_type
from sage.symbolic.constants import NaN
from sage.misc.functional import sqrt

def mean(v):
    """
    Return the mean of the elements of ``v``.

    For empty data ``v``, a warning is raised and the (symbolic) ``NaN`` is
    returned.

    EXAMPLES::

        sage: import sage.stats.statistics as statistics
        sage: statistics.mean([pi, e])
        1/2*pi + 1/2*e
        sage: statistics.mean([I, sqrt(2), 3/5])
        1/3*sqrt(2) + 1/3*I + 1/5
        sage: statistics.mean([RIF(1.0103,1.0103), RIF(2)])
        1.5051500000000000?
        sage: statistics.mean(range(4))
        1.5
        sage: statistics.mean(srange(4))
        3/2

        sage: import numpy
        sage: x = numpy.array([1,2,3,4,5])
        sage: statistics.mean(x)
        3.0

        sage: statistics.mean([])
        doctest:warning
        ...
        RuntimeWarning: mean of empty data
        NaN
    """
    if is_numpy_type(type(v)):
        import numpy
        if isinstance(v, numpy.ndarray):
            return v.mean()

    if not v:
        warn("mean of empty data", RuntimeWarning)
        return NaN

    return sum(v) / len(v)


def stdev(v):
    r"""
    Return the standard deviation of the elements of ``v``.

    For empty data ``v``, a warning is raised and the (symbolic) ``NaN`` is
    returned.

    EXAMPLES::

        sage: import sage.stats.statistics as statistics
        sage: statistics.stdev([1..6])
        sqrt(7/2)
        sage: statistics.stdev([e, pi])
        sqrt(1/2)*abs(pi - e)

        sage: import numpy
        sage: x = numpy.array([1,2,3,4,5])
        sage: statistics.stdev(x)
        1.5811388300841898
        sage: statistics.stdev([1.,2.,3.,4.,5.])
        1.58113883008419

        sage: statistics.stdev([])
        doctest:warning
        ...
        RuntimeWarning: standard deviation of empty data
        NaN
    """
    if is_numpy_type(type(v)):
        import numpy
        if isinstance(v, numpy.ndarray):
            return v.std(ddof=1)

    if not v:
        warn("standard deviation of empty data", RuntimeWarning)
        return NaN

    return sqrt(variance(v))


def pstdev(v):
    """
    Return the population standard deviation of the elements of ``v``.

    For empty data ``v``, a warning is raised and the (symbolic) ``NaN`` is
    returned.

    EXAMPLES::

        sage: import sage.stats.statistics as statistics
        sage: statistics.pstdev([1..6])
        1/2*sqrt(35/3)
        sage: statistics.pstdev([e, pi])
        1/2*abs(pi - e)
        sage: statistics.pstdev([I, sqrt(2), 3/5])
        1/15*sqrt(1/3)*sqrt((10*sqrt(2) - 5*I - 3)^2 +
        (5*sqrt(2) - 10*I + 3)^2 + (5*sqrt(2) + 5*I - 6)^2)
        sage: statistics.pstdev([RIF(1.0103, 1.0103), RIF(2)])
        0.4948500000000000?

        sage: import numpy
        sage: x = numpy.array([1,2,3,4,5])
        sage: statistics.pstdev(x)
        1.4142135623730951
        sage: statistics.pstdev([1.,2.,3.,4.,5.])
        1.41421356237310

        sage: statistics.pstdev([])
        doctest:warning
        ...
        RuntimeWarning: population standard deviation of empty data
        NaN

    TESTS::
    
        sage: import sage.stats.statistics as sage_statistics
        sage: import statistics as python_statistics
        sage: v = [float(1), float(2), float(3)]
        sage: sage_statistics.pstdev(v) == python_statistics.pstdev(v)
        True
    """
    if is_numpy_type(type(v)):
        import numpy
        if isinstance(v, numpy.ndarray):
            return v.std()

    if not v:
        warn("population standard deviation of empty data", RuntimeWarning)
        return NaN

    return sqrt(pvariance(v))


def variance(v):
    r"""
    Return the variance of the elements of `v`.

    For empty data ``v``, a warning is raised and the (symbolic) ``NaN`` is
    returned.

    EXAMPLES::

        sage: import sage.stats.statistics as statistics
        sage: statistics.variance([1..6])
        7/2
        sage: statistics.variance([e, pi])
        1/2*(pi - e)^2
        sage: statistics.variance([I, sqrt(2), 3/5])
        1/450*(10*sqrt(2) - 5*I - 3)^2 + 1/450*(5*sqrt(2) - 10*I + 3)^2
        + 1/450*(5*sqrt(2) + 5*I - 6)^2

        sage: import numpy
        sage: x = numpy.array([1,2,3,4,5])
        sage: statistics.variance(x)
        2.5
        sage: statistics.variance([1.,2.,3.,4.,5.])
        2.50000000000000

        sage: statistics.variance([])
        doctest:warning
        ...
        RuntimeWarning: variance of empty data
        NaN

    TESTS::

        sage: import sage.stats.statistics as sage_statistics
        sage: import statistics as python_statistics
        sage: v = [float(1), float(2), float(3)]
        sage: sage_statistics.variance(v) == python_statistics.variance(v)
        True
    """
    if is_numpy_type(type(v)):
        import numpy
        if isinstance(v, numpy.ndarray):
            return v.var(ddof=1)

    if not v:
        warn("variance of empty data", RuntimeWarning)
        return NaN

    mu = mean(v)
    return sum( (x - mu)**2 for x in v) / (len(v) - 1)


def pvariance(v):
    r"""
    Return the population variance of the elements of `v`.

    For empty data ``v``, a warning is raised and the (symbolic) ``NaN`` is
    returned.

    EXAMPLES::

        sage: import sage.stats.statistics as statistics
        sage: statistics.pvariance([1..6])
        35/12
        sage: statistics.pvariance([e, pi])
        1/4*(pi - e)^2

        sage: import numpy
        sage: x = numpy.array([1,2,3,4,5])
        sage: statistics.pvariance(x)
        2.0
        sage: statistics.pvariance([1.,2.,3.,4.,5.])
        2.00000000000000

        sage: statistics.pvariance([])
        doctest:warning
        ...
        RuntimeWarning: population variance of empty data
        NaN

    TESTS::

        sage: import sage.stats.statistics as sage_statistics
        sage: import statistics as python_statistics
        sage: v = [float(1), float(2), float(3)]
        sage: sage_statistics.pvariance(v) == python_statistics.pvariance(v)
        True
    """
    if is_numpy_type(type(v)):
        import numpy
        if isinstance(v, numpy.ndarray):
            return v.var()

    if not v:
        warn("population variance of empty data", RuntimeWarning)
        return NaN

    mu = mean(v)
    return sum( (x - mu)**2 for x in v) / len(v)
