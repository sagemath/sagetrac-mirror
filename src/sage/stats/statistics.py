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

.. SEEALSO::

    - `numpy statistics routine <https://numpy.org/doc/stable/reference/routines.statistics.html>`_
    - `scipy statistical functions <https://docs.scipy.org/doc/scipy/reference/stats.html>`_
    - `pandas module <https://pandas.pydata.org/docs/index.html>`_
"""
# ***********************************************************************
#          Copyright (C) 2009, Andrew Hou <amhou@uw.edu>
#                        2021, Vincent Delecroix <20100.delecroix@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from statistics import mode, multimode, median, StatisticsError
from warnings import warn

from sage.structure.coerce import is_numpy_type, py_scalar_to_element
from sage.misc.functional import sqrt

def mean(data):
    """
    Return the sample mean of ``data``.

    For empty ``data``, a ``StatisticsError`` is raised.

    If ``data`` is a numpy array, then the specialized method ``mean`` is
    called. Otherwise, the elements are converted to sage objects and their
    mean is computed.

    EXAMPLES::

        sage: import sage.stats.statistics as statistics

        sage: statistics.mean([pi, e])
        1/2*pi + 1/2*e
        sage: statistics.mean([I, sqrt(2), 3/5])
        1/3*sqrt(2) + 1/3*I + 1/5
        sage: statistics.mean([RIF(1.0103,1.0103), RIF(2)])
        1.5051500000000000?
        sage: statistics.mean(range(4))
        3/2

        sage: import numpy
        sage: x = numpy.array([1,2,3,4,5])
        sage: statistics.mean(x)
        3.0

        sage: statistics.mean([])
        Traceback (most recent call last):
        ...
        statistics.StatisticsError: mean requires at least one data point

    TESTS::

        sage: import sage.stats.statistics as sage_statistics
        sage: import statistics as python_statistics
        sage: v = [float(1), float(2), float(3)]
        sage: sage_statistics.mean(v) == python_statistics.mean(v)
        True
    """
    if iter(data) is data:
        data = list(data)
    if not len(data):
        raise StatisticsError(" mean requires at least one data point")
    if is_numpy_type(type(data)):
        import numpy
        if isinstance(data, numpy.ndarray):
            return data.mean()
    return sum(map(py_scalar_to_element, data)) / len(data)


def stdev(data, xbar=None):
    r"""
    Return the square root of the sample variance of ``data``.

    See :func:`variance` for arguments.

    If ``data`` is a numpy array, then the specialized method ``std`` is
    called.

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
        Traceback (most recent call last):
        ...
        statistics.StatisticsError: stdev requires at least two data points
        sage: statistics.stdev(numpy.array([]))
        Traceback (most recent call last):
        ...
        statistics.StatisticsError: stdev requires at least two data points

    TESTS::

        sage: import sage.stats.statistics as sage_statistics
        sage: import statistics as python_statistics
        sage: v = [float(1), float(2), float(3)]
        sage: sage_statistics.stdev(v) == python_statistics.stdev(v)
        True
    """
    if iter(data) is data:
        data = list(data)
    if len(data) < 2:
        raise StatisticsError("stdev requires at least two data points")
    if is_numpy_type(type(data)):
        import numpy
        if isinstance(data, numpy.ndarray):
            return data.std(ddof=1)
    return sqrt(variance(data, xbar))


def pstdev(data, mu=None):
    """
    Return the square root of the population variance of ``data``.

    See :func:`pvariance` for arguments.

    If ``data`` is a numpy array, then the specialized method ``std`` is
    called.

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
        Traceback (most recent call last):
        ...
        statistics.StatisticsError: pstdev requires at least one data point
        sage: statistics.pstdev(numpy.array([]))
        Traceback (most recent call last):
        ...
        statistics.StatisticsError: pstdev requires at least one data point

    TESTS::
    
        sage: import sage.stats.statistics as sage_statistics
        sage: import statistics as python_statistics
        sage: v = [float(1), float(2), float(3)]
        sage: sage_statistics.pstdev(v) == python_statistics.pstdev(v)
        True
    """
    if iter(data) is data:
        data = list(data)
    if not len(data):
        raise StatisticsError("pstdev requires at least one data point")
    if is_numpy_type(type(data)):
        import numpy
        if isinstance(data, numpy.ndarray):
            return data.std()
    return sqrt(pvariance(data, mu))


def variance(data, xbar=None):
    r"""
    Return the sample variance of ``data``.

    If you have already calculated the mean of your data, you can pass it as
    the optional second argument ``xbar`` to avoid recalculating it.

    If ``data`` is a numpy array, then the specialized method ``var`` is
    called. Otherwise, the elements are converted to sage objects and their
    population variance is computed.

    EXAMPLES::

        sage: import sage.stats.statistics as statistics
        sage: statistics.variance([1..6])
        7/2
        sage: m = statistics.mean([1..6])
        sage: statistics.variance([1..6], m)
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

        sage: statistics.variance([1])
        Traceback (most recent call last):
        ...
        statistics.StatisticsError: variance requires at least two data points
        sage: statistics.variance(numpy.array([1]))
        Traceback (most recent call last):
        ...
        statistics.StatisticsError: variance requires at least two data points

    TESTS::

        sage: import sage.stats.statistics as sage_statistics
        sage: import statistics as python_statistics
        sage: v = [float(1), float(2), float(3)]
        sage: sage_statistics.variance(v) == python_statistics.variance(v)
        True
    """
    if iter(data) is data:
        data = list(data)
    if len(data) < 2:
        raise StatisticsError('variance requires at least two data points')
    if is_numpy_type(type(data)):
        import numpy
        if isinstance(data, numpy.ndarray):
            return data.var(ddof=1)
    if xbar is None:
        xbar = mean(data)
    xbar = py_scalar_to_element(xbar)
    return sum( (py_scalar_to_element(x) - xbar)**2 for x in data) / (len(data) - 1)


def pvariance(data, mu=None):
    r"""
    Return the population variance of ``data``.

    The optional argument ``mu``, if given, should be the mean of the data. If
    it is missing or None, the mean is automatically calculated.

    If ``data`` is a numpy array, then the specialized method ``var`` is
    called. Otherwise, the elements are converted to sage objects and their
    population variance is computed.

    EXAMPLES::

        sage: import sage.stats.statistics as statistics
        sage: statistics.pvariance([1..6])
        35/12
        sage: mu = statistics.mean([1..6])
        sage: statistics.pvariance([1..6], mu)
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
        Traceback (most recent call last):
        ...
        statistics.StatisticsError: pvariance requires at least one data point
        sage: statistics.pvariance(numpy.array([]))
        Traceback (most recent call last):
        ...
        statistics.StatisticsError: pvariance requires at least one data point

    TESTS::

        sage: import sage.stats.statistics as sage_statistics
        sage: import statistics as python_statistics
        sage: v = [float(1), float(2), float(3)]
        sage: sage_statistics.pvariance(v) == python_statistics.pvariance(v)
        True
    """
    if iter(data) is data:
        data = list(data)
    if len(data) < 1:
        raise StatisticsError('pvariance requires at least one data point')
    if is_numpy_type(type(data)):
        import numpy
        if isinstance(data, numpy.ndarray):
            return data.var()
    if mu is None:
        mu = mean(data)
    mu = py_scalar_to_element(mu)
    return sum( (py_scalar_to_element(x) - mu)**2 for x in data) / len(data)
