"""
Long tests for GAP

These stress test the garbage collection inside GAP
"""

from sage.libs.gap.libgap import libgap


def _test_loop_1():
    """
    EXAMPLES::

        sage: from sage.libs.gap._test_long import _test_loop_1
        sage: _test_loop_1()  # long time (up to 25s on sage.math, 2013)
    """
    libgap.collect()
    for i in range(10000):
        G = libgap.CyclicGroup(2)


def _test_loop_2():
    """
    EXAMPLES::

        sage: from sage.libs.gap._test_long import _test_loop_2
        sage: _test_loop_2()  # long time (10s on sage.math, 2013)
    """
    G =libgap.FreeGroup(2)
    a,b = G.GeneratorsOfGroup()
    for i in range(100):
        rel = libgap([a**2, b**2, a*b*a*b])
        H = G / rel
        H1 = H.GeneratorsOfGroup()[0]
        n = H1.Order()
        assert n.sage() == 2

    for i in range(300000):
        n = libgap.Order(H1)


def _test_loop_3():
    """
    EXAMPLES::

        sage: from sage.libs.gap._test_long import _test_loop_3
        sage: _test_loop_3()  # long time (31s on sage.math, 2013)
    """
    G = libgap.FreeGroup(2)
    (a,b) = G.GeneratorsOfGroup()
    for i in range(300000):
        lis=libgap([])
        lis.Add(a ** 2)
        lis.Add(b ** 2)
        lis.Add(b * a)



