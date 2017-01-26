"""
test the zip behaviour in doc::

    sage: zip([1,2,3,4],['a','b','c','d'])
    <itertools.izip object at ...>
"""


def test(a, b):
    """
    test the zip behaviour in code::

        sage: from sage.combinat.test_zip import test
        sage: test([1,2,3,4],['a','b','c','d'])
        <itertools.izip object at ...>
    """
    return zip(a, b)

         
