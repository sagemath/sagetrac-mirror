"""
Field of Arbitrary Precision Real Number Intervals
"""

from sage.rings.real_mpfi import RealIntervalField_class, RealIntervalFieldElement

def is_RealIntervalField(x):
    """
    Check if ``x`` is a :class:`RealIntervalField_class`.

    EXAMPLES::

        sage: from sage.rings.real_interval_field import is_RealIntervalField as is_RIF
        sage: is_RIF(RIF)
        True
    """
    return isinstance(x, RealIntervalField_class)

def is_RealIntervalFieldElement(x):
    """
    Check if ``x`` is a :class:`RealIntervalFieldElement`.

    EXAMPLES::

        sage: from sage.rings.real_interval_field import is_RealIntervalFieldElement as is_RIFE
        sage: is_RIFE(RIF(2.5))
        True
    """
    return isinstance(x, RealIntervalFieldElement)
