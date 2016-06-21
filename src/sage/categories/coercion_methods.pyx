"""
Coercion methods for categories

The purpose of this Cython module is to hold special coercion methods,
which are inserted by their respective categories.
"""

from __future__ import absolute_import, division, print_function

from sage.structure.element cimport Element
cimport cython


@cython.binding
def _mul_parent(self, other):
    r"""
    Return the product of the two elements, calculated using
    the ``product`` method of the parent.

    This is the default implementation of ``_mul_`` if
    ``product`` is implemented in the parent.

    INPUT:

    - ``other`` -- an element of the parent of ``self``

    OUTPUT:

    - an element of the parent of ``self``

    EXAMPLES::

        sage: S = Semigroups().example("free")
        sage: x = S('a'); y = S('b')
        sage: x._mul_parent(y)
        'ab'

    .. SEEALSO::

        - :meth:`Magmas.ElementMethods._mul_`
        - :meth:`Magmas.ElementMethods._mul_parent`
        - :meth:`Magmas.ParentMethods.product`

    This is :meth:`Magmas.ElementMethods._mul_parent`, implemented as
    a Cython method in :mod:`sage.categories.coercion_methods`::

        sage: x._mul_parent.im_func is Magmas.ElementMethods._mul_parent.im_func
        True
        sage: x._mul_parent.im_func is sage.categories.coercion_methods._mul_parent
        True
    """
    return (<Element>self)._parent.product(self, other)
