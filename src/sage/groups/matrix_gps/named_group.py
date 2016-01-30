"""
Base for Classical Matrix Groups

This module implements the base class for matrix groups that have
various famous names, like the general linear group.

EXAMPLES::

    sage: SL(2, ZZ)
    Special Linear Group of degree 2 over Integer Ring
    sage: G = SL(2,GF(3)); G
    Special Linear Group of degree 2 over Finite Field of size 3
    sage: G.is_finite()
    True
    sage: G.conjugacy_class_representatives()
    (
    [1 0]  [0 2]  [0 1]  [2 0]  [0 2]  [0 1]  [0 2]
    [0 1], [1 1], [2 1], [0 2], [1 2], [2 2], [1 0]
    )
    sage: G = SL(6,GF(5))
    sage: G.gens()
    (
    [2 0 0 0 0 0]  [4 0 0 0 0 1]
    [0 3 0 0 0 0]  [4 0 0 0 0 0]
    [0 0 1 0 0 0]  [0 4 0 0 0 0]
    [0 0 0 1 0 0]  [0 0 4 0 0 0]
    [0 0 0 0 1 0]  [0 0 0 4 0 0]
    [0 0 0 0 0 1], [0 0 0 0 4 0]
    )
"""

##############################################################################
#       Copyright (C) 2006 David Joyner and William Stein <wstein@gmail.com>
#       Copyright (C) 2013 Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
##############################################################################

from sage.structure.unique_representation import UniqueRepresentation
from sage.groups.matrix_gps.matrix_group import (
    MatrixGroup_generic, MatrixGroup_gap )


def normalize_args_vectorspace(*args, **kwds):
    """
    Normalize the arguments that relate to a vector space.

    INPUT:

    Something that defines an affine space. For example

    * An affine space itself:

      - ``A`` -- affine space

    * A vector space:

      - ``V`` -- a vector space

    * Degree and base ring:

      - ``degree`` -- integer. The degree of the affine group, that
        is, the dimension of the affine space the group is acting on.

      - ``ring`` -- a ring or an integer. The base ring of the affine
        space. If an integer is given, it must be a prime power and
        the corresponding finite field is constructed.

      - ``var='a'`` -- optional keyword argument to specify the finite
        field generator name in the case where ``ring`` is a prime
        power.

    OUTPUT:

    A pair ``(degree, ring)``.

    TESTS::

        sage: from sage.groups.matrix_gps.named_group import normalize_args_vectorspace
        sage: A = AffineSpace(2, GF(4,'a'));  A
        Affine Space of dimension 2 over Finite Field in a of size 2^2
        sage: normalize_args_vectorspace(A)
        (2, Finite Field in a of size 2^2)

        sage: normalize_args_vectorspace(2,4)   # shorthand
        (2, Finite Field in a of size 2^2)

        sage: V = ZZ^3;  V
        Ambient free module of rank 3 over the principal ideal domain Integer Ring
        sage: normalize_args_vectorspace(V)
        (3, Integer Ring)

        sage: normalize_args_vectorspace(2, QQ)
        (2, Rational Field)
    """
    from sage.rings.all import ZZ
    if len(args) == 1:
        V = args[0]
        try:
            degree = V.dimension_relative()
        except AttributeError:
            degree = V.dimension()
        ring = V.base_ring()
    if len(args) == 2:
        degree, ring = args
        from sage.rings.integer import is_Integer
        try:
            ring = ZZ(ring)
            from sage.rings.finite_rings.finite_field_constructor import FiniteField
            var = kwds.get('var', 'a')
            ring = FiniteField(ring, var)
        except (ValueError, TypeError):
            pass
    return (ZZ(degree), ring)


class NamedMatrixGroup_generic(UniqueRepresentation, MatrixGroup_generic):

    def __init__(self, degree, base_ring, special, sage_name, latex_string):
        """
        Base class for "named" matrix groups

        INPUT:

        - ``degree`` -- integer. The degree (number of rows/columns of matrices).

        - ``base_ring`` -- rinrg. The base ring of the matrices.

        - ``special`` -- boolean. Whether the matrix group is special,
          that is, elements have determinant one.

        - ``latex_string`` -- string. The latex representation.

        EXAMPLES::

            sage: G = GL(2, QQ)
            sage: from sage.groups.matrix_gps.named_group import NamedMatrixGroup_generic
            sage: isinstance(G, NamedMatrixGroup_generic)
            True
        """
        MatrixGroup_generic.__init__(self, degree, base_ring)
        self._special = special
        self._name_string = sage_name
        self._latex_string = latex_string

    def _an_element_(self):
        """
        Return an element.

        OUTPUT:

        A group element.

        EXAMPLES::

            sage: GL(2, QQ)._an_element_()
            [1 0]
            [0 1]
        """
        return self(1)

    def _repr_(self):
        """
        Return a string representation.

        OUTPUT:

        String.

        EXAMPLES::

            sage: GL(2, QQ)._repr_()
            'General Linear Group of degree 2 over Rational Field'
        """
        return self._name_string

    def _latex_(self):
        """
        Return a LaTeX representation

        OUTPUT:

        String.

        EXAMPLES::

            sage: GL(2, QQ)._latex_()
            'GL(2, \\Bold{Q})'
        """
        return self._latex_string

    def __eq__(self, other):
        """
        Override comparison.

        We need to override the comparison since the named groups
        derive from
        :class:`~sage.structure.unique_representation.UniqueRepresentation`,
        which compare by identity.

        EXAMPLES::

            sage: G = GL(2,3)
            sage: G == MatrixGroup(G.gens())
            True
        """
        return self.__cmp__(other) == 0


class NamedMatrixGroup_gap(NamedMatrixGroup_generic, MatrixGroup_gap):

    def __init__(self, degree, base_ring, special, sage_name, latex_string, gap_command_string):
        """
        Base class for "named" matrix groups using LibGAP

        INPUT:

        - ``degree`` -- integer. The degree (number of rows/columns of matrices).

        - ``base_ring`` -- rinrg. The base ring of the matrices.

        - ``special`` -- boolean. Whether the matrix group is special,
          that is, elements have determinant one.

        - ``latex_string`` -- string. The latex representation.

        - ``gap_command_string`` -- string. The GAP command to construct the matrix group.

        EXAMPLES::

            sage: G = GL(2, GF(3))
            sage: from sage.groups.matrix_gps.named_group import NamedMatrixGroup_gap
            sage: isinstance(G, NamedMatrixGroup_gap)
            True
        """
        from sage.libs.gap.libgap import libgap
        group = libgap.eval(gap_command_string)
        MatrixGroup_gap.__init__(self, degree, base_ring, group)
        self._special = special
        self._gap_string = gap_command_string
        self._name_string = sage_name
        self._latex_string = latex_string

