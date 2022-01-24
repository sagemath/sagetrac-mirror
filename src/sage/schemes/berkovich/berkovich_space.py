r"""
Berkovich Space over `\CC_p`

The Berkovich affine line is the set of seminorms on `\CC_p[x]`,
with the weakest topology that makes the map `| \cdot | \to |f|` continuous
for all `f \in \CC_p[x]`. The Berkovich projective line is the
one-point compactification of the Berkovich affine line.

The two main classes are :class:`Berkovich_Cp_Affine` and
:class:`Berkovich_Cp_Projective`, which implement the affine and
projective lines, respectively.

:class:`Berkovich_Cp_Affine` and :class:`Berkovich_Cp_Projective`
take as input one of the following: the prime `p`, a finite
extension of `\QQ_p`, or a number field and a place.

For an exposition of Berkovich space over `\CC_p`, see Chapter 6 of [Ben2019]_. For a more
involved exposition, see Chapter 1 and 2 of [BR2010]_.

AUTHORS:

 - Alexander Galarraga (2020-06-22): initial implementation

"""

#*****************************************************************************
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.schemes.berkovich.berkovich_cp_element import (Berkovich_Element_Cp_Affine,
                                                         Berkovich_Element_Cp_Projective)
from sage.structure.parent import Parent
from sage.schemes.affine.affine_space import is_AffineSpace
from sage.schemes.projective.projective_space import is_ProjectiveSpace, ProjectiveSpace
from sage.structure.unique_representation import UniqueRepresentation
from sage.categories.number_fields import NumberFields
import sage.rings.abc
from sage.rings.integer_ring import ZZ
from sage.rings.padics.factory import Qp
from sage.rings.rational_field import QQ
from sage.rings.number_field.number_field_ideal import NumberFieldFractionalIdeal
from sage.categories.topological_spaces import TopologicalSpaces

def is_Berkovich(space):
    """
    Checks if ``space`` is a Berkovich space.

    OUTPUT:

    - ``True`` if ``space`` is a Berkovich space.
    - ``False`` otherwise.

    EXAMPLES::

        sage: B = Berkovich_Cp_Projective(3)
        sage: from sage.schemes.berkovich.berkovich_space import is_Berkovich
        sage: is_Berkovich(B)
        True
    """
    return isinstance(space, Berkovich)

def is_Berkovich_Cp(space):
    """
    Checks if ``space`` is a Berkovich space over ``Cp``.

    OUTPUT:

    - ``True`` if ``space`` is a Berkovich space over ``Cp``.
    - ``False`` otherwise.

    EXAMPLES::

        sage: B = Berkovich_Cp_Projective(3)
        sage: from sage.schemes.berkovich.berkovich_space import is_Berkovich
        sage: is_Berkovich(B)
        True
    """
    return isinstance(space, Berkovich_Cp)

class Berkovich(UniqueRepresentation, Parent):
    """
    The parent class for any Berkovich space
    """
    pass

class Berkovich_Cp(Berkovich):
    """
    Abstract parent class for Berkovich space over ``Cp``.
    """

    def residue_characteristic(self):
        """
        The residue characteristic of the ``base``.

        EXAMPLES::

            sage: B = Berkovich_Cp_Projective(3)
            sage: B.prime()
            3

        ::

            sage: R.<x> = QQ[]
            sage: A.<a> = NumberField(x^3 + 20)
            sage: ideal = A.ideal(-1/2*a^2 + a - 3)
            sage: B = Berkovich_Cp_Affine(A, ideal)
            sage: B.residue_characteristic()
            7
        """
        return self._p

    prime = residue_characteristic

    def is_padic_base(self):
        """
        Return ``True`` if this Berkovich space is backed by a p-adic field.

        OUTPUT:

        - ``True`` if this Berkovich space was created with a p-adic field.
        - ``False`` otherwise.

        EXAMPLES::

            sage: B = Berkovich_Cp_Affine(Qp(3))
            sage: B.is_padic_base()
            True

       ::

            sage: B = Berkovich_Cp_Affine(QQ, 3)
            sage: B.is_padic_base()
            False
        """
        return self._base_type == 'padic field'

    def is_number_field_base(self):
        """
        Return ``True`` if this Berkovich space is backed by a p-adic field.

        OUTPUT:

        - ``True`` if this Berkovich space was created with a number field.
        - ``False`` otherwise.

        EXAMPLES::

            sage: B = Berkovich_Cp_Affine(Qp(3))
            sage: B.is_number_field_base()
            False

        ::

            sage: B = Berkovich_Cp_Affine(QQ, 3)
            sage: B.is_number_field_base()
            True
        """
        return self._base_type == 'number field'

    def ideal(self):
        r"""
        The ideal which defines an embedding of the ``base_ring`` into `\CC_p`.

        If this Berkovich space is backed by a p-adic field, then an embedding is
        already specified, and this returns ``None``.

        OUTPUT:

        - An ideal of a ``base_ring`` if ``base_ring`` is a number field.

        - A prime of `\QQ` if ``base_ring`` is `\QQ`.

        - ``None`` if ``base_ring`` is a p-adic field.

        EXAMPLES::

            sage: R.<z> = QQ[]
            sage: A.<a> = NumberField(z^2 + 1)
            sage: ideal = A.prime_above(5)
            sage: B = Berkovich_Cp_Projective(A, ideal)
            sage: B.ideal()
            Fractional ideal (-a - 2)

        ::

            sage: B = Berkovich_Cp_Projective(QQ, 3)
            sage: B.ideal()
            3

        ::

            sage: B = Berkovich_Cp_Projective(Qp(3))
            sage: B.ideal() is None
            True
        """
        return self._ideal

    def __eq__(self,right):
        """
        Equality operator.

        EXAMPLES::

            sage: B = Berkovich_Cp_Affine(3)
            sage: A.<x> = Qq(27)
            sage: C = Berkovich_Cp_Affine(A)
            sage: B == C
            True

        ::

            sage: R.<x> = QQ[]
            sage: A.<a> = NumberField(x^2 + 1)
            sage: A_ideal = A.prime_above(2)
            sage: B.<b> = NumberField(x^4 + 1)
            sage: B_ideal = B.prime_above(2)
            sage: C = Berkovich_Cp_Projective(A, A_ideal)
            sage: D = Berkovich_Cp_Projective(B, B_ideal)
            sage: C == D
            False

        ::

            sage: C = Berkovich_Cp_Affine(A, A_ideal)
            sage: D = Berkovich_Cp_Affine(B, B_ideal)
            sage: C == D
            False

        ::

            sage: A_ideal_2 = A.prime_above(5)
            sage: E = Berkovich_Cp_Affine(A, A_ideal_2)
            sage: C == E
            False
        """
        if not isinstance(right, Berkovich_Cp):
            return False
        if self._base_type != right._base_type:
            return False
        if self._base_type == 'padic field':
            return self.prime() == right.prime()
        else:
            return self.base() == right.base() and self.ideal() == right.ideal()

    def __ne__(self,right):
        """
        Inequality operator.

        EXAMPLES::

            sage: B = Berkovich_Cp_Affine(5)
            sage: A.<x> = Qq(25)
            sage: C = Berkovich_Cp_Affine(A)
            sage: B != C
            False
        """
        return not(self == right)

    def __hash__(self):
        """
        Hash function.

        EXAMPLES::

            sage: hash(Berkovich_Cp_Projective(3))
            3

        ::

            sage: R.<z> = QQ[]
            sage: A.<a> = NumberField(z^2 + 1)
            sage: B = Berkovich_Cp_Projective(A, A.primes_above(5)[0])
            sage: C = Berkovich_Cp_Projective(A, A.primes_above(5)[1])
            sage: hash(B) != hash(C)
            True
        """
        if self._base_type == 'padic field':
            return hash(self.prime())
        return hash(self.ideal())

class Berkovich_Cp_Affine(Berkovich_Cp):
    r"""
    The Berkovich affine line over `\CC_p`.

    The Berkovich affine line is the set of seminorms on `\CC_p[x]`,
    with the weakest topology such that the map `| \cdot | \to |f|` is continuous
    for all `f \in \CC_p[x]`.

    We can represent the Berkovich affine line in two separate ways:
    either using a p-adic field to represent elements or using
    a number field to represent elements while storing an ideal
    of the ring of integers of the number field, which specifies
    an embedding of the number field into `\CC_p`. See the examples.

    INPUT:

    - ``base`` -- Three cases:

      * a prime number `p`. Centers of elements are then represented
        as points of `\QQ_p`.

      * `\QQ_p` or a finite extension of `\QQ_p`. Centers of elements
        are then represented as points of ``base``.

      * A number field `K`. Centers of elements are then represented
        as points of `K`.

    - ``ideal`` -- (optional) a prime ideal of ``base``. Must be
      specified if a number field is passed to ``base``, otherwise
      it is ignored.

    EXAMPLES::

        sage: B = Berkovich_Cp_Affine(3); B
        Affine Berkovich line over Cp(3) of precision 20

    We can create elements::

        sage: B(-2)
        Type I point centered at 1 + 2*3 + 2*3^2 + 2*3^3 + 2*3^4 + 2*3^5
        + 2*3^6 + 2*3^7 + 2*3^8 + 2*3^9 + 2*3^10 + 2*3^11 + 2*3^12 + 2*3^13
        + 2*3^14 + 2*3^15 + 2*3^16 + 2*3^17 + 2*3^18 + 2*3^19 + O(3^20)

    ::

        sage: B(1, 2)
        Type III point centered at 1 + O(3^20) of radius 2.00000000000000

    For details on element creation, see the documentation
    of :class:`Berkovich_Element_Cp_Affine`. Initializing by
    passing in `\QQ_p` looks the same::

        sage: B = Berkovich_Cp_Affine(Qp(3)); B
        Affine Berkovich line over Cp(3) of precision 20

    However, this method allows for more control over behind-the-scenes conversion::

        sage: B = Berkovich_Cp_Affine(Qp(3, 1)); B
        Affine Berkovich line over Cp(3) of precision 1

        sage: B(1/2)
        Type I point centered at 2 + O(3)

    Note that this point has very low precision, as ``B`` was initialized
    with a p-adic field of capped-relative precision one. For high precision,
    pass in a high precision p-adic field::

        sage: B = Berkovich_Cp_Affine(Qp(3, 1000)); B
        Affine Berkovich line over Cp(3) of precision 1000

    Points of Berkovich space can be created from points of
    extensions of `\QQ_p`::

        sage: B = Berkovich_Cp_Affine(3)
        sage: A.<a> = Qp(3).extension(x^3 - 3)
        sage: B(a)
        Type I point centered at a + O(a^61)

    For exact computation, a number field can be used::

        sage: R.<x> = QQ[]
        sage: A.<a> = NumberField(x^3 + 20)
        sage: ideal = A.prime_above(3)
        sage: B = Berkovich_Cp_Affine(A, ideal); B
        Affine Berkovich line over Cp(3), with base Number
        Field in a with defining polynomial x^3 + 20

    Number fields have a major advantage of exact computation.

    Number fields also have added functionality. Arbitrary extensions of
    `\QQ` are supported, while there is currently limited functionality
    for extensions of `\QQ_p`. As seen above, constructing a Berkovich
    space backed by a number field requires specifying an ideal of the
    ring of integers of the number field. Specifying the ideal uniquely
    specifies an embedding of the number field into `\CC_p`.

    Unlike in the case where Berkovich space is backed by a p-adic
    field, any point of a Berkovich space backed by a number field
    must be centered at a point of that number field::

        sage: R.<x> = QQ[]
        sage: A.<a> = NumberField(x^3 + 20)
        sage: ideal = A.prime_above(3)
        sage: B = Berkovich_Cp_Affine(A, ideal)
        sage: C.<c> = NumberField(x^2 + 1)
        sage: B(c)
        Traceback (most recent call last):
        ...
        ValueError: could not convert c to Number Field in a
        with defining polynomial x^3 + 20

    TESTS::

        sage: A.<x> = AffineSpace(Qp(3), 1)
        sage: Berkovich_Cp_Affine(A)
        Affine Berkovich line over Cp(3) of precision 20

    ::

        sage: B = Berkovich_Cp_Projective(3)
        sage: TestSuite(B).run()
    """

    Element = Berkovich_Element_Cp_Affine

    def __init__(self, base, ideal=None):
        """
        The Python constructor.

        EXAMPLES::

            sage: Berkovich_Cp_Affine(3)
            Affine Berkovich line over Cp(3) of precision 20
        """
        if base in ZZ:
            if base.is_prime():
                base = Qp(base)  # change to Qpbar
            else:
                raise ValueError("non-prime passed into Berkovich space")
        if is_AffineSpace(base):
            base = base.base_ring()
        if base in NumberFields():
            if ideal is None:
                raise ValueError('passed a number field but not an ideal')
            if base is not QQ:
                if not isinstance(ideal, NumberFieldFractionalIdeal):
                    raise ValueError('ideal was not an ideal of a number field')
                if ideal.number_field() != base:
                    raise ValueError('passed number field ' + \
                        '%s but ideal was an ideal of %s' %(base, ideal.number_field()))
                prime = ideal.smallest_integer()
            else:
                if ideal not in QQ:
                    raise ValueError('ideal was not an element of QQ')
                prime = ideal
            if not ideal.is_prime():
                raise ValueError('passed non prime ideal')
            self._base_type = 'number field'
        elif isinstance(base, sage.rings.abc.pAdicField):  # change base to Qpbar
            prime = base.prime()
            ideal = None
            self._base_type = 'padic field'
        else:
            raise ValueError("base of Berkovich Space must be a padic field " + \
                "or a number field")
        self._ideal = ideal
        self._p = prime
        Parent.__init__(self, base=base, category=TopologicalSpaces())

    def _repr_(self):
        """
        String representation of this Berkovich Space.

        EXAMPLES::

            sage: B = Berkovich_Cp_Affine(3)
            sage: B
            Affine Berkovich line over Cp(3) of precision 20

        ::

            sage: R.<z> = QQ[]
            sage: A.<a> = NumberField(z^2 + 1)
            sage: ideal = A.prime_above(3)
            sage: Berkovich_Cp_Affine(A, ideal)
            Affine Berkovich line over Cp(3), with base Number Field
            in a with defining polynomial z^2 + 1
        """
        if self._base_type == 'padic field':
            return "Affine Berkovich line over Cp(%s) of precision %s" %(self.prime(),\
                self.base().precision_cap())
        else:
            return  "Affine Berkovich line over Cp(%s), with base %s" %(self.prime(),\
                self.base())

    def _latex_(self):
        r"""
        LaTeX representation of this Berkovich Space.

        EXAMPLES:

            sage: B = Berkovich_Cp_Affine(3)
            sage: latex(B)
            \text{Affine Berkovich line over } \Bold{C}_{3}
        """
        return r"\text{Affine Berkovich line over } \Bold{C}_{%s}" % (self.prime())


class Berkovich_Cp_Projective(Berkovich_Cp):
    r"""
    The Berkovich projective line over `\CC_p`.

    The Berkovich projective line is the one-point compactification
    of the Berkovich affine line.

    We can represent the Berkovich projective line in two separate ways:
    either using a p-adic field to represent elements or using
    a number field to represent elements while storing an ideal
    of the ring of integers of the number field, which specifies
    an embedding of the number field into `\CC_p`. See the examples.

    INPUT:

    - ``base`` -- Three cases:

      * a prime number `p`. Centers of elements are then represented
        as points of projective space of dimension 1 over `\QQ_p`.

      * `\QQ_p` or a finite extension of `\QQ_p`. Centers of elements
        are then represented as points of projective space of dimension 1
        over ``base``.

      * A number field `K`. Centers of elements are then represented
        as points of projective space of dimension 1 over ``base``.

    - ``ideal`` -- (optional) a prime ideal of ``base``. Must be
      specified if a number field is passed to ``base``, otherwise
      it is ignored.

    EXAMPLES::

        sage: B = Berkovich_Cp_Projective(3); B
        Projective Berkovich line over Cp(3) of precision 20

    Elements can be constructed::

        sage: B(1/2)
        Type I point centered at (2 + 3 + 3^2 + 3^3 + 3^4 + 3^5
        + 3^6 + 3^7 + 3^8 + 3^9 + 3^10 + 3^11 + 3^12 + 3^13 + 3^14
        + 3^15 + 3^16 + 3^17 + 3^18 + 3^19 + O(3^20) : 1 + O(3^20))

    ::

        sage: B(2, 1)
        Type II point centered at (2 + O(3^20) : 1 + O(3^20)) of radius 3^0

    For details about element construction, see the documentation of
    :class:`Berkovich_Element_Cp_Projective`. Initializing a Berkovich projective
    line by passing in a p-adic space looks the same::

        sage: B = Berkovich_Cp_Projective(Qp(3)); B
        Projective Berkovich line over Cp(3) of precision 20

    However, this method allows for more control over
    behind-the-scenes conversion::

        sage: S = Qp(3, 1)
        sage: B = Berkovich_Cp_Projective(S); B
        Projective Berkovich line over Cp(3) of precision 1

        sage: Q1 = B(1/2); Q1
        Type I point centered at (2 + O(3) : 1 + O(3))

    Note that this point has very low precision, as S has low
    precision cap. Berkovich space can also be created over
    a number field, as long as an ideal is specified::

        sage: R.<x> = QQ[]
        sage: A.<a> = NumberField(x^2 + 1)
        sage: ideal = A.prime_above(2)
        sage: B = Berkovich_Cp_Projective(A, ideal); B
        Projective Berkovich line over Cp(2), with base
        Number Field in a with defining polynomial x^2 + 1

    Number fields have the benefit that computation is exact,
    but lack support for all of `\CC_p`.

    Number fields also have the advantage of added functionality,
    as arbitrary extensions of `\QQ` can be constructed while
    there is currently limited functionality for extensions of `\QQ_p`.
    As seen above, constructing a Berkovich space backed by a number
    field requires specifying an ideal of the ring of integers
    of the number field. Specifying the ideal uniquely specifies
    an embedding of the number field into `\CC_p`.

    Unlike in the case where Berkovich space is backed by a p-adic
    field, any point of a Berkovich space backed by a number field
    must be centered at a point of that number field::

        sage: R.<x> = QQ[]
        sage: A.<a> = NumberField(x^3 + 20)
        sage: ideal = A.prime_above(3)
        sage: B = Berkovich_Cp_Projective(A, ideal)
        sage: C.<c> = NumberField(x^2 + 1)
        sage: B(c)
        Traceback (most recent call last):
        ...
        TypeError: could not convert c to Projective Space
        of dimension 1 over Number Field in a with defining polynomial x^3 + 20

    TESTS::

        sage: B = Berkovich_Cp_Projective(3)
        sage: TestSuite(B).run()
    """

    Element = Berkovich_Element_Cp_Projective

    def __init__(self, base, ideal=None):
        """
        The Python constructor.

        EXAMPLES::

            sage: Berkovich_Cp_Projective(3)
            Projective Berkovich line over Cp(3) of precision 20
        """
        if base in ZZ:
            if base.is_prime():
                base = ProjectiveSpace(Qp(base), 1)
            else:
                raise ValueError("non-prime passed into Berkovich space")
        if base in NumberFields() or isinstance(base, sage.rings.abc.pAdicField):
            base = ProjectiveSpace(base, 1)
        if not is_ProjectiveSpace(base):
            try:
                base = ProjectiveSpace(base)
            except:
                raise ValueError("base of projective Berkovich space must be projective space")
        if not isinstance(base.base_ring(), sage.rings.abc.pAdicField):
            if base.base_ring() not in NumberFields():
                raise ValueError("base of projective Berkovich space must be " + \
                    "projective space over Qp or a number field")
            else:
                if ideal is None:
                    raise ValueError('passed a number field but not an ideal')
                if base.base_ring() is not QQ:
                    if not isinstance(ideal, NumberFieldFractionalIdeal):
                        raise ValueError('ideal was not a number field ideal')
                    if ideal.number_field() != base.base_ring():
                        raise ValueError('passed number field ' + \
                            '%s but ideal was an ideal of %s' %(base.base_ring(), ideal.number_field()))
                    prime = ideal.smallest_integer()
                else:
                    if ideal not in QQ:
                        raise ValueError('ideal was not an element of QQ')
                    prime = ideal
                if not ideal.is_prime():
                    raise ValueError('passed non prime ideal')
                self._base_type = 'number field'
        else:
            prime = base.base_ring().prime()
            ideal = None
            self._base_type = 'padic field'
        if base.dimension_relative() != 1:
            raise ValueError("base of projective Berkovich space must be " + \
                "projective space of dimension 1 over Qp or a number field")
        self._p = prime
        self._ideal = ideal
        Parent.__init__(self, base = base, category=TopologicalSpaces())

    def base_ring(self):
        r"""
        The base ring of this Berkovich Space.

        OUTPUT: A field.

        EXAMPLES::

            sage: B = Berkovich_Cp_Projective(3)
            sage: B.base_ring()
            3-adic Field with capped relative precision 20

        ::

            sage: C = Berkovich_Cp_Projective(ProjectiveSpace(Qp(3, 1), 1))
            sage: C.base_ring()
            3-adic Field with capped relative precision 1

        ::

            sage: R.<x> = QQ[]
            sage: A.<a> = NumberField(x^3 + 20)
            sage: ideal = A.prime_above(3)
            sage: D = Berkovich_Cp_Projective(A, ideal)
            sage: D.base_ring()
            Number Field in a with defining polynomial x^3 + 20
        """
        return self.base().base_ring()

    def _repr_(self):
        """
        String representation of this Berkovich Space.

        EXAMPLES::

            sage: B = Berkovich_Cp_Projective(3)
            sage: B
            Projective Berkovich line over Cp(3) of precision 20

        ::

            sage: R.<x> = QQ[]
            sage: A.<a> = NumberField(x^2 + 1)
            sage: v = A.ideal(a + 1)
            sage: Berkovich_Cp_Projective(A, v)
            Projective Berkovich line over Cp(2), with base Number Field in a with defining polynomial x^2 + 1
        """
        if self._base_type == 'padic field':
            return "Projective Berkovich line over Cp(%s) of precision %s" %(self.prime(),\
                self.base().base_ring().precision_cap())
        else:
            return  "Projective Berkovich line over Cp(%s), with base %s" %(self.prime(),\
                self.base().base_ring())

    def _latex_(self):
        r"""
        LaTeX representation of this Berkovich Space.

        EXAMPLES:

            sage: B = Berkovich_Cp_Projective(3)
            sage: latex(B)
            \text{Projective Berkovich line over } \Bold{C}_{3}
        """
        return r"\text{Projective Berkovich line over } \Bold{C}_{%s}" %(self.prime())

    def convex_hull(self, points, labels=None):
        """
        Return a bipartite graph containing the convex hull and options for plotting the graph.

        The convex hull of a set of points is the smallest
        path-connected space that contains all the points,
        or equivalently, the union of all intervals starting
        and ending at points in the set.

        For reference, 0, the Gauss point, and infinity are added
        to the bipartite graph. The convex hull is highlighted in red.

        By default, the vertices are labeled ``P1`` through ``Pn``,
        where ``n`` is the number of points in ``points``. To change the
        labels for the vertices, use the optional parameter ``labels``.

        INPUT:

        - points -- A list of points of this Berkovich space.

        - labels -- (optional) A list of strings to use as labels for the vertices
          bipartite graph.

        OUTPUT: A bipartite graph and options for plotting the graph

        EXAMPLES:

        To plot the graph showing the convex hull, use the specified options.
        Note that the options are a dictionary, and so must be passed to the plot
        function using **::

            sage: B = Berkovich_Cp_Projective(3)
            sage: Q1 = B(1/3)
            sage: Q2 = B(1/3, 3)
            sage: Q3 = B(3)
            sage: G, options = B.convex_hull([Q1, Q2, Q3])
            sage: G.plot(**options)
            Graphics object consisting of 12 graphics primitives

        ::

            sage: B = Berkovich_Cp_Projective(3)
            sage: Q1 = B(1)
            sage: Q2 = B(1/81,1)
            sage: Q4 = B(1/243, 1/3)
            sage: Q3 = B(1/3)
            sage: Q5 = B(9, 1/3)
            sage: Q6 = B(1/2)
            sage: Q7 = B(7)
            sage: G, options = B.convex_hull([Q1, Q2, Q3, Q4, Q5, Q6, Q7])
            sage: G.plot(**options)
            Graphics object consisting of 28 graphics primitives

        We can use ``labels`` to specify the naming of the vertices::

            sage: B = Berkovich_Cp_Projective(5)
            sage: Q1 = B(1/5)
            sage: Q2 = B(5, 1/5)
            sage: Q3 = B(25, 1/25)
            sage: convex_hull, options = B.convex_hull([Q1, Q2, Q3], ['Q1', 'Q2', 'Q3'])
            sage: convex_hull.vertices()
            ['0', 'Gauss', 'Q1', 'Q1 ^ Q2', 'Q2', 'Q3', 'oo']

        TESTS::

            sage: B = Berkovich_Cp_Projective(3)
            sage: Q1 = B(1/3)
            sage: Q2 = B(1/3, 3)
            sage: Q3 = B(3)
            sage: convex_hull =  B.convex_hull([Q1, Q2, Q3])
            sage: convex_hull[0].adjacency_matrix()
            [0 1 0 0 0 0]
            [0 0 0 1 0 0]
            [0 0 0 1 0 0]
            [0 0 0 0 0 1]
            [0 1 0 0 0 0]
            [0 0 0 0 0 0]

        ::

            sage: B = Berkovich_Cp_Projective(3)
            sage: Q1 = B(1)
            sage: Q2 = B(1/81,1)
            sage: Q4 = B(1/243, 1/3)
            sage: Q3 = B(1/3)
            sage: Q5 = B(9, 1/3)
            sage: Q6 = B(1/2)
            sage: Q7 = B(7)
            sage: convex_hull = B.convex_hull([Q1, Q2, Q3, Q4, Q5, Q6, Q7])
            sage: convex_hull[0].adjacency_matrix()
            [0 0 0 0 0 0 0 0 0 0 1 0 0 0]
            [0 0 0 0 1 0 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 1 0 0 0 0 0 0 0]
            [0 0 0 0 0 1 0 0 0 0 0 0 0 0]
            [0 0 0 1 0 0 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0 0 0 0 0 0 1]
            [0 1 0 0 0 0 0 0 0 0 0 0 0 0]
            [0 0 0 1 0 0 0 0 0 0 0 0 0 0]
            [0 0 0 0 1 0 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 1 0 0 0 0 0 0 0 0]
            [0 1 0 0 0 0 0 0 0 0 0 0 0 0]
            [0 1 0 0 0 0 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 1 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0 0 0 0 0 0 0]
        """
        if not (isinstance(points, list) or isinstance(points, tuple)):
            raise ValueError("points must be a list")
        if not (isinstance(labels, list) or isinstance(labels, tuple) or labels == None):
            raise ValueError("labels must be a list")
        if labels != None:
            if len(labels) != len(points):
                raise ValueError("the number of labels and the number of points must be the same")
            for i in labels:
                if not isinstance(i, str):
                    raise ValueError("labels must be a list of strings")
        for i in points:
            if not isinstance(i, Berkovich_Element_Cp_Projective):
                raise ValueError("input to convex_hull must be a list of points of Berkovich space")
            if i.parent() != self:
                raise ValueError("input to convex_hull must be a list of points of this Berkovich space")
        if labels is None:
            labels = ["P" + str(i+1) for i in range(len(points))]
        point_to_name = {}
        name_to_highlight = {}
        vertices = points[:]
        for i in range(len(points)):
            point_to_name[points[i]] = labels[i]
            name_to_highlight[point_to_name[points[i]]] = True
        zero = self(0)
        point_to_name[zero] = '0'
        infinity = self((1,0))
        point_to_name[infinity] = 'oo'
        gauss = self(ZZ(0), ZZ(1))
        point_to_name[gauss] = 'Gauss'

        # we add zero, infinity, and the gauss point for reference
        for point in [zero, infinity, gauss]:
            if point not in vertices:
                vertices.append(point)
                name_to_highlight[point_to_name[point]] = False
            else:
                name_to_highlight[point_to_name[point]] = True

        # we add joins if necessary
        for i in range(len(points)):
            for j in range(i+1,len(points)):
                join = points[i].join(points[j])
                if join not in vertices:
                    vertices.append(join)
                    point_to_name[join] = point_to_name[points[i]] + " ^ " + point_to_name[points[j]]
                    name_to_highlight[point_to_name[join]] = True

        points_in_below_gauss = [point for point in vertices if point.lt(gauss)] + [gauss]
        points_in_above_gauss = [point for point in vertices if not point.lt(gauss)]

        # if there are highlighted points both above and below the gauss point we must highlight the gauss point
        if len(points_in_below_gauss) != 2 and len(points_in_above_gauss) != 2:
            name_to_highlight[point_to_name[gauss]] = True

        E = {}
        red_colored_edges = []
        for start in vertices:
            outgoing = {}
            end = infinity
            if start == end:
                continue
            # we now make better guesses for end until we find the unique end which is 'closest' to start -
            # for which the interval [start, end] contains no other points of our tree
            good_end = False
            while not good_end:
                good_end = True
                for v in vertices:
                    if v != start and v != end:
                        if v.contained_in_interval(start, end):
                            good_end = False
                            end = v
                            break
            outgoing[point_to_name[end]] = point_to_name[start] + ' to ' + point_to_name[end]
            E[point_to_name[start]] = outgoing
            if name_to_highlight[point_to_name[start]] and name_to_highlight[point_to_name[end]] and outgoing:
                t1 = point_to_name[start]
                t2 = point_to_name[end]
                t3 = outgoing[point_to_name[end]]
                red_colored_edges.append((t1, t2, t3))
        edge_colors = {'#FF0000':red_colored_edges}
        from sage.graphs.digraph import DiGraph
        G = DiGraph(E)
        red_colored_vertices = []
        for vertex in G:
            if name_to_highlight[vertex]:
                red_colored_vertices.append(vertex)
        vertex_colors = {'#FF0000': red_colored_vertices}
        options = {'layout':'tree',
                'edge_colors':edge_colors,
                'vertex_colors':vertex_colors,
                'tree_orientation': 'down',
                'tree_root': 'oo'}
        return G, options
