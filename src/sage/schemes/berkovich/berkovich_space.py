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
from sage.schemes.projective.projective_space import is_ProjectiveSpace, ProjectiveSpace
from sage.structure.unique_representation import UniqueRepresentation
from sage.categories.number_fields import NumberFields
from sage.rings.integer_ring import ZZ
from sage.rings.padics.factory import Qp
from sage.rings.rational_field import QQ
from sage.rings.number_field.number_field_ideal import NumberFieldFractionalIdeal
from sage.rings.padics.generic_nodes import is_pAdicField
from sage.categories.topological_spaces import TopologicalSpaces

def is_Berkovich(space):
    """
    Checks if ``space`` is a Berkovich space.

    OUTPUT:

    - ``True`` if ``space`` is a Berkovich space.
    - ``False`` otherwise.
    """
    return isinstance(space, Berkovich)

def is_Berkovich_Cp(space):
    """
    Checks if ``space`` is a Berkovich space over ``Cp``.

    OUTPUT:

    - ``True`` if ``space`` is a Berkovich space over ``Cp``.
    - ``False`` otherwise.
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
            sage: A.<a> = NumberField(x^3+20)
            sage: ideal = A.ideal(-1/2*a^2 + a - 3)
            sage: B = Berkovich_Cp_Affine(A, ideal)
            sage: B.residue_characteristic()
            7
        """
        return self._p

    prime = residue_characteristic

    def is_padic_base(self):
        """
        Returns ``True`` if this Berkovich space is backed by a p-adic field.

        OUTPUT:

        - ``True`` if this Berkovich space was created with a p-adic field
        - ``False`` otherwise

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
        Returns ``True`` if this Berkovich space is backed by a p-adic field.

        OUTPUT:

        - ``True`` if this Berkovich space was created with a number field
        - ``False`` otherwise

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
            sage: A.<a> = NumberField(z^2+1)
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
            sage: print(B.ideal())
            None
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
            sage: A.<a> = NumberField(x^2+1)
            sage: A_ideal = A.prime_above(2)
            sage: B.<b> = NumberField(x^4+1)
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
        """
        if not isinstance(right, Berkovich_Cp):
            return False
        if self._base_type != right._base_type:
            return False
        if self._base_type == 'padic field':
            return self.prime() == right.prime()
        else:
            return self.base() == right.base()

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

class Berkovich_Cp_Affine(Berkovich_Cp):
    r"""
    The Berkovich Affine line over `\CC_p`.
    
    The Berkovich Affine line is the set of seminorms on `\CC_p[x]`,
    with the weakest topology such that the map `| \cdot | \to |f|` is continuous
    for all `f \in \CC_p[x]`.

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

    Initializing by passing in ``Qp`` looks the same::

        sage: B = Berkovich_Cp_Affine(Qp(3)); B
        Affine Berkovich line over Cp(3) of precision 20

    However, this method allows for more control over behind-the-scenes conversion::

        sage: B = Berkovich_Cp_Affine(Qp(3, 1)); B
        Affine Berkovich line over Cp(3) of precision 1

        sage: Q1 = B(1/2); Q1
        Type I point centered at 2 + O(3)

    Note that this point has very low precision, as B was initialized
    with a padic field of capped-relative precision one. For high precision,
    pass in a high precision padic field::

        sage: B = Berkovich_Cp_Affine(Qp(3, 1000)); B
        Affine Berkovich line over Cp(3) of precision 1000

    For exact computation, a number field can be used::

        sage: R.<x> = QQ[]
        sage: A.<a> = NumberField(x^3 + 20)
        sage: ideal = A.prime_above(3)
        sage: B = Berkovich_Cp_Affine(A, ideal); B
        Affine Berkovich line over Cp(3), with base Number
        Field in a with defining polynomial x^3 + 20
    """

    Element = Berkovich_Element_Cp_Affine

    def __init__(self, base, ideal=None):
        if base in ZZ:
            if base.is_prime():
                base = Qp(base) #TODO chance to Qpbar
            else:
                raise ValueError("non-prime pased into Berkovich space")
        if base in NumberFields():
            if ideal == None:
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
        elif is_pAdicField(base): #TODO change base to Qpbar(prime)
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
            sage: B = Berkovich_Cp_Affine(A, ideal); B
            Affine Berkovich line over Cp(3), with base Number Field
            in a with defining polynomial z^2 + 1
        """
        if self._base_type == 'padic field':
            return "Affine Berkovich line over Cp(%s) of precision %s" %(self.prime(),\
                self.base().precision_cap())
        else:
            return  "Affine Berkovich line over Cp(%s), with base %s" %(self.prime(),\
                self.base())

    def _projective_space(self):
        """
        Creates a projective Berkovich space with the same characteristics as this space.

        OUTPUT: A projective Berkovich space over ``Cp``.

        EXAMPLES::

            sage: B = Berkovich_Cp_Affine(3)
            sage: B._projective_space()
            Projective Berkovich line over Cp(3) of precision 20

        ::

            sage: R.<z> = QQ[]
            sage: A.<a> = NumberField(z^3 + 20)
            sage: ideal = A.prime_above(3)
            sage: B = Berkovich_Cp_Affine(A, ideal)
            sage: B._projective_space()
            Projective Berkovich line over Cp(3), with base Number Field
            in a with defining polynomial z^3 + 20
        """
        return Berkovich_Cp_Projective(self.base_ring(), self.ideal())

    def _latex_(self):
        r"""
        LaTeX representation of this Berkovich Space.

        EXAMPLES:

            sage: B = Berkovich_Cp_Affine(3)
            sage: latex(B)
            \text{Affine Berkovich line over } \Bold{C}_{3}
        """
        return r"\text{Affine Berkovich line over } \Bold{C}_{%s}" %(self.prime())

class Berkovich_Cp_Projective(Berkovich_Cp):
    r"""
    The Berkovich projective line over `\CC_p`.

    The Berkovich projective line is the one-point compactification
    of the Berkovich Affine line.

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

    Initializing by passing in a padic space looks the same::

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
        sage: A.<a> = NumberField(x^2+1)
        sage: ideal = A.prime_above(2)
        sage: B = Berkovich_Cp_Projective(A, ideal); B
        Projective Berkovich line over Cp(2), with base
        Number Field in a with defining polynomial x^2 + 1

    Number fields have the benefit that computation is exact,
    but lack support for all of `\CC_p`.
    """

    Element = Berkovich_Element_Cp_Projective

    def __init__(self, base, ideal=None):
        if base in ZZ:
            if base.is_prime():
                base = ProjectiveSpace(Qp(base), 1)
            else:
                raise ValueError("non-prime pased into Berkovich space")
        if base in NumberFields() or is_pAdicField(base):
            base = ProjectiveSpace(base, 1)
        if not is_ProjectiveSpace(base):
            raise ValueError("base of projective Berkovich space must be projective space")
        if not (is_pAdicField(base.base_ring())):
            if base.base_ring() not in NumberFields():
                raise ValueError("base of projective Berkovich space must be " + \
                    "projective space over Qp or a number field")
            else:
                if ideal == None:
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
            sage: A.<a> = NumberField(x^3+20)
            sage: ideal = A.prime_above(3)
            sage: D = Berkovich_Cp_Projective(A, ideal)
            sage: D.base_ring()
            Number Field in a with defining polynomial x^3 + 20
        """
        return self.base().base_ring()

    def _affine_space(self):
        """
        Creates an affine Berkovich space with the same characteristics as this space.

        OUTPUT: An affine Berkovich space over ``Cp``.

        EXAMPLES::

            sage: B = Berkovich_Cp_Projective(3)
            sage: B._affine_space()
            Affine Berkovich line over Cp(3) of precision 20

        ::

            sage: R.<z> = QQ[]
            sage: A.<a> = NumberField(z^3+20)
            sage: ideal = A.prime_above(3)
            sage: B = Berkovich_Cp_Projective(A, ideal)
            sage: B._affine_space()
            Affine Berkovich line over Cp(3), with base Number Field in
            a with defining polynomial z^3 + 20
        """
        return Berkovich_Cp_Affine(self.base_ring(), self.ideal())

    def _repr_(self):
        """
        String representation of this Berkovich Space.

        EXAMPLES::

            sage: B = Berkovich_Cp_Projective(3)
            sage: B
            Projective Berkovich line over Cp(3) of precision 20

        ::

            sage: R.<x> = QQ[]
            sage: A.<a> = NumberField(x^2+1)
            sage: v = A.ideal(a+1)
            sage: B = Berkovich_Cp_Projective(A, v); B
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