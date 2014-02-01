"""
Lie Algebra Elements

AUTHORS:

- Travis Scrimshaw (2005-05-04): Initial implementation
"""

#*****************************************************************************
#  Copyright (C) 2013 Travis Scrimshaw <tscrim@ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty
#    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
#  See the GNU General Public License for more details; the full text
#  is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.abstract_method import abstract_method
from sage.misc.classcall_metaclass import ClasscallMetaclass, typecall
from sage.misc.misc import repr_lincomb
from copy import copy
from sage.structure.element import RingElement, coerce_binop
from sage.structure.sage_object import SageObject
from sage.combinat.free_module import CombinatorialFreeModuleElement

class LieGenerator(SageObject):
    """
    A wrapper around a string so it can compare with :class:`LieBracket`.
    """
    __metaclass__ = ClasscallMetaclass

    @staticmethod
    def __classcall_private__(cls, name):
        """
        Return ``name`` if it is a :class:`LieGenerator`, otherwise construct
        a new object.

        EXAMPLES::
        """
        if isinstance(name, LieGenerator):
            return name
        return typecall(cls, name)

    def __init__(self, name):
        """
        Initalize ``self``.

        EXAMPLES::
        """
        self._name = name

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::
        """
        return self._name

    def __eq__(self, rhs):
        """
        Compare equals.

        EXAMPLES::
        """
        return isinstance(rhs, LieGenerator) and self._name == rhs._name

    def __lt__(self, rhs):
        """
        Compare less than.

        EXAMPLES::
        """
        if isinstance(rhs, LieGenerator):
            return self._name < rhs._name
        if isinstance(rhs, LieBracket):
            return not rhs.__lt__(self) # Clearly self != rhs
        return False

    def __le__(self, rhs):
        """
        Compare less than or equals.

        EXAMPLES::
        """
        return self.__lt__(rhs) or self.__eq__(rhs)

    def _im_gens_(self, codomain, im_gens, names):
        """
        Return the image of ``self``.

        EXAMPLES::
        """
        x = im_gens[names.index(self._name)]
        return im_gens[names.index(self._name)]

    def to_word(self):
        """
        Return ``self`` as a word in the variable names.

        EXAMPLES::
        """
        return [self._name]

class LieBracket(SageObject):
    """
    A Lie bracket. This is the building blocks for Lie algebra elements.
    """
    def __init__(self, l, r):
        """
        Initialize ``self``.

        EXAMPLES::
        """
        self._left = l
        self._right = r

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::
        """
        return "[%s, %s]"%(self._left, self._right)

    def _latex_(self):
        r"""
        Return a `\LaTeX` representation of ``self``.

        EXAMPLES::
        """
        from sage.misc.latex import latex
        return "\\left[%s,%s\\right]"%(latex(self._left), latex(self._right))

    def __getitem__(self, i):
        r"""
        Return the `i`-th item of ``self``.

        EXAMPLES::
        """
        if i == 0:
            return self._left
        if i == 1:
            return self._right
        raise IndexError("i must be either 0 or 1")

    def __eq__(self, rhs):
        """
        Check equality.

        EXAMPLES::
        """
        if not isinstance(rhs, LieBracket):
            return False
        return self._left == rhs._left and self._right == rhs._right

    def __ne__(self, rhs):
        """
        Check inequality.

        EXAMPLES::
        """
        return not self.__eq__(rhs)

    def __lt__(self, rhs):
        """
        Check less than.

        EXAMPLES::
        """
        if not isinstance(rhs, LieBracket):
            return False
        if self._left < rhs._left:
            return True
        elif self._left == rhs._left:
            return self._right < rhs._right
        return False

    def __le__(self, rhs):
        """
        Check less than or equal.

        EXAMPLES::
        """
        return self.__lt__(rhs) or self.__eq__(rhs)

    def __gt__(self, rhs):
        """
        Check greater than.

        EXAMPLES::
        """
        return rhs.__lt__(self)

    def __ge__(self, rhs):
        """
        Check greater than or equal.

        EXAMPLES::
        """
        return self.__gt__(rhs) or self.__eq__(rhs)

    def __hash__(self):
        """
        Return the hash value of ``self``.

        EXAMPLES::
        """
        return hash((self._left, self._right))

    def _im_gens_(self, codomain, im_gens, names):
        """
        Return the image of ``self``.

        EXAMPLES::
        """
        return codomain.bracket(self._left._im_gens_(codomain, im_gens, names),
                                self._right._im_gens_(codomain, im_gens, names))

    def lift(self, UEA_gens_dict):
        """
        Lift ``self`` to the UEA.

        EXAMPLES::
        """
        if isinstance(self._left, LieBracket):
            l = self._left.lift(UEA_gens_dict)
        else:
            l = UEA._gens_dict[self._left]

        if isinstance(self._right, LieBracket):
            r = self._right.lift(UEA_gens_dict)
        else:
            r = UEA_gens_dict[self._right]

        return l*r - r*l

    def to_word(self):
        """
        Return ``self`` as a word expressed in the variable names.

        EXAMPLES::
        """
        return self._left.to_word() + self._right.to_word()

class GradedLieBracket(LieBracket):
    """
    A Lie bracket in a graded Lie algebra.
    """
    def __init__(self, l, r, grade):
        """
        Initialize ``self``.

        EXAMPLES::
        """
        self._grade = grade
        LieBracket.__init__(self, l, r)

    def __lt__(self, rhs):
        """
        Check less than.

        EXAMPLES::
        """
        if isinstance(rhs, GradedLieBracket) and self._grade != rhs._grade:
            return self._grade < rhs._grade
        if isinstance(rhs, LieGenerator):
            return False # Lie generators have grade 1 and our grade > 1
        return LieBracket.__lt__(self, rhs)

    def __hash__(self):
        """
        Return the hash value of ``self``.

        EXAMPLES::
        """
        return hash((self._grade, self._left, self._right))

# TODO: factor out parts of CombinatorialFreeModuleElement into a SparseFreeModuleElement
class LieAlgebraElement(CombinatorialFreeModuleElement):
    """
    A Lie algebra element.
    """
    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::
        """
        if len(self._monomial_coefficients) == 0:
            return '0'
        return repr_lincomb(self.list())

    # Default implementation
    def _latex_monomial(self, m):
        """
        Return a `\LaTeX` representation of the monomial ``m``.

        EXAMPLES::
        """
        from sage.misc.latex import latex
        return latex(m)

    def _latex_(self):
        r"""
        Return a `\LaTeX` representation of ``self``.

        EXAMPLES::
        """
        if len(self._monomial_coefficients) == 0:
            return '0'
        return repr_lincomb(self.list(), repr_monomial=self._latex_monomial, is_latex=True)

    def __mul__(self, y):
        """
        If we are multiplying two non-zero elements, automatically
        lift up to the universal enveloping algebra.

        EXAMPLES::
        """
        if self == 0 or y == 0:
            return self.parent().zero()
        # Otherwise we lift to the UEA
        return self.lift() * y

    def _im_gens_(self, codomain, im_gens):
        """
        Return the image of ``self`` in ``codomain`` under the map that sends
        the images of the generators of the parent of ``self`` to the
        tuple of elements of ``im_gens``.

        EXAMPLES::
        """
        s = codomain.zero()
        if not self: # If we are 0
            return s
        names = self.parent().variable_names()
        return codomain.sum(c * t._im_gens_(codomain, im_gens, names)
                            for t, c in self._monomial_coefficients.iteritems())

    def lift(self):
        """
        Lift ``self`` to the universal enveloping algebra.

        EXAMPLES::
        """
        UEA = self.parent().universal_enveloping_algebra()
        gen_dict = UEA.gens_dict()
        s = UEA.zero()
        if not self:
            return s
        for t, c in self._monomial_coefficients.iteritems():
            if isinstance(t, LieBracket):
                s += c * t.lift(gen_dict)
            else:
                s += c * gen_dict[t._name]
        return s

    def is_constant(self):
        """
        Check if ``self`` is a constant (i.e. zero).

        EXAMPLES::
        """
        return len(self._monomial_coefficients) == 0

    def dict(self):
        """
        Return ``self`` as a dictionary mapping monomials to coefficients.

        EXAMPLES::
        """
        return copy(self._monomial_coefficients)

    def list(self):
        """
        Return ``self`` as a list of pairs ``(m, c)`` where ``m`` is a
        monomial and ``c`` is the coefficient.

        EXAMPLES::
        """
        L = self._monomial_coefficients.items()
        L.sort()
        return L

