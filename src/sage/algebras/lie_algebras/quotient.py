"""
Quotient Lie Algebras

AUTHORS:

- Travis Scrimshaw (2013-05-03): Initial version
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

from sage.misc.cachefunc import cached_method
from sage.algebras.lie_algebras.lie_algebra import FinitelyGeneratedLieAlgebra
#from sage.algebras.lie_algebras.lie_algebra_element import LieAlgebraElement
from sage.structure.category_object import check_default_category
from sage.structure.element_wrapper import ElementWrapper
from sage.categories.lie_algebras import LieAlgebras

class QuotientLieAlgebraElement(ElementWrapper):
    r"""
    An element of a quotient Lie algebra.
    """
    def __init__(self, parent, rep, reduce=True):
        """
        An element of a quotient ring `R/I`.  See
        ``QuotientRingElement`` for full documentation.

        EXAMPLES::
        
            sage: R.<x> = PolynomialRing(ZZ)
            sage: S.<xbar> = R.quo((4 + 3*x + x^2, 1 + x^2)); S
            Quotient of Univariate Polynomial Ring in x over Integer Ring by the ideal (x^2 + 3*x + 4, x^2 + 1)
            sage: v = S.gens(); v
            (xbar,)
        """
        ElementWrapper.__init__(self, parent, rep)
        if reduce:
            self._reduce_()

    def _repr_(self):
        """
        Return a string representation of ``self``.
        """
        # We assume "good" variable names, nothing like 'x' and 'xx'
        #   so we can do a straight string replacement.
        # This is a hack until a better idea comes along...
        orig_names = self.parent().cover_lie_algebra().variable_names()
        new_names = self.parent().variable_names()
        ret = repr(self.value)
        for i,n in enumerate(new_names):
            ret = ret.replace(orig_names[i], n)
        return ret

    def __hash__(self):
        """
        Return the hash value of ``self``.
        """
        # Perform the reduction to guarantee the hash value is constant
        self._reduce_()
        return hash(self.value)

    def _reduce_(self):
        """
        Reduce the element modulo the defining ideal of the quotient
        Lie algebra.  This internal method replaces the cached representative
        by one in reduced form.

        (Note that this has nothing to do with pickling.)

        TESTS::
        """
        I = self.parent().defining_ideal()
        self.value = I.reduce(self.value)

    def _add_(self, rhs):
        """
        Add ``self`` to ``rhs``.
        """
        return self.__class__(self.parent(), self.value + rhs.value)

    def __neg__(self):
        """
        Return the negation of ``self``.
        """
        return self.__class__(self.parent(), -self.value)

    def _bracket_(self, rhs):
        """
        Return the bracket of ``self`` with ``rhs``.
        """
        return self.__class__(self.parent(), self.value._bracket_(rhs.value))

class QuotientLieAlgebra(FinitelyGeneratedLieAlgebra):
    r"""
    The quotient Lie algebra of `\mathfrak{g}` by an ideal `I`.

    INPUT:
    
    - ``lie`` -- a Lie algebra

    - ``I`` -- an ideal of ``lie``

    - ``names`` -- a list of generator names
    """
    def __init__(self, lie, I, names=None, index_set=None, category=None):
        r"""
        Create the quotient Lie algebra of `\mathfrak{g}` by the ideal `I`.

        EXAMPLES::
        """
        self.__lie = lie
        self.__I = I
        R = lie.base_ring()
        ##
        # Unfortunately, computing the join of categories, which is done in
        # check_default_category, is very expensive.
        # However, we don't just want to use the given category without mixing in
        # some quotient stuff - unless Parent.__init__ was called
        # previously, in which case the quotient ring stuff is just
        # a vaste of time. This is the case for FiniteField_prime_modn.
        if not self._is_category_initialized():
            if category is None:
                category = check_default_category(LieAlgebras(R).Quotients(), category)
            FinitelyGeneratedLieAlgebra.__init__(self, R, names, index_set, category)
        # self._populate_coercion_lists_([lie]) # we don't want to do this, since subclasses will often implement improved coercion maps.

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::
        """
        return "Quotient of {0} by the ideal {1}".format(self.cover_lie_algebra(), self.defining_ideal()._repr_short())

    def _element_constructor_(self, x, coerce=True):
        """
        Construct an element with ``self`` as the parent.

        EXAMPLES::
        """
        if isinstance(x, list):
            return FinitelyGeneratedLieAlgebra._element_constructor_(self, x)

        if isinstance(x, QuotientLieAlgebraElement):
            if x.parent() is self:
                return x
            x = x.lift()
        if coerce:
            g = self.cover_lie_algebra()
            x = g(x)
        return self.element_class(self, x)

    Element = QuotientLieAlgebraElement

    def gen(self, i=0):
        """
        Return the ``i``-th generator of ``self``.
        """
        return self.element_class(self, self.__lie.gen(i))

    @cached_method
    def zero(self):
        """
        Return the element `0`.
        """
        return self.element_class(self, self.__lie.zero())

    zero_element = zero

    def defining_ideal(self):
        r"""
        Return the ideal generating this quotient Lie algebra.

        EXAMPLES::
        """
        return self.__I

    def cover_lie_algebra(self):
        r"""
        Return the cover Lie algebra of the quotient Lie algebra: that is,
        the original Lie algebra `\mathfrak{g}` from which we modded out
        an ideal, `I`.

        EXAMPLES::
        """
        return self.__lie

    def _construct_UEA(self):
        """
        Construct the universal enveloping algebra of ``self``.

        EXAMPLES::
        """
        UEA = self.__lie.UEA()
        I = UEA.ideal([UEA(x) for x in self.__I.gens()])
        return UEA.quotient(I)

    # This is to make the category framework happy
    ambient = cover_lie_algebra

