"""
Hecke algebras for Weyl groups with support for multiple parameters

AUTHOR:

- Mark Shimozono 

"""
#*****************************************************************************
#  Copyright (C) 2014 Mark Shimozono <mshimo at math.vt.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.structure.element import parent
from sage.categories.all import AlgebrasWithBasis, WeylGroups
from sage.combinat.root_system.weyl_group import WeylGroup
from sage.combinat.family import Family, FiniteFamily
from sage.combinat.free_module import CombinatorialFreeModule, CombinatorialFreeModuleElement
from sage.combinat.root_system.cartan_type import CartanType
from sage.rings.all import QQ, ZZ

class MultiParameterHeckeAlgebraElement(CombinatorialFreeModuleElement):
    r"""
    The element class of :class:`MultiParameterHeckeAlgebra`.

    EXAMPLES::

        sage: H = MultiParameterHeckeAlgebra("B3")
        sage: T1,T2,T3 = H.algebra_generators()
        sage: T1+2*T2*T3
        2*T[2,3] + T[1]
        sage: T1*T1
        ((v^2-1)/v)*T[1] + 1

        sage: H = MultiParameterHeckeAlgebra("A2",prefix="x")
        sage: sum(H.algebra_generators())^2
        x[1,2] + x[2,1] + ((v^2-1)/v)*x[1] + ((v^2-1)/v)*x[2] + 2

        sage: H = MultiParameterHeckeAlgebra("A2",prefix="t")
        sage: t1,t2 = H.algebra_generators()
        sage: (t1-t2)^3
        ((v^4+v^2+1)/v^2)*t[1] + ((-v^4-v^2-1)/v^2)*t[2]

        sage: H = MultiParameterHeckeAlgebra("G2")
        sage: T1, T2 = H.algebra_generators()
        sage: T1*T2*T1*T2*T1*T2 == T2*T1*T2*T1*T2*T1
        True
        sage: T1*T2*T1 == T2*T1*T2
        False

    """

    def _repr_term(self, t):
        r"""
        Return the string representation of the term indexed by ``t``.

        EXAMPLES::

            sage: H = MultiParameterHeckeAlgebra("A3")
            sage: W = H.weyl_group()
            sage: H._repr_term(W.from_reduced_word([1,2,3]))
            'T[1,2,3]'
        """
        redword = t.reduced_word()
        if len(redword) == 0:
            return "1"
        return self._print_options['prefix'] + '[%s]'%','.join('%d'%i for i in redword)

    def _latex_term(self, t):
        r"""
        Return latex for the term indexed by ``t``.

        EXAMPLES::

            sage: H = MultiParameterHeckeAlgebra("A3")
            sage: W = H.weyl_group()
            sage: H._latex_term(W.from_reduced_word([1,2,3]))
            'T_{1}T_{2}T_{3}'
        """
        redword = t.reduced_word()
        if len(redword) == 0:
            return '1'
        return ''.join("%s_{%d}"%(self._print_options['prefix'], i) for i in redword)

    def inverse(self):
        r"""
        Return the inverse if ``self`` is a basis element.

        An element is a basis element if it is `T_w` where `w` is in
        the Weyl group. The base ring must be a field or Laurent
        polynomial ring. Other elements of the ring have inverses but
        the inverse method is only implemented for the basis elements.

        EXAMPLES::

            sage: H = MultiParameterHeckeAlgebra("A2")
            sage: T1,T2 = H.algebra_generators()
            sage: x = (T1*T2).inverse(); x
            T[2,1] + ((-v^2+1)/v)*T[1] + ((-v^2+1)/v)*T[2] + ((v^4-2*v^2+1)/v^2)
            sage: x*T1*T2
            1

        TESTS:

        We check some alternative forms of input for inverting
        an element::

            sage: H = MultiParameterHeckeAlgebra("A2")
            sage: T1,T2 = H.algebra_generators()
            sage: ~(T1*T2)
            T[2,1] + ((-v^2+1)/v)*T[1] + ((-v^2+1)/v)*T[2] + ((v^4-2*v^2+1)/v^2)
            sage: (T1*T2)^(-1)
            T[2,1] + ((-v^2+1)/v)*T[1] + ((-v^2+1)/v)*T[2] + ((v^4-2*v^2+1)/v^2)
        """

        if len(self.support()) != 1:
            raise NotImplementedError("inverse only implemented for basis elements (monomials in the generators)"%self)
        H = self.parent()
        w = self.support_of_term()

        return H.prod(H.inverse_generator(i) for i in reversed(w.reduced_word()))

    __invert__ = inverse


class MultiParameterHeckeAlgebra(CombinatorialFreeModule):
    r"""
    Hecke algebra for a Weyl group with support for multiple parameters.

    INPUT:

    - ``weyl_group`` -- A Weyl group or Cartan type
    - ``q1, q2`` -- (default: None) Parameters
    - ``prefix`` -- (default: None) String for printing basis elements
    - ``category`` -- (default: None) category

    It is assumed that two parameters are equal if the corresponding simple roots are
    conjugate by the Weyl group. 

    - None -- `q1` and `q2` are set to Families with constant value `v` and `-1/v` for `v` the generator of `QQ[v].fraction_field()`.
    - ring elements `q1, q2` -- The result is a pair of constant families with values `q1, q2`.
    - ring element `q1` -- Same as with two ring elements except the second is set to `-1/q1`.
    - Family `q1` -- In this case `q2` is the Family such that `q2[i]` has value `-1/q1[i]` for all `i`.
    - Families `q1, q2` -- No adjustment needed.

    EXAMPLES::

        sage: H = MultiParameterHeckeAlgebra("A2"); H
        Hecke algebra of type ['A', 2]
        sage: K = QQ['v,vl'].fraction_field()
        sage: v,vl=K.gens()
        sage: H = MultiParameterHeckeAlgebra("B2", Family(dict({1:vl,2:v}))); H
        Hecke algebra of type ['B', 2]
        sage: L = QQ['v,vl,v0'].fraction_field()
        sage: v,vl,v0=L.gens()
        sage: H = MultiParameterHeckeAlgebra(['D',3,2],Family(dict({0:v0,1:vl,2:v}))); H
        Hecke algebra of type ['C', 2, 1]^*

    """

    @staticmethod
    def __classcall_private__(cls, weyl_group, q1=None, q2=None, prefix=None, category=None):
        if not weyl_group in WeylGroups():
            cartan_type = CartanType(weyl_group)
            weyl_group=WeylGroup(cartan_type,prefix="s")
        if isinstance(q1, dict):
            q1 = Family(q1)
        if isinstance(q2, dict):
            q2 = Family(q2)
        return super(MultiParameterHeckeAlgebra, cls).__classcall__(cls, weyl_group, q1, q2, prefix, category)

    def __init__(self, weyl_group, q1, q2, prefix, category):
        r"""
        Initialize Hecke algebra

        INPUT:

        - ``weyl_group`` -- Either a Weyl group or a Cartan type
        - ``q1,q2`` -- parameters
        - ``prefix`` -- (default: None) A prefix for the basis elements
        - ``category`` -- (default: None) A category for the algebra

        Each of `q1` and `q2` can be None, or a single ring element, or a family of
        ring elements with keys the Dynkin node set `I`. They are used to create
        a pair of families `f1` and `f2` from `I` to a ring `R`.

        If `q1` is a dictionary/Family then `f1` is set to `q1`.
        If `q1` is a single variable then `f1` is the constant Family on `I`
        with all values equal to `q1`. If `q1` is None then a ring `QQ['v'].fraction_field()`
        is created and `f1` is the constant Family on `I` with value `v`.
        If `q2` is a Family or single variable then `f2` is created from `q2`
        in the same way `f1` is formed from `q1`. If `q2` is None then `f2` is the family
        on `I` defined by `f2[i] = -1/f1[i]` for all `i\in I`.
        """
        if not weyl_group in WeylGroups():
            raise TypeError, "%s is not a Weyl group"%weyl_group
        self._weyl_group = weyl_group
        I = self._weyl_group.index_set()
        if q1 and isinstance(q1, FiniteFamily):
            self._q1 = q1
            try:
                # the 0 is sleazy
                F = q1[I[0]].parent()
            except AttributeError:
                raise TypeError, "Values in %s should be ring elements"%q1
        else:
            if q1:
                try:
                    F = q1.parent()
                except AttributeError:
                    raise TypeError, "%s does not have a parent"%q1
            else:
                F = QQ['v'].fraction_field()
                q1 = F.gen(0)
            self._q1 = Family(dict([[i,q1] for i in I]))
        self._base_ring = F
        if q2 and isinstance(q2, FiniteFamily):
            self._q2 = q2
        else:
            if q2:
                try:
                    F = q2.parent()
                except AttributeError:
                    raise TypeError, "%s should be a ring element"%q2
                if F != self._base_ring:
                    raise TypeError, "%s should have parent"%(q1,self._base_ring)
                self._q2 = Family(dict([[i,q2] for i in I]))
            else:
                self._q2 = Family(dict([[i,-1/self._q1[i]] for i in I]))

        # Used when multiplying generators: minor speed-up as it avoids the
        # need to constantly add and multiply the parameters when applying the
        # quadratic relation

        self._eigensum = Family(dict([[i,self._q1[i]+self._q2[i]] for i in I]))
        self._eigenproduct = Family(dict([[i,self._q1[i]*self._q2[i]] for i in I]))

        if not prefix:
            prefix = "T"
        if not category:
            category=AlgebrasWithBasis(F)
        from sage.algebras.iwahori_hecke_algebra import index_cmp
        CombinatorialFreeModule.__init__(self, F, weyl_group, category=category, monomial_cmp = index_cmp, prefix=prefix)

    def _repr_(self):
        r"""
        EXAMPLES::

            sage: MultiParameterHeckeAlgebra(['A',2,1])
            Hecke algebra of type ['A', 2, 1]

        """
        return "Hecke algebra of type %s"%self.cartan_type()

    @cached_method
    def one_basis(self):
        return self.weyl_group().one()

    @cached_method
    def cartan_type(self):
        r"""
        Return the Cartan type of ``self``.

        EXAMPLES::

            sage: MultiParameterHeckeAlgebra("B3").cartan_type()
            ['B', 3]

        """
        return self._weyl_group.cartan_type()

    def weyl_group(self):
        r"""
        Return the Weyl group of ``self``.

        EXAMPLES::

            sage: F = QQ['v'].fraction_field()
            sage: v = F.gen(0)
            sage: MultiParameterHeckeAlgebra("B2",Family(dict({1:v,2:v}))).weyl_group()
            Weyl Group of type ['B', 2] (as a matrix group acting on the ambient space)

        """
        return self._weyl_group

    @cached_method
    def index_set(self):
        r"""
        Return the index set of ``self``.

        EXAMPLES::

            sage: MultiParameterHeckeAlgebra(['A',2,1]).index_set()
            (0, 1, 2)

        """
        return self.weyl_group().index_set()

    @cached_method
    def algebra_generators(self):
        r"""
        Return the algebra generators.

        EXAMPLES::

            sage: MultiParameterHeckeAlgebra(['A',2,1]).algebra_generators()
            Finite family {0: T[0], 1: T[1], 2: T[2]}

        """
        return Family(dict([[i,self.monomial(self.weyl_group().simple_reflection(i))] for i in self.index_set()]))

    def algebra_generator(self, i):
        r"""
        Return the `i`-th generator of ``self``.

        EXAMPLES::

            sage: H = MultiParameterHeckeAlgebra("B2").algebra_generator(1)

        """
        return self.algebra_generators()[i]


    def _repr_term(self, t):
        r"""
        Return the string representation of the term indexed by ``t``.

        EXAMPLES::

            sage: H = MultiParameterHeckeAlgebra("A3")
            sage: W = H.weyl_group()
            sage: H._repr_term(W.from_reduced_word([1,2,3]))
            'T[1,2,3]'
        """
        redword = t.reduced_word()
        if len(redword) == 0:
            return "1"
        return self._print_options['prefix'] + '[%s]'%','.join('%d'%i for i in redword)

    def _latex_term(self, t):
        r"""
        Return latex for the term indexed by ``t``.

        EXAMPLES::

            sage: H = MultiParameterHeckeAlgebra("A3")
            sage: W = H.weyl_group()
            sage: H._latex_term(W.from_reduced_word([1,2,3]))
            'T_{1}T_{2}T_{3}'

        """
        redword = t.reduced_word()
        if len(redword) == 0:
            return '1'
        return ''.join("%s_{%d}"%(self._print_options['prefix'], i) for i in redword)

    def inverse_generator(self, i):
        r"""
        Return the inverse of the `i`-th generator.

        EXAMPLES::

            sage: H = MultiParameterHeckeAlgebra("A2")
            sage: H.base_ring()
            Fraction Field of Univariate Polynomial Ring in v over Rational Field
            sage: H.inverse_generator(1)
            T[1] + ((-v^2+1)/v)
            sage: H.inverse_generator(2)
            T[2] + ((-v^2+1)/v)
            sage: K = QQ['v,vl'].fraction_field()
            sage: v,vl=K.gens()
            sage: H1 = MultiParameterHeckeAlgebra("B2", Family(dict({1:vl,2:v})))
            sage: H1.base_ring()
            Fraction Field of Multivariate Polynomial Ring in v, vl over Rational Field
            sage: H1.inverse_generator(2)
            T[2] + ((-v^2+1)/v)           
            sage: H2 = MultiParameterHeckeAlgebra("C2", Family(dict({1:v,2:vl})))
            sage: H2.inverse_generator(2)
            T[2] + ((-vl^2+1)/vl)           

        """
        return (1/-self._eigenproduct[i])*self.algebra_generator(i) +  (-self._eigensum[i]/-self._eigenproduct[i])*self.one()

    @cached_method
    def inverse_generators(self):
        r"""
        Return the inverses of all the generators, if they exist.

        EXAMPLES::

            sage: K = QQ['q1,q2'].fraction_field()
            sage: q1,q2=K.gens()
            sage: H = MultiParameterHeckeAlgebra("A2",q1,q2)
            sage: T1,T2 = H.algebra_generators()
            sage: U1,U2 = H.inverse_generators()
            sage: U1*T1,T1*U1
            (1, 1)
            sage: H1 = MultiParameterHeckeAlgebra("A2",prefix="V")
            sage: V1,V2 = H1.algebra_generators()
            sage: W1,W2 = H1.inverse_generators()
            sage: [W1,W2]
            [V[1] + ((-v^2+1)/v), V[2] + ((-v^2+1)/v)]
            sage: V1*W1, W2*V2
            (1, 1)
        """
        return Family(self.index_set(), self.inverse_generator)

    def product_on_basis(self, w1, w2):
        r"""
        Return `T_{w_1} T_{w_2}`, where `w_1` and `w_2` are in the  Weyl group.

        EXAMPLES::

            sage: H = MultiParameterHeckeAlgebra("A2")
            sage: s1,s2 = H.weyl_group().simple_reflections()
            sage: [H.product_on_basis(s1,x) for x in [s1,s2]]
            [((v^2-1)/v)*T[1] + 1, T[1,2]]
        """
        result = self.monomial(w1)
        for i in w2.reduced_word():
            result = self.product_by_generator(result, i)
        return result

    def product_by_generator_on_basis(self, w, i, side="right"):
        r"""
        Return the product `T_w T_i` (resp. `T_i T_w`) if ``side`` is
        ``'right'`` (resp. ``'left'``).

        INPUT:

        - ``w`` -- an element of the Weyl group
        - ``i`` -- an element of the index set
        - ``side`` -- ``'right'`` (default) or ``'left'``

        EXAMPLES::

            sage: K = QQ['q1,q2'].fraction_field()
            sage: q1,q2=K.gens()
            sage: H = MultiParameterHeckeAlgebra("A2",q1,q2)
            sage: s1,s2 = H.weyl_group().simple_reflections()
            sage: [H.product_by_generator_on_basis(w, 1) for w in [s1,s2,s1*s2]]
            [(q1+q2)*T[1] + (-q1*q2), T[2,1], T[1,2,1]]
            sage: [H.product_by_generator_on_basis(w, 1, side="left") for w in [s1,s2,s1*s2]]
            [(q1+q2)*T[1] + (-q1*q2), T[1,2], (q1+q2)*T[1,2] + (-q1*q2)*T[2]]
        """
        wi = w.apply_simple_reflection(i, side = side)
        if w.has_descent(i, side = side):
            # 10% faster than a plain addition on the example of #12528
            return self.sum_of_terms(((w , self._eigensum[i]), (wi, -self._eigenproduct[i])), distinct=True)
        else:
            return self.monomial(wi)

    def product_by_generator(self, x, i, side="right"):
        r"""
        Return `T_i \cdot x`, where `T_i` is the `i`-th generator. This is
        coded individually for use in ``x._mul_()``.

        EXAMPLES::

            sage: K = QQ['q1,q2'].fraction_field()
            sage: q1,q2=K.gens()
            sage: H = MultiParameterHeckeAlgebra("A2",q1,q2)
            sage: T1, T2 = H.algebra_generators()
            sage: [H.product_by_generator(x, 1) for x in [T1,T2]]
            [(q1+q2)*T[1] + (-q1*q2), T[2,1]]
            sage: [H.product_by_generator(x, 1, side = "left") for x in [T1,T2]]
            [(q1+q2)*T[1] + (-q1*q2), T[1,2]]

        """
        return self.linear_combination((self.product_by_generator_on_basis(w, i, side), c)
                                       for (w,c) in x)

    def product_by_inverse_generator_on_basis(self, w, i, side="right"):
        r"""
        Return the product `T_w T_i^{-1}` (resp. `T_i^{-1} T_w`) if ``side`` is
        ``'right'`` (resp. ``'left'``).

        INPUT:

        - ``w`` -- an element of the Weyl group
        - ``i`` -- an element of the index set
        - ``side`` -- ``'right'`` (default) or ``'left'``

        EXAMPLES::

            sage: H = MultiParameterHeckeAlgebra("A2")
            sage: s1,s2 = H.weyl_group().simple_reflections()
            sage: [H.product_by_inverse_generator_on_basis(w, 1) for w in [s1,s2,s1*s2]]
            [1, T[2,1] + ((-v^2+1)/v)*T[2], T[1,2,1] + ((-v^2+1)/v)*T[1,2]]
            sage: [H.product_by_inverse_generator_on_basis(w, 1, side="left") for w in [s1,s2,s1*s2]]
            [1, T[1,2] + ((-v^2+1)/v)*T[2], T[2]]

        """
        wi = w.apply_simple_reflection(i, side = side)
        if w.has_descent(i, side = side):
            # 10% faster than a plain addition on the example of #12528
            return self.monomial(wi)
        else:
            return self.sum_of_terms(((w , -self._eigensum[i]), (wi, 1)), distinct=True)

    def product_by_inverse_generator(self, x, i, side="right"):
        r"""
        Return `x \cdot T_i^{-1}` where `T_i` is the `i`-th generator.

        If ``side`` is "left" return `T_i^{-1} \cdot x`.

        EXAMPLES::

            sage: H = MultiParameterHeckeAlgebra("A2")
            sage: T1, T2 = H.algebra_generators()
            sage: [H.product_by_inverse_generator(x, 1) for x in [T1,T2]]
            [1, T[2,1] + ((-v^2+1)/v)*T[2]]            
            sage: [H.product_by_inverse_generator(x, 1, side = "left") for x in [T1,T2]]
            [1, T[1,2] + ((-v^2+1)/v)*T[2]]

        """
        return self.linear_combination((self.product_by_inverse_generator_on_basis(w, i, side), c)
                                       for (w,c) in x)

    def product_by_signed_generator_sequence(self, x, word, signs, side="right"):
        r"""
        Returns `x T_{i_1}^{e_1} T_{i_2}^{e_2} \dotsm T_{i_l}^{e_l}` where the indices `i_k`
        are given by `word` and the signs `e_k = \pm 1` are given by `signs`.

        If ``side`` is "left", returns the product `T_{i_1}^{e_1}\dotsm T_{i_l}^{e_l} x`.

        EXAMPLES::

            sage: H = MultiParameterHeckeAlgebra(CartanType(['A',2,1]))
            sage: y = H.product_by_signed_generator_sequence(H.one(), [0,1,2,1], [1,1,1,1]); y
            T[0,1,2,1]
            sage: z = H.product_by_signed_generator_sequence(H.one(), [1,2,1,0], [-1,-1,-1,-1]); z
            T[1,2,1,0] + ((-v^2+1)/v)*T[1,2,0] + ((-v^2+1)/v)*T[2,1,0] + ((-v^2+1)/v)*T[1,2,1] + ((v^4-2*v^2+1)/v^2)*T[1,0] + ((v^4-2*v^2+1)/v^2)*T[2,0] + ((v^4-2*v^2+1)/v^2)*T[2,1] + ((v^4-2*v^2+1)/v^2)*T[1,2] + ((-v^6+2*v^4-2*v^2+1)/v^3)*T[0] + ((-v^6+3*v^4-3*v^2+1)/v^3)*T[2] + ((-v^6+3*v^4-3*v^2+1)/v^3)*T[1] + ((v^8-3*v^6+4*v^4-3*v^2+1)/v^4)
            sage: y*z
            1
            sage: y == H.product_by_signed_generator_sequence(H.one(), [0,1,2,1], [1,1,1,1], side="left")
            True
            sage: z == H.product_by_signed_generator_sequence(H.one(), [1,2,1,0], [-1,-1,-1,-1], side="left")
            True

        """
        seq = range(len(word))
        if side == "left":
            seq.reverse()

        for i in seq:
            if signs[i] == 1:
                x = self.product_by_generator(x, word[i], side=side)
            else:
                x = self.product_by_inverse_generator(x, word[i], side=side)
        return x

    def __getitem__(self, i):
        r"""
        Return the basis element indexed by ``i``.

        INPUT:

        - ``i`` -- either an element of the Weyl group or a
          reduced word

        .. WARNING::

            If `i`` is not a reduced expression then the basis element
            indexed by the corresponding element of the algebra is
            returned rather than the corresponding product of the
            generators::

                sage: H = MultiParameterHeckeAlgebra("A2")
                sage: H[1,1]
                1

        EXAMPLES::

            sage: H = MultiParameterHeckeAlgebra(['B',2,1],prefix="Tx")
            sage: W = H.weyl_group()
            sage: H[W.one()]
            1
            sage: H[W.simple_reflection(1)]
            Tx[1]
            sage: H[W.from_reduced_word([1,2,1])]
            Tx[1,2,1]
            sage: H[[]]
            1
            sage: H[1]
            Tx[1]
            sage: H[1,2,1]
            Tx[1,2,1]
        """
        W = self.weyl_group()
        if i in ZZ:
            return self(W.simple_reflection(i))
        if i in W:
            return self(i)
        if i == []:
            return self.one()
        return self(W.from_reduced_word(i))

MultiParameterHeckeAlgebra.Element = MultiParameterHeckeAlgebraElement

