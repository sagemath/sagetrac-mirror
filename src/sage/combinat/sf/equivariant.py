r"""
Equivariant symmetric functions
"""
from sage.categories.modules_with_basis import ModulesWithBasis
from sage.categories.morphism import SetMorphism
from sage.categories.homset import Hom
from sage.rings.polynomial.zpoly import FunctionFieldIntegerGenerators
from sage.structure.unique_representation import UniqueRepresentation
from sage.sets.family import Family
from sage.misc.cachefunc import cached_method
from sage.combinat.partition import Partition, Partitions
from sage.matrix.matrix_space import MatrixSpace
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.categories.tensor import tensor
import sfa

class EquivariantSymmetricFunctions(UniqueRepresentation):
    r"""
    A container class for the family of bases for the algebra of equivariant symmetric functions.

    This is analogous to the class that implements `Sym.macdonald()` where `Sym` is a symmetric function algebra.

    The main basis is the equivariant Schur basis, sometimes called double or factorial Schur functions.
    These are to Schur functions what double Schubert polynomials are to Schubert polynomials.

    INPUT:

    - ``Sym`` -- A symmetric function algebra
    - ``equivariant_parameter`` -- a function from the integers to the base ring `R` of `Sym`
    - `shift_auto` -- a function of two variables `(r,z)` where `r` is an integer and `z` an element of `R`
    - ``positive`` -- (default: True) Use the positive integers to index the variables in symmetric functions.
    If False, use nonpositive indices.
    - ``positive_roots`` -- (default: True) In Schubert class localizations, use the convention
    for which Graham positivity means a nonnegative integer polynomial in the simple roots.
    If False, the same holds but for negatives of simple roots.

    The function value ``equivariant_parameter(i)`` supplies an element `a_i` of `R` for every integer `i`.
    The function ``shift_auto`` should compute the image of `z` under the automorphism of `R` that
    sends `a_i` to `a_{i+r}` for all integers `i`. Such methods are built into the class
    :class:`FunctionFieldIntegerGenerators` which implements a polynomial ring in variables indexed
    by integers.

    All usual bases of symmetric functions such as power sums, are interpreted as being in the difference of
    variables ``x_i - a_i`` for positive integers `i` where the `x_i` are the usual symmetric
    function variables, where `i` runs over positive or nonpositive integers depending on the parameter ``positive``.

    EXAMPLES::

            sage: from sage.rings.polynomial.zpoly import PolynomialRingIntegerGenerators
            sage: F = PolynomialRingIntegerGenerators(QQ,('a','b'))
            sage: from sage.combinat.sf.sf import SymmetricFunctions
            sage: Sym = SymmetricFunctions(F)
            sage: Eq = Sym.equivariant(F.igen, F.shift_auto_on_element)

    The basis of equivariant Schur functions is denoted `es`.

            sage: es = Eq.es()
            sage: es[2,1]
            EQs[2, 1]
            sage: s = Sym.s()
            sage: ans = s(es[2]); ans
            (-b_0)*s[1] + s[2]
            sage: s(es[2,1])
            (-a_1*b_0)*s[1] + (-b_0)*s[1, 1] + a_1*s[2] + s[2, 1]

    TODO:

    Elements of this ring can be localized at torus-fixed points, which are bijective
    with partitions.
    """
    def __init__(self, Sym, equivariant_parameter, shift, positive=True, positive_roots=True):
        r"""
        Initialize ``self``.
        """
        self._sym = Sym
        self._equivariant_parameter = equivariant_parameter
        self._shift = shift
        self._positive=positive
        self._positive_roots=positive_roots
        self._s = Sym.s()
        self._p = Sym.p()

    def __repr__(self):
        r"""
        A string representing ``self``.

        EXAMPLES::

            sage: from sage.rings.polynomial.zpoly import PolynomialRingIntegerGenerators
            sage: F = PolynomialRingIntegerGenerators(QQ,('a','b'))
            sage: from sage.combinat.sf.sf import SymmetricFunctions
            sage: SymmetricFunctions(F).equivariant(F.igen, F.shift_auto_on_element)
            Equivariant Symmetric Functions over Polynomial Ring over Rational Field with variables indexed by integers
        """
        return "Equivariant Symmetric Functions over {}".format(self.base_ring())

    def symmetric_function_ring(self):
        r"""
        The symmetric function ring underlying ``self``.

        EXAMPLES::

            sage: from sage.rings.polynomial.zpoly import FunctionFieldIntegerGenerators
            sage: F = FunctionFieldIntegerGenerators(QQ,('a','b'))
            sage: from sage.combinat.sf.sf import SymmetricFunctions
            sage: SymmetricFunctions(F).equivariant(F.igen, F.shift_auto_on_element).symmetric_function_ring()
            Symmetric Functions over Function Field over Rational Field with generators indexed by the integers
        """
        return self._sym

    def base_ring(self):
        r"""
        The base ring of the symmetric function ring underlying ``self``.

        EXAMPLES::

            sage: from sage.rings.polynomial.zpoly import FunctionFieldIntegerGenerators
            sage: F = FunctionFieldIntegerGenerators(QQ,('a','b'))
            sage: from sage.combinat.sf.sf import SymmetricFunctions
            sage: SymmetricFunctions(F).equivariant(F.igen, F.shift_auto_on_element).base_ring()
            Function Field over Rational Field with generators indexed by the integers
        """
        return self._sym.base_ring()

    def equivariant_parameter(self, i):
        r"""
        The `i`-th equivariant parameter.

        INPUT:

        - `i` -- An integer.

        EXAMPLES::

            sage: from sage.rings.polynomial.zpoly import FunctionFieldIntegerGenerators
            sage: F = FunctionFieldIntegerGenerators(QQ,('a','b'))
            sage: from sage.combinat.sf.sf import SymmetricFunctions
            sage: Eq = SymmetricFunctions(F).equivariant(F.igen, F.shift_auto_on_element)
            sage: Eq.equivariant_parameter(-2)
            b_2
            sage: Eq.equivariant_parameter(4)
            a_4
        """
        return self._equivariant_parameter(i)

    def shift_auto(self, r, z):
        r"""
        The image of `z` under shifting the base ring by `r`.

        INPUT:

        - `r` -- an integer
        - `z` -- an element of the base ring.

        EXAMPLES::

            sage: from sage.rings.polynomial.zpoly import FunctionFieldIntegerGenerators
            sage: F = FunctionFieldIntegerGenerators(QQ,('a','b'))
            sage: from sage.combinat.sf.sf import SymmetricFunctions
            sage: Eq = SymmetricFunctions(F).equivariant(F.igen, F.shift_auto_on_element)
            sage: Eq.shift_auto(4, Eq.equivariant_parameter(-2))
            a_2
        """
        return self._shift(r, z)

    def do_plethysm(self, f, func, codomain):
        r"""
        Replace power sums in `f` by values of ``func``.

        INPUT:

        - `f` -- A symmetric function in the p basis
        - ``func`` -- Function from first several positive integers to ``codomain``.
        - ``codomain`` -- The codomain of the plethysm

        EXAMPLES::

            sage: from sage.rings.polynomial.zpoly import PolynomialRingIntegerGenerators
            sage: F = PolynomialRingIntegerGenerators(QQ,('a','b'))
            sage: from sage.combinat.sf.sf import SymmetricFunctions
            sage: Sym = SymmetricFunctions(F)
            sage: Eq = Sym.equivariant(F.igen, F.shift_auto_on_element)
            sage: func = lambda i: Eq._p[i] + F.igen(1)**i
            sage: Eq._s(Eq.do_plethysm(Eq.equivariant_schur_function([2]), func , F))
            (a_1^2-a_1*b_0)*s[] + (a_1-b_0)*s[1] + s[2]
        """
        f = self._p(f)
        max = 0
        for part in f.support():
            if part.length() > 0:
                if part[0] > max:
                    max = part[0]
        return codomain.sum([f[part]*codomain.prod([func(i) for i in part]) for part in f.support()])

    def localize_symmetric_function(self, f, w):
        r"""
        Localize the symmetric function `f` at the partition ``part``.

        INPUT:

        - `f` -- a symmetric function
        - `w` -- an instance of :class:`PermutationOfIntegers`

        If `a_i` is sent to `a_{1-i}` for all `i` then
        the value of the equivariant Schur indexed by `\lambda` at the partition `\mu`,
        is the localization such that the diagonal terms are the product of left inversions.

        EXAMPLES::

            sage: from sage.rings.polynomial.zpoly import PolynomialRingIntegerGenerators
            sage: F = PolynomialRingIntegerGenerators(QQ,('a','b'))
            sage: from sage.combinat.sf.sf import SymmetricFunctions
            sage: Sym = SymmetricFunctions(F)
            sage: Eq = Sym.equivariant(F.igen, F.shift_auto_on_element)
            sage: from sage.combinat.permutation_of_integers import PermutationOfIntegers                          
            sage: w = PermutationOfIntegers(-1,3,[1,3,-1,0,2])
            sage: Eq.localize_symmetric_function(Eq.equivariant_schur_function([3,2]),w).factor()
            (-1) * (-b_2 + b_1) * (-a_1 + b_0) * (a_1 - b_2) * (-a_2 + b_0) * (a_2 - b_2)
            sage: w = PermutationOfIntegers(0,2,[2,0,1])
            sage: Eq.equivariant_schur_function([2])
            (-b_0)*s[1] + s[2]
            sage: Eq.localize_symmetric_function(Eq.equivariant_schur_function([2]),w).factor()
            (-b_1 + b_0) * (a_1 - b_1)

            sage: from sage.rings.polynomial.zpoly import PolynomialRingIntegerGenerators
            sage: F = PolynomialRingIntegerGenerators(QQ,('a','b'))
            sage: from sage.combinat.sf.sf import SymmetricFunctions
            sage: Sym = SymmetricFunctions(F)
            sage: Eq = Sym.equivariant(F.igen, F.shift_auto_on_element, False)
            sage: Eq.localize_symmetric_function(Eq.equivariant_schur_function([1]),PermutationOfIntegers(0,1,[1,0]))
            a_1 - b_0
            sage: Eq.localize_symmetric_function(Eq.equivariant_schur_function([2]),PermutationOfIntegers(0,2,[2,0,1])).factor()
            (-a_2 + b_0) * (-a_2 + a_1)
        """
        par = f.parent()
        F = self.base_ring()
        f = self._p(f)
        if self._positive:
            pos = [self.equivariant_parameter(1-w.value(i)) for i in w.nonpos_to_pos_set()]
            neg = [self.equivariant_parameter(1-w.value(i)) for i in w.pos_to_nonpos_set()]
        else:
            pos = [self.equivariant_parameter(w.value(i)) for i in w.nonpos_to_pos_set()]
            neg = [self.equivariant_parameter(w.value(i)) for i in w.pos_to_nonpos_set()]
        func = lambda i: F.sum([v**i for v in pos])-F.sum([v**i for v in neg])
        return self.do_plethysm(f, func, F)

    def shift_symmetric_function(self, r, f):
        r"""
        Apply the `r`-th shift to a symmetric function `f`.

        INPUT:

        - `r` -- an integer
        - `f` -- a symmetric function

        EXAMPLES::

            sage: from sage.rings.polynomial.zpoly import PolynomialRingIntegerGenerators
            sage: F = PolynomialRingIntegerGenerators(QQ,('a','b'))
            sage: from sage.combinat.sf.sf import SymmetricFunctions
            sage: Sym = SymmetricFunctions(F)
            sage: Eq = Sym.equivariant(F.igen, F.shift_auto_on_element)
            sage: s = Sym.s()
            sage: Eq.shift_symmetric_function(-3, s[2])
            (b_2*b_1+b_2*b_0+b_1*b_0)*s[] + (-b_2-b_1-b_0)*s[1] + s[2]
            sage: Eq.shift_symmetric_function(2, Eq.equivariant_parameter(4)**2*s[2])
            (a_6^2*a_2^2+a_6^2*a_2*a_1+a_6^2*a_1^2)*s[] + (a_6^2*a_2+a_6^2*a_1)*s[1] + a_6^2*s[2]
            sage: Eq = Sym.equivariant(F.igen, F.shift_auto_on_element,positive=False)
            sage: s = Sym.s()
            sage: Eq.shift_symmetric_function(-3, s[2])
            (b_2^2+b_2*b_1+b_1^2+b_2*b_0+b_1*b_0+b_0^2)*s[] + (b_2+b_1+b_0)*s[1] + s[2]
            sage: Eq.shift_symmetric_function(2, Eq.equivariant_parameter(4)**2*s[2])
            a_6^2*a_2*a_1*s[] + (-a_6^2*a_2-a_6^2*a_1)*s[1] + a_6^2*s[2]
        """
        if r == 0:
            return f
        par = f.parent()
        # apply the shift automorphism to the coefficients
        # change to power sums and do the usual kind of thing for plethysm
        # the built-in plethysm doesn't work for some reason
        f = f.map_coefficients(lambda y: self.shift_auto(r, y))
        if r > 0:
            sign = 1
            include_vars = [self.equivariant_parameter(i) for i in range(1,r+1)]
        else:
            sign = -1
            include_vars = [self.equivariant_parameter(-i) for i in range(0,-r)]
        if not self._positive:
            sign = - sign
        func = lambda i: self._p[i] + sign * self._p.one() * sum([v**i for v in include_vars])
        return self.do_plethysm(f, func, par)

    def equivariant_schur_function(self, part):
        r"""
        The equivariant Schur function indexed by the partition ``part``, expressed in the Schur basis.

        EXAMPLES::

            sage: from sage.rings.polynomial.zpoly import PolynomialRingIntegerGenerators
            sage: F = PolynomialRingIntegerGenerators(QQ,('a','b'))
            sage: from sage.combinat.sf.sf import SymmetricFunctions
            sage: Eq = SymmetricFunctions(F).equivariant(F.igen, F.shift_auto_on_element)
            sage: Eq.equivariant_schur_function([2])
            (-b_0)*s[1] + s[2]
            sage: Eq.equivariant_schur_function([2,1])
            (-a_1*b_0)*s[1] + (-b_0)*s[1, 1] + a_1*s[2] + s[2, 1]
            sage: Eq = SymmetricFunctions(F).equivariant(F.igen, F.shift_auto_on_element,positive=False)
            sage: Eq.equivariant_schur_function([2])
            (-a_1)*s[1] + s[2]
            sage: Eq.equivariant_schur_function([2,1])
            (-a_1*b_0)*s[1] + (-a_1)*s[1, 1] + b_0*s[2] + s[2, 1]
        """
        return self._esf(Partition(part))

    @cached_method
    def _esf(self, part):
        if part.length() == 0:
            return self._s.one()
        k = part.length()
        def entry(i,j,p):
            r"""
            Entry of determinant for factorial Schur function.

            Indices are zero-based.
            """
            v = p-i+j
            if v < 0:
                return self._s.zero()
            if v == 0:
                return self._s.one()
            sh = 1 + i - p
            if not self._positive:
                sh = - sh
            return self.shift_symmetric_function(sh, self._s[v])
        m = [[entry(i,j,part[i]) for j in range(k)] for i in range(k)]
        return MatrixSpace(self._s, k, k)(m).det()

    @cached_method
    def _s_to_es(self, part):
        r"""
        The family expanding a Schur function into equivariant Schur functions.

        INPUT:

        - ``part`` -- a partition

        OUTPUT:

        A family whose keys are the partitions contained in ``part`` and whose
        values are the coefficients of the expansion of the Schur function indexed by ``part``
        into equivariant Schur functions.

        EXAMPLES::

            sage: from sage.rings.polynomial.zpoly import PolynomialRingIntegerGenerators
            sage: F = PolynomialRingIntegerGenerators(QQ,('a','b'))
            sage: from sage.combinat.sf.sf import SymmetricFunctions
            sage: Sym = SymmetricFunctions(F)
            sage: Eq = Sym.equivariant(F.igen, F.shift_auto_on_element)
            sage: Eq._s_to_es(Partition([2,1]))
            Finite family {[2]: -a_1, [1, 1]: b_0, [2, 1]: 1, [1]: -a_1*b_0}
            sage: Eq._s_to_es(Partition([2]))
            Finite family {[2]: 1, [1]: b_0}
        """
        rest_esf = self._s[part] - self.equivariant_schur_function(part)
        exp = self._s[part] # faking the expansion using the Schur basis as placeholders
        for pt in rest_esf.support():
            cf = rest_esf[pt]
            fam = self._s_to_es(pt)
            for ptn in fam.keys():
                exp = exp + cf * fam[ptn] * self._s[ptn]
        return Family(exp.support(), lambda pt: exp[pt])

    def equivariant_q_function(self, part):
        r"""
        An equivariant analogue of the power sum basis, expressed in equivariant Schurs.

        EXAMPLES::

            sage: from sage.rings.polynomial.zpoly import PolynomialRingIntegerGenerators
            sage: F = PolynomialRingIntegerGenerators(QQ,('a','b'))
            sage: from sage.combinat.sf.sf import SymmetricFunctions
            sage: Eq = SymmetricFunctions(F).equivariant(F.igen, F.shift_auto_on_element)
            sage: Eq.equivariant_q_function([2,1])
            (a_2-b_0)*EQs[1, 1] - EQs[1, 1, 1] + (-a_1+b_1)*EQs[2] + EQs[3]
            sage: Eq = SymmetricFunctions(F).equivariant(F.igen, F.shift_auto_on_element,positive=False)
            sage: Eq.equivariant_q_function([2,1])
            (-a_1+b_1)*EQs[1, 1] - EQs[1, 1, 1] + (a_2-b_0)*EQs[2] + EQs[3]
        """
        return self._eqq(Partition(part))

    @cached_method
    def _eqq_single_part(self, r):
        r"""
        The alternating sum of equivariant Schurs for hooks of size `r`.

        EXAMPLES::

            sage: from sage.rings.polynomial.zpoly import PolynomialRingIntegerGenerators
            sage: F = PolynomialRingIntegerGenerators(QQ,('a','b'))
            sage: from sage.combinat.sf.sf import SymmetricFunctions
            sage: Eq = SymmetricFunctions(F).equivariant(F.igen, F.shift_auto_on_element)
            sage: [Eq._eqq_single_part(i) for i in range(4)]
            [EQs[], EQs[1], -EQs[1, 1] + EQs[2], EQs[1, 1, 1] - EQs[2, 1] + EQs[3]]
        """
        es = self.es()
        if r == 0:
            return es.one()
        ans = es[r]
        sign = 1
        for i in reversed(range(1,r)):
            sign = -sign
            ans = ans + sign * es(Partition([i]+[1]*(r-i)))
        return ans

    @cached_method
    def _eqq(self, part):
        r"""
        Basis of product of the above `q_i`, expanded into equivariant Schurs.

            sage: from sage.rings.polynomial.zpoly import PolynomialRingIntegerGenerators
            sage: F = PolynomialRingIntegerGenerators(QQ,('a','b'))
            sage: from sage.combinat.sf.sf import SymmetricFunctions
            sage: Eq = SymmetricFunctions(F).equivariant(F.igen, F.shift_auto_on_element)
            sage: [[part, Eq._eqq(part)] for part in Partitions(3)]
            [[[3], EQs[1, 1, 1] - EQs[2, 1] + EQs[3]], [[2, 1], (a_2-b_0)*EQs[1, 1] - EQs[1, 1, 1] + (-a_1+b_1)*EQs[2] + EQs[3]], [[1, 1, 1], (a_1^2-2*a_1*b_0+b_0^2)*EQs[1] + (-a_2-a_1+2*b_0)*EQs[1, 1] + EQs[1, 1, 1] + (-2*a_1+b_1+b_0)*EQs[2] + 2*EQs[2, 1] + EQs[3]]]
        """
        es = self.es()
        if part.length() == 0:
            return es.one()
        return es.prod([self._eqq_single_part(i) for i in part])

    def es(self):
        r"""
        Equivariant (factorial) Schur basis.

        EXAMPLES::

            sage: from sage.rings.polynomial.zpoly import PolynomialRingIntegerGenerators
            sage: F = PolynomialRingIntegerGenerators(QQ,('a','b'))
            sage: from sage.combinat.sf.sf import SymmetricFunctions
            sage: SymmetricFunctions(F).equivariant(F.igen, F.shift_auto_on_element).es()
            Symmetric Functions over Polynomial Ring over Rational Field with variables indexed by integers in the Equivariant s basis
        """
        return Equivariant_s(self)

    def eq(self):
        r"""
        Equivariant q basis.

        See :meth:`EquivariantSymmetricFunctions.equivariant_q_function`.

        EXAMPLES::

            sage: from sage.rings.polynomial.zpoly import PolynomialRingIntegerGenerators
            sage: F = PolynomialRingIntegerGenerators(QQ,('a','b'))
            sage: from sage.combinat.sf.sf import SymmetricFunctions
            sage: SymmetricFunctions(F).equivariant(F.igen, F.shift_auto_on_element).eq()
            Symmetric Functions over Polynomial Ring over Rational Field with variables indexed by integers in the Equivariant q basis
        """
        return Equivariant_q(self)

    def rees_ring(self):
        r"""
        Associated rees ring.

        EXAMPLES::

            sage: from sage.rings.polynomial.zpoly import PolynomialRingIntegerGenerators
            sage: F = PolynomialRingIntegerGenerators(QQ,('a','b'))
            sage: from sage.combinat.sf.sf import SymmetricFunctions
            sage: Eq = SymmetricFunctions(F).equivariant(F.igen, F.shift_auto_on_element)
            sage: Eq.rees_ring()
            Symmetric Functions over Fraction Field of Univariate Polynomial Ring in d over Rational Field
        """
        from sage.combinat.sf.sf import SymmetricFunctions
        return SymmetricFunctions(PolynomialRing(self.base_ring().base_ring(), ['d']).fraction_field())

    def rees_gen(self):
        r"""
        Associated rees ring generator.

        EXAMPLES::

            sage: from sage.rings.polynomial.zpoly import PolynomialRingIntegerGenerators
            sage: F = PolynomialRingIntegerGenerators(QQ,('a','b'))
            sage: from sage.combinat.sf.sf import SymmetricFunctions
            sage: SymmetricFunctions(F).equivariant(F.igen, F.shift_auto_on_element).rees_gen()
            d
        """
        return self.rees_ring().base_ring().gen()

    @cached_method
    def _to_rees_on_basis(self, part):
        r"""
        The image of the equivariant q function indexed by the partition ``part``
        in the rees ring.

        EXAMPLES::

            sage: from sage.rings.polynomial.zpoly import PolynomialRingIntegerGenerators
            sage: F = PolynomialRingIntegerGenerators(QQ,('a','b'))
            sage: from sage.combinat.sf.sf import SymmetricFunctions
            sage: Eq = SymmetricFunctions(F).equivariant(F.igen, F.shift_auto_on_element)
            sage: Eq._to_rees_on_basis(Partition([2]))
            -1/d*s[1, 1, 1]
            sage: Eq._to_rees_on_basis(Partition([2,1]))
            -1/d^2*s[1, 1, 1, 1, 1] - 1/d^2*s[2, 1, 1, 1] - 1/d^2*s[2, 2, 1]
        """
        rees = self.rees_ring()
        d = self.rees_gen()
        s = rees.s()
        if self._positive:
            e = rees.e()
            return s(e.prod([(-1)**(i+1)*(1/d)*e[i+1] for i in part]))
        else:
            return s.prod([(1/d)*s[i+1] for i in part])

    def to_rees(self, f):
        r"""
        Image of equivariant symmetric function `f` in the rees ring.

        EXAMPLES::

            sage: from sage.rings.polynomial.zpoly import PolynomialRingIntegerGenerators
            sage: F = PolynomialRingIntegerGenerators(QQ,('a','b'))
            sage: from sage.combinat.sf.sf import SymmetricFunctions
            sage: Eq = SymmetricFunctions(F).equivariant(F.igen, F.shift_auto_on_element)
            sage: Eq.to_rees(Eq.eq()[2])
            -1/3/d*p[3]
            sage: Eq.to_rees(Eq.es()[2])
            -1/4*p[2] + 1/8/d^2*p[2, 2] - 1/6/d*p[3]
        """
        eq = self.eq()
        f = eq(f)
        rees = self.rees_ring()
        d = self.rees_gen()
        s = rees.s()
        Q = self.base_ring()
        ans = s.zero()
        for part in f.support():
            cf = Q.forget_to_loop_rotation(d,Q(f[part]))
            ans = ans + cf * self._to_rees_on_basis(part)
        p = rees.p()
        ans = p(ans)
        return p.sum([ans[part]*p[part] for part in ans.support() if part.length()==0 or part[part.length()-1]>1])

    @cached_method
    def _single_power_to_rees(self, k):
        r"""
        Image of power sum in Rees ring.

        MATH::

            \tilde{p}_k = \sum_{j=1}^{k-1} \binomial{k}{j} (-\delta)^{k-j} p_j

        where `\tilde{p}` is the power sum polynomial in the Rees ring and p is the usual power sum.

        EXAMPLES::

            sage: from sage.rings.polynomial.zpoly import PolynomialRingIntegerGenerators
            sage: F = PolynomialRingIntegerGenerators(QQ,('a','b'))
            sage: from sage.combinat.sf.sf import SymmetricFunctions
            sage: Eq = SymmetricFunctions(F).equivariant(F.igen, F.shift_auto_on_element, positive=False)
            sage: Eq._single_power_to_rees(1)
            -1/(2*d)*p[2]
            sage: Eq._single_power_to_rees(2)
            -1/2*p[2] - 1/(3*d)*p[3]
            sage: Eq._single_power_to_rees(3)
            -1/4*d*p[2] - 1/2*p[3] - 1/(4*d)*p[4]
            sage: Eq = SymmetricFunctions(F).equivariant(F.igen, F.shift_auto_on_element)
            sage: Eq._single_power_to_rees(1)
            1/(2*d)*p[2]
            sage: Eq._single_power_to_rees(2)
            -1/2*p[2] + 1/(3*d)*p[3]
        """
        from sage.arith.misc import binomial
        rees = self.rees_ring()
        d = self.rees_gen()
        p = rees.p()
        ans = p[k+1]
        for j in range(1,k):
            ans = ans - (-d)**(k+1-j) * binomial(k+1,j) * self._single_power_to_rees(j)
        ans = ans/((k+1)*d)
        if self._positive:
            return ans
        return -ans

    @cached_method
    def _to_rees_on_basis_powers(self, ptn):
        r"""
        The image of the power sum indexed by the partition ``ptn``
        in the rees ring.

        EXAMPLES::

            sage: from sage.rings.polynomial.zpoly import PolynomialRingIntegerGenerators
            sage: F = PolynomialRingIntegerGenerators(QQ,('a','b'))
            sage: from sage.combinat.sf.sf import SymmetricFunctions
            sage: Eq = SymmetricFunctions(F).equivariant(F.igen, F.shift_auto_on_element, positive=False)
            sage: Eq._to_rees_on_basis_powers(Partition([2]))
            -1/2*p[2] - 1/(3*d)*p[3]
            sage: Eq._to_rees_on_basis_powers(Partition([2,1]))
            1/2/(2*d)*p[2, 2] + 1/(6*d^2)*p[3, 2]
        """
        rees = self.rees_ring()
        d = self.rees_gen()
        p = rees.p()
        return p.prod([self._single_power_to_rees(i) for i in ptn])

    def to_rees_by_powers(self, f):
        r"""
        Image of equivariant symmetric function `f` in the rees ring.

        EXAMPLES::

            sage: from sage.rings.polynomial.zpoly import PolynomialRingIntegerGenerators
            sage: F = PolynomialRingIntegerGenerators(QQ,('a','b'))
            sage: from sage.combinat.sf.sf import SymmetricFunctions
            sage: Eq = SymmetricFunctions(F).equivariant(F.igen, F.shift_auto_on_element, positive=False)
            sage: rees = Eq.rees_ring()
            sage: p = rees.p()
            sage: s = rees.s()
            sage: d = Eq.rees_gen()
            sage: d*Eq.to_rees_by_powers(Eq.eq()[2])
            -1/3*p[3]
            sage: p(s[1,1,1])
            1/6*p[1, 1, 1] - 1/2*p[2, 1] + 1/3*p[3]
            sage: d*Eq.to_rees_by_powers(Eq.eq()[3])
            1/8*p[2, 2] - 1/4*p[4]
            sage: p(s[1,1,1,1])
            1/24*p[1, 1, 1, 1] - 1/4*p[2, 1, 1] + 1/8*p[2, 2] + 1/3*p[3, 1] - 1/4*p[4]            
        """
        rees = self.rees_ring()
        d = self.rees_gen()
        p = rees.p()
        Q = self.base_ring()
        ans = p.zero()
        peq = self._p
        f = peq(f)
        for part in f.support():
            cf = Q.forget_to_loop_rotation(d,Q(f[part]),positive=self._positive)
            ans = ans + cf * self._to_rees_on_basis_powers(part)
        return ans

class EquivariantSymmetricFunctions_generic(sfa.SymmetricFunctionAlgebra_generic):
    r"""
    A class for a general basis of the equivariant symmetric functions.

    INPUT:

    - ``self`` -- an equivariant basis of symmetric functions
    - ``equivariant`` -- a family of equivariant symmetric function bases

    EXAMPLES::

        sage: from sage.rings.polynomial.zpoly import PolynomialRingIntegerGenerators
        sage: F = PolynomialRingIntegerGenerators(QQ,('a','b'))
        sage: from sage.combinat.sf.sf import SymmetricFunctions
        sage: SymmetricFunctions(F).equivariant(F.igen, F.shift_auto_on_element).es()
        Symmetric Functions over Polynomial Ring over Rational Field with variables indexed by integers in the Equivariant s basis
    """
    def __init__(self, equivariant):
        name = self.__class__.__name__[len("Equivariant_"):] # note the naming convention being used
        sfa.SymmetricFunctionAlgebra_generic.__init__(
            self, equivariant._sym,
            basis_name = "Equivariant " + name,
            prefix = "EQ"+name)
        self._equivariant = equivariant
        self._s = self._equivariant._s
        self._p = self._equivariant._p

    def equivariant_family(self):
        return self._equivariant

    class Element(sfa.SymmetricFunctionAlgebra_generic.Element):
        def separate(self):
            r"""
            Separate an element into nonequivariant symmetric functions and equivariant terms

            EXAMPLES::

                sage: from sage.rings.polynomial.zpoly import PolynomialRingIntegerGenerators
                sage: F = PolynomialRingIntegerGenerators(QQ,('a','b'))
                sage: from sage.combinat.sf.sf import SymmetricFunctions
                sage: Eq = SymmetricFunctions(F).equivariant(F.igen, F.shift_auto_on_element)
                sage: es = Eq.es()
                sage: es[2].separate()
                b_0*s[] # s[1] + s[] # s[1, 1] + (-b_0)*s[1] # s[] - s[1] # s[1] + s[2] # s[]
                sage: eq = Eq.eq()
                sage: eq[1].separate()
                -s[] # s[1] + s[1] # s[]
                sage: eq[2].separate()
                (a_1+b_0)*s[] # s[1] + s[] # s[1, 1] - s[] # s[2] + (-a_1-b_0)*s[1] # s[] - s[1, 1] # s[] + s[2] # s[]
                sage: eq[3].separate()
                (-a_2*a_1-a_1*b_0-b_1*b_0)*s[] # s[1] + (-a_1-b_1-b_0)*s[] # s[1, 1] - s[] # s[1, 1, 1] + (a_2+a_1+b_0)*s[] # s[2] + s[] # s[2, 1] - s[] # s[3] + (a_2*a_1+a_1*b_0+b_1*b_0)*s[1] # s[] + (-a_2+b_1)*s[1] # s[1] + (a_2+a_1+b_0)*s[1, 1] # s[] + s[1, 1, 1] # s[] + (-a_1-b_1-b_0)*s[2] # s[] - s[2, 1] # s[] + s[3] # s[]
            """
            s = self.parent().equivariant_family().symmetric_function_ring().s()
            return s(self).coproduct().apply_multilinear_morphism(lambda a,b: tensor([a,b.antipode()]))

class Equivariant_s(EquivariantSymmetricFunctions_generic):
    r"""
    The equivariant Schur basis.

    INPUT:

    - ``self`` -- The basis of equivariant Schur functions
    - ``equivariant`` -- A family of bases of equivariant symmetric functions

    EXAMPLES::

        sage: from sage.rings.polynomial.zpoly import PolynomialRingIntegerGenerators
        sage: F = PolynomialRingIntegerGenerators(QQ,('a','b'))
        sage: from sage.combinat.sf.sf import SymmetricFunctions
        sage: Eq = SymmetricFunctions(F).equivariant(F.igen, F.shift_auto_on_element)
        sage: es = Eq.es()
        sage: es[2]
        EQs[2]
        sage: s = Eq.symmetric_function_ring().s()
        sage: s(es[2])
        (-b_0)*s[1] + s[2]
    """
    def __init__(self, equivariant):
        r"""
        Initialize ``self``.
        """
        EquivariantSymmetricFunctions_generic.__init__(self, equivariant)
        # set up coercion between self and the s basis
        category = ModulesWithBasis(self.base_ring())
        self._self_to_s = self.module_morphism(on_basis = lambda part: equivariant.equivariant_schur_function(part), codomain=self._s)
        self._s_to_self = self._s.module_morphism(on_basis = self._s_to_self_on_basis, codomain=self)
        self._s.register_coercion(self._self_to_s)
        self   .register_coercion(self._s_to_self)

    def _s_to_self_on_basis(self, part):
        r"""
        Given a partition, return the expansion of the corresponding Schur function in the equivariant
        Schur basis.

        EXAMPLES::

            sage: from sage.rings.polynomial.zpoly import PolynomialRingIntegerGenerators
            sage: F = PolynomialRingIntegerGenerators(QQ,('a','b'))
            sage: from sage.combinat.sf.sf import SymmetricFunctions
            sage: Eq = SymmetricFunctions(F).equivariant(F.igen, F.shift_auto_on_element)
            sage: es = Eq.es()
            sage: es._s_to_self_on_basis(Partition([2,1]))
            (-a_1*b_0)*EQs[1] + b_0*EQs[1, 1] + (-a_1)*EQs[2] + EQs[2, 1]
        """
        fam = self._equivariant._s_to_es(part)
        return self.sum_of_terms([(pt,fam[pt]) for pt in fam.keys()], distinct=True)

    class Element(EquivariantSymmetricFunctions_generic.Element):
        pass

class Equivariant_q(EquivariantSymmetricFunctions_generic):
    r"""
    This is the basis of products of elements `q_i` defined above.

    See :meth:`EquivariantSymmetricFunctions.equivariant_q_function`.

    INPUT:

    - ``self`` -- The basis of equivariant q functions
    - ``equivariant`` -- A family of bases of equivariant symmetric functions

    EXAMPLES::

        sage: from sage.rings.polynomial.zpoly import PolynomialRingIntegerGenerators
        sage: F = PolynomialRingIntegerGenerators(QQ,('a','b'))
        sage: from sage.combinat.sf.sf import SymmetricFunctions
        sage: Eq = SymmetricFunctions(F).equivariant(F.igen, F.shift_auto_on_element)
        sage: eq = Eq.eq()
        sage: eq[2]
        EQq[2]
        sage: es = Eq.es()
        sage: es(eq[2])
        -EQs[1, 1] + EQs[2]
    """
    def __init__(self, equivariant):
        r"""
        Initialize ``self``.
        """
        EquivariantSymmetricFunctions_generic.__init__(self, equivariant)
        # set up coercion between self and the equivariant Schur basis
        self._es = equivariant.es()
        category = ModulesWithBasis(self.base_ring())
        self._self_to_es = self.module_morphism(on_basis = lambda part: equivariant.equivariant_q_function(part), codomain=self._es)
        self._es_to_self = self._es.module_morphism(on_basis = self._es_to_eq, codomain=self)
        self._es.register_coercion(self._self_to_es)
        self    .register_coercion(self._es_to_self)

    @cached_method
    def _es_to_eq_degree(self, d):
        r"""
        The expansion of an equivariant Schur function into equivariant q-functions,
        for every partition of the integer `d`.

        EXAMPLES::

            sage: from sage.rings.polynomial.zpoly import PolynomialRingIntegerGenerators
            sage: F = PolynomialRingIntegerGenerators(QQ,('a','b'))
            sage: from sage.combinat.sf.sf import SymmetricFunctions
            sage: Sym = SymmetricFunctions(F)
            sage: Eq = Sym.equivariant(F.igen, F.shift_auto_on_element)
            sage: eq = Eq.eq()
            sage: eq._es_to_eq_degree(2)
            Finite family {[2]: (1/2*a_1-1/2*b_0)*EQq[1] + 1/2*EQq[1, 1] + 1/2*EQq[2], [1, 1]: (1/2*a_1-1/2*b_0)*EQq[1] + 1/2*EQq[1, 1] - 1/2*EQq[2]}
            sage: eq._es_to_eq_degree(3)
            Finite family {[1, 1, 1]: (1/3*a_2*a_1-1/6*a_1^2+1/6*a_1*b_1-1/3*a_2*b_0-1/6*a_1*b_0-1/6*b_1*b_0+1/3*b_0^2)*EQq[1] + (1/3*a_2+1/6*b_1-1/2*b_0)*EQq[1, 1] + 1/6*EQq[1, 1, 1] + (-1/3*a_2-1/6*a_1+1/6*b_1+1/3*b_0)*EQq[2] - 1/2*EQq[2, 1] + 1/3*EQq[3], [3]: (-1/6*a_2*a_1+1/3*a_1^2-1/3*a_1*b_1+1/6*a_2*b_0-1/6*a_1*b_0+1/3*b_1*b_0-1/6*b_0^2)*EQq[1] + (-1/6*a_2+1/2*a_1-1/3*b_1)*EQq[1, 1] + 1/6*EQq[1, 1, 1] + (1/6*a_2+1/3*a_1-1/3*b_1-1/6*b_0)*EQq[2] + 1/2*EQq[2, 1] + 1/3*EQq[3], [2, 1]: (1/6*a_2*a_1+1/6*a_1^2-1/6*a_1*b_1-1/6*a_2*b_0-1/3*a_1*b_0+1/6*b_1*b_0+1/6*b_0^2)*EQq[1] + (1/6*a_2+1/2*a_1-1/6*b_1-1/2*b_0)*EQq[1, 1] + 1/3*EQq[1, 1, 1] + (-1/6*a_2+1/6*a_1-1/6*b_1+1/6*b_0)*EQq[2] - 1/3*EQq[3]}
        """
        if d == 0:
            return Family(dict([[Partition([]),self.one()]]))
        ptns = [ptn for ptn in Partitions(d)]
        nptns = len(ptns)
        q_to_es = Family(ptns, lambda ptn: self._equivariant.equivariant_q_function(ptn))
        F = self._equivariant.base_ring()
        Mat = MatrixSpace(F, len(ptns), len(ptns))
        topdeg = Mat([[q_to_es[a][b] for b in ptns] for a in ptns])
        topdeginv = topdeg.inverse()
        topdeginv = Mat(topdeginv)
        dc = dict()
        for i in range(nptns):
            dc[ptns[i]] = self.sum_of_terms([(ptns[j], topdeginv[i][j]) for j in range(nptns)], distinct=True) + self(self._es[ptns[i]] - self._es.sum([topdeginv[i][j]*q_to_es[ptns[j]]  for j in range(nptns)]))
        return Family(dc)

    @cached_method
    def _es_to_eq(self, part):
        r"""
        An equivariant Schur function expressed in the equivariant q basis.

        INPUT:

        - ``part`` -- a partition

        EXAMPLES::

            sage: from sage.rings.polynomial.zpoly import PolynomialRingIntegerGenerators
            sage: F = PolynomialRingIntegerGenerators(QQ,('a','b'))
            sage: from sage.combinat.sf.sf import SymmetricFunctions
            sage: Sym = SymmetricFunctions(F)
            sage: Eq = Sym.equivariant(F.igen, F.shift_auto_on_element)
            sage: eq = Eq.eq()
            sage: eq._es_to_eq(Partition([2,1]))
            (1/6*a_2*a_1+1/6*a_1^2-1/6*a_1*b_1-1/6*a_2*b_0-1/3*a_1*b_0+1/6*b_1*b_0+1/6*b_0^2)*EQq[1] + (1/6*a_2+1/2*a_1-1/6*b_1-1/2*b_0)*EQq[1, 1] + 1/3*EQq[1, 1, 1] + (-1/6*a_2+1/6*a_1-1/6*b_1+1/6*b_0)*EQq[2] - 1/3*EQq[3]
            sage: Eq.es()(_)            
            EQs[2, 1]
        """
        d = part.size()
        return self._es_to_eq_degree(d)[part]
