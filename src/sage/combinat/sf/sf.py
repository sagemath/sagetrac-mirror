"""
Symmetric functions, with their multiple realizations
"""
from __future__ import absolute_import
# ****************************************************************************
#       Copyright (C) 2007 Mike Hansen <mhansen@gmail.com>
#                     2009-2012 Jason Bandlow <jbandlow@gmail.com>
#                     2012 Anne Schilling <anne at math.ucdavis.edu>
#                     2009-2012 Nicolas M. Thiery <nthiery at users.sf.net>
#                     2012 Mike Zabrocki <mike.zabrocki@gmail.com>
#                     2017-2019 Pauline Hubert <hubert.pauline@courrier.uqam.ca>
#                     2017-2019 Melodie Lapointe <lapointe.melodie@courrier.uqam.ca>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.categories.graded_hopf_algebras import GradedHopfAlgebras
from sage.categories.fields import Fields
from sage.categories.rings import Rings
from sage.combinat.partition import Partitions
from sage.combinat.free_module import CombinatorialFreeModule
from sage.rings.rational_field import QQ

from . import schur
from . import monomial
from . import powersum
from . import elementary
from . import homogeneous
from . import hall_littlewood
from . import jack
from . import macdonald
from . import llt

class SymmetricFunctions(UniqueRepresentation, Parent):
    r"""
    The abstract algebra of commutative symmetric functions

    .. rubric:: Symmetric Functions in Sage

    .. MODULEAUTHOR:: Jason Bandlow, Pauline Hubert, Melodie Lapointe, Anne Schilling, Nicolas M. Thiery, Mike Zabrocki

    The first part of this document consists of a tutorial on symmetric functions.
    It is then followed by the documentation.
    The aim of this tutorial is to present what can be done in Sage with symmetric functions.
    We suppose that the reader knows the basics about symmetric functions;
    more can be learned from the excellent resources
    [Mac1995]_ and [EnumComb2]_, Chapter 7.

    **Caveat:** In this tutorial, the term "symmetric functions" will
    mostly stand for "abstract" symmetric functions, in which variables
    are not made explicit. Indeed, for most practical calculations,
    variables need not appear. Moreover, one may show that this does not
    cause any trouble in the calculations.

    *Outputs printed in latex mode.*

    ::

        sage: %display latex           # not tested


    .. rubric:: For the impatient

    Before going into details with symmetric functions in Sage, here is
    a quick example of what we can do in Sage.

    We recall that the **complete homogeneous** symmetric function
    :math:`h_d` (for a given nonnegative integer :math:`d`)
    can be expressed in terms of the **power sum** symmetric
    functions :math:`p_{\mu}` (corresponding to the partitions
    :math:`\mu` of :math:`d`) by the formula:

    .. MATH:: h_d = \sum\limits_{\mu \vdash d} \dfrac{1}{z_{\mu}} p_{\mu} ,

    where :math:`z_\mu` is the number of "automorphisms" of a permutation
    having cycle structure :math:`\mu` (that is, the size of the
    centralizer of such a permutation).

    Here is how to obtain both sides of this equality in the ring
    :math:`\operatorname{Sym}` of symmetric functions over
    :math:`\QQ` ::

        sage: Sym = SymmetricFunctions(QQ)
        sage: Sym.inject_shorthands()
        Defining e as shorthand for Symmetric Functions over Rational Field in the elementary basis
        Defining f as shorthand for Symmetric Functions over Rational Field in the forgotten basis
        Defining h as shorthand for Symmetric Functions over Rational Field in the homogeneous basis
        Defining m as shorthand for Symmetric Functions over Rational Field in the monomial basis
        Defining p as shorthand for Symmetric Functions over Rational Field in the powersum basis
        Defining s as shorthand for Symmetric Functions over Rational Field in the Schur basis

    ::

        sage: p(h[6])
        1/720*p[1, 1, 1, 1, 1, 1] + 1/48*p[2, 1, 1, 1, 1] + 1/16*p[2, 2, 1, 1] + 1/48*p[2, 2, 2] + 1/18*p[3, 1, 1, 1] + 1/6*p[3, 2, 1] + 1/18*p[3, 3] + 1/8*p[4, 1, 1] + 1/8*p[4, 2] + 1/5*p[5, 1] + 1/6*p[6]
        sage: sum((1/Partition(i).aut())*p(i) for i in Partitions(6).list())
        1/720*p[1, 1, 1, 1, 1, 1] + 1/48*p[2, 1, 1, 1, 1] + 1/16*p[2, 2, 1, 1] + 1/48*p[2, 2, 2] + 1/18*p[3, 1, 1, 1] + 1/6*p[3, 2, 1] + 1/18*p[3, 3] + 1/8*p[4, 1, 1] + 1/8*p[4, 2] + 1/5*p[5, 1] + 1/6*p[6]

    On the first line, we defined ``Sym`` to be the ring
    :math:`\operatorname{Sym}`. On the second, we used
    ``Sym.inject_shorthands()`` to introduce the standard
    notations for its classical bases -- ``e`` for the
    elementary symmetric functions, ``p`` for the power-sum
    symmetric functions, etc. (We could have achieved the
    same by defining each basis manually, using
    ``e = Sym.e()``, ``f = Sym.f()``, etc.).
    Then, we expanded :math:`h_6`
    in the power-sum basis (by writing ``p(h[6])``; here,
    ``h[6]`` computes `h_6`, and then the ``p`` converts
    it into the  power-sum basis). Finally, we computed
    the right hand side directly in the power-sum basis.

    .. rubric:: Abstract symmetric functions

    There are two equivalent ways to define the algebra
    :math:`\operatorname{Sym}` of symmetric functions over a
    given (commutative) base ring :math:`R`.
    The first way defines it as a certain subalgebra of the ring
    :math:`R[[x_1, x_2, x_3, \ldots]]` of formal power series
    in countably many variables :math:`x_1, x_2, x_3, \ldots`
    over :math:`R`; namely, it is the subalgebra consisting of
    all power series that are symmetric and of bounded degree.
    The second way is to define it as the unique free commutative
    graded connected algebra over :math:``R` with one generator
    in each positive degree (i.e., as a polynomial ring in
    countably many variables :math:`e_1, e_2, e_3, \ldots` over
    :math:`R`, where each variable :math:`e_d` is assigned
    degree :math:`d`).
    The equivalence of these two definitions follows from the
    "fundamental theorem of symmetric functions", which says
    that if we define :math:`\operatorname{Sym}` in the first
    way (i.e., as a ring of symmetric power series of bounded
    degree), then the *elementary symmetric functions*
    :math:`e_d = \sum_{i_1 < i_2 < \cdots < i_d} x_{i_1} x_{i_2} \cdots x_{i_d}`
    form an algebraically independent generating set of the
    :math:`R`-algebra :math:`\operatorname{Sym}` (whence
    :math:`\operatorname{Sym}` can be identified with the
    polynomial ring in :math:`e_1, e_2, e_3, \ldots`).

    Yet another equivalent definition of the ring
    :math:`\operatorname{Sym}` constructs it as the inverse
    limit (in the category of graded algebras) of the algebra
    of symmetric polynomials in :math:`n` variables as
    :math:`n \rightarrow \infty`. This is due to the fact that
    a symmetric power series :math:`f(x_1, x_2, x_3, \ldots)`
    of bounded degree can always be evaluated at finitely many
    variables :math:`x_1, x_2, \ldots, x_n` (by setting all
    "higher" variables :math:`x_{n+1}, x_{n+2}, x_{n+3}, \ldots`
    to zero), yielding a symmetric polynomial
    :math:`f(x_1, x_2, \ldots, x_n, 0, 0, 0, \ldots)` in
    :math:`x_1, x_2, \ldots, x_n`.

    We first describe how to manipulate symmetric functions
    without using the variables :math:`x_1, x_2, x_3, \ldots`.
    Such functions are linear combinations of one of the six classical
    bases of symmetric functions; all indexed by integer partitions
    :math:`\mu = \left(\mu_1, \mu_2, \ldots \mu_k\right)`
    (with :math:`\mu_1 \geq \mu_2 \geq \cdots \geq \mu_k > 0`).

    -   The **power sum** symmetric functions :math:`p_\mu=p_{\mu_1}p_{\mu_2}\cdots p_{\mu_k}`
    -   The **(complete) homogeneous** symmetric functions :math:`h_\mu=h_{\mu_1}h_{\mu_2}\cdots h_{\mu_k}`
    -   The **elementary** symmetric functions :math:`e_\mu=e_{\mu_1}e_{\mu_2}\cdots e_{\mu_k}`
    -   The **monomial** symmetric functions :math:`m_{\mu}`
    -   The **Schur** functions :math:`s{\mu}`
    -   The **forgotten** symmetric functions :math:`f_{\mu}`

    (or various other bases).

    As our base ring, we take the field :math:`\QQ (q,t)`
    of rational functions in two variables :math:`q` and
    :math:`t` with rational coefficients. ::

        sage: F = QQ['q','t'].fraction_field()
        sage: F.inject_variables()
        Defining q, t
        sage: Symqt = SymmetricFunctions(F)
        sage: Symqt.inject_shorthands(verbose=False)

    In the first line here, we have defined ``F`` to be the
    field :math:`\QQ (q, t)`.
    In the second, we have "injected" the variables :math:`q` and 
    :math:`t` in order to make them available as ``q`` and ``t``.
    Then, we have defined ``Symqt`` to be the ring of symmetric
    functions over :math:`\QQ (q,t)`. Finally, the
    ``Symqt.inject_shorthands()`` command makes the "usual"
    short names (as in Macdonald's book [Mac1995]_) for the bases
    of :math:`\operatorname{Sym}` available.
    The keyword `verbose` allows you to make the injection quiet
    (i.e., skip the output).
    Instead of using ``Symqt.inject_shorthands()``, you could
    also manually define the bases you need: e.g.,
    ``s = Symqt.s()`` for the Schur basis,
    ``p = Symqt.p()`` for the powersum basis, etc.

    Now that we have access to all the bases we need,
    we can start manipulating them.
    The basis elements are indexed by partitions::

        sage: s[101,14,13,11]
        s[101, 14, 13, 11]

    ::

        sage: e[3,2,1]
        e[3, 2, 1]

    ::

        sage: h[2,1],p[2,1],m[1]
        (h[2, 1], p[2, 1], m[1])

    .. NOTE::

        There are several ways to have Sage produce a basis element
        corresponding to a given partition::

             sage: p[2, 1, 1]
             p[2, 1, 1]
             sage: p[[2, 1, 1]]
             p[2, 1, 1]
             sage: p[Partition([2, 1, 1])]
             p[2, 1, 1]
             sage: p(Partition([2, 1, 1]))
             p[2, 1, 1]
             sage: p.basis()[Partition([2,1,1])] # the "conceptual" way
             p[2, 1, 1]
             sage: p[(i for i in [2, 1, 1])] # any iterable works here
             p[2, 1, 1]

        The first is the most convenient for direct use; the others
        are more suited for programming.
    
        In the special case of the empty partition, due to a limitation in
        Python syntax, one cannot use::

            sage: p[]       # todo: not implemented

        Please use instead::

            sage: p[[]]
            p[]

    .. NOTE::

        When elements are constructed using the ``p[something ]`` syntax,
        an error will be raised if the input cannot be interpreted as a partition.
        This is *not* the case when ``p.basis()`` is used::

            sage: p['something']
            Traceback (most recent call last):
            ...
            ValueError: all parts of 'something' should be nonnegative integers
            sage: p.basis()['something']
            p'something'

    .. WARNING::

        Do not confuse ``p[4]`` (the power-sum symmetric function
        :math:`p_4`) with ``p(4)`` (the scalar :math:`4`, regarded
        as a symmetric function in the ``p``-basis)::

            sage: p(4) # This is the constant 4, not the partition 4.
            4*p[]
            sage: p([4]) # This is the partition 4.
            p[4]

    Note that ``Sym`` is an *abstract* algebra.  This reflects the
    fact that it has multiple natural bases, and no basis is
    "privileged".  SageMath implements the algebra in each of
    these bases separately, calling it a *realization* of ``Sym``.
    For example, the realization ``Sym.m()`` is the algebra of
    symmetric functions in the basis of the monomial symmetric
    functions.  Elements can be easily combined, converted and
    compared between different realizations (e.g., you can write
    ``e[2] - h[1,1]`` or ``e[1] == h[1]``; see below for
    conversion).

    Thus, strictly speaking, ``p`` is not the power-sum *basis* in the
    linear-algebraic sense (the latter can be obtained by writing
    ``p.basis()``), but rather the ring of symmetric functions written in
    the power-sum basis (a realization of :math:`\operatorname{Sym}`) ::

        sage: p # The ring
        Symmetric Functions over Fraction Field of Multivariate
         Polynomial Ring in q, t over Rational Field in the
         powersum basis
        sage: p.basis() # The basis
        Lazy family (Term map from Partitions to Symmetric
         Functions over Fraction Field of Multivariate
         Polynomial Ring in q, t over Rational Field in the
         powersum basis(i))_{i in Partitions}
        sage: p.basis().keys() # The partitions indexing the basis
        Partitions

    (This explains the ``p.basis()[Partition([2,1,1])]`` syntax
    above.)

    Elements of ``p`` (that is, symmetric functions written
    in the power-sum basis) are linear combinations of basis
    elements::

        sage: p.an_element()
        2*p[] + 2*p[1] + 3*p[2]

    Likewise for the other bases::

        sage: (q+t)*s[2,1,1]
        (q+t)*s[2, 1, 1]
        sage: q*s[2,1,1] + t*s[2,1,1]
        (q+t)*s[2, 1, 1]
        sage: q*s[2,1,1] + t*s[2]
        t*s[2] + q*s[2, 1, 1]

    .. rubric:: The ring structure

    In the multiplicative bases (that is, :math:`e`, :math:`h` and :math:`p`),
    the product of two basis elements is obtained simply by combining the
    corresponding partitions (i.e., by concatenating them and sorting the
    result in nondecreasing order)::

        sage: p([2,1,1])*p([5,2]) == p([5,2,2,1,1])
        True
        sage: h[4, 2] * h[3, 2]
        h[4, 3, 2, 2]
        sage: e.one()
        e[]
        sage: p[2,1,1] + 2 * p[1] * (p[4] + p[2,1])
        3*p[2, 1, 1] + 2*p[4, 1]
        sage: (p.one() + 2 * p[3,1]) * p[4, 2]
        p[4, 2] + 2*p[4, 3, 2, 1]

    For the non-multiplicative bases, such as the Schur functions, products
    are expanded as linear combinations in the same (linear) basis
    according to combinatorial rules (such as the Littlewood-Richardson
    rule for the Schur basis)::

        sage: s([5])^2*s([1,1,1])
        s[5, 5, 1, 1, 1] + s[6, 4, 1, 1, 1] + 2*s[6, 5, 1, 1] + s[6, 6, 1]
         + s[7, 3, 1, 1, 1] + 2*s[7, 4, 1, 1] + s[7, 5, 1]
         + s[8, 2, 1, 1, 1] + 2*s[8, 3, 1, 1] + s[8, 4, 1]
         + s[9, 1, 1, 1, 1] + 2*s[9, 2, 1, 1] + s[9, 3, 1]
         + 2*s[10, 1, 1, 1] + s[10, 2, 1] + s[11, 1, 1]

        sage: m([3,1])*m([2,2])
        m[3, 2, 2, 1] + 2*m[3, 3, 2] + m[5, 2, 1] + m[5, 3]

    These calculations are relatively fast as illustrated in the following,
    showing only the length of the output rather than printing it out in all its glory::

        sage: len(s[10,5,5,3]*s[12,5,2])
        2986

    When we mix different bases, the result will, by default, be expressed in one of
    the bases, usually the first basis encountered in the expression::

        sage: s([2,1])*m([1,1])+p([2,2])
        s[1, 1, 1, 1] - s[2, 1, 1] + s[2, 1, 1, 1] + 2*s[2, 2] + s[2, 2, 1] - s[3, 1] + s[3, 1, 1] + s[3, 2] + s[4]

        sage: m([1,1])*s([2,1])+p([2,2])
        20*m[1, 1, 1, 1, 1] + 9*m[2, 1, 1, 1] + 2*m[2, 2] + 4*m[2, 2, 1] + 2*m[3, 1, 1] + m[3, 2] + m[4]

        sage: p([2,2])+m([1,1])*s([2,1])
        1/6*p[1, 1, 1, 1, 1] - 1/6*p[2, 1, 1, 1] + p[2, 2] - 1/6*p[3, 1, 1] + 1/6*p[3, 2]

    We can force the answer to be expressed with respect to a
    particular basis as follows::

        sage: s(m([1,1])*s([2,1])+p([2,2]))
        s[1, 1, 1, 1] - s[2, 1, 1] + s[2, 1, 1, 1] + 2*s[2, 2] + s[2, 2, 1] - s[3, 1] + s[3, 1, 1] + s[3, 2] + s[4]

    In general, if ``f`` is a symmetric function expressed in some
    basis, then ``p(f)`` will return ``f`` expressed in the power-sum
    basis; ``m(f)`` will return ``f`` expressed in the monomial basis;
    and so on. We will see more of this below in the "Change of bases"
    section.

    Internally, each of the bases ``s``, ``m``, ... of
    :math:`\operatorname{Sym}` is an algebra in its own right; all
    of these algebras are isomorphic, and the syntax ``p(f)`` or
    ``m(f)`` used above applies these isomorphisms to an element ``f``.

    .. NOTE::

        These isomorphisms need to be called in order to
        convert a symmetric function into a different basis; but they
        need not to be called to compare two symmetric functions written
        in different bases. Sage knows how to do it already::

            sage: m[1,1] - m[2] == 2*s[1,1] - s[2]
            True

    .. NOTE::

        The multiplicative identity of ``p`` is the
        basis element indexed by the empty partition::

            sage: p.one()
            p[]

        The same holds for all other bases of
        :math:`\operatorname{Sym}`::

            sage: p.one() == p[[]] == s.one() == s[[]] == m.one() == m[[]] == e.one() == e[[]] == h.one() == h[[]]
            True

    .. NOTE::

        Sage knows that :math:`\operatorname{Sym}` is not just an algebra,
        but has many other structures implemented in Sage::

            sage: Sym = SymmetricFunctions(QQ)
            sage: Sym
            Symmetric Functions over Rational Field
            sage: Sym.category()
            Join of Category of hopf algebras over Rational Field
                and Category of graded algebras over Rational Field
                and Category of monoids with realizations
                and Category of graded coalgebras over Rational Field
                and Category of coalgebras over Rational Field with realizations

        These extra structures are reflected on each realization
        (i.e., basis) of :math:`\operatorname{Sym}`.  For example, let us
        explore them on ``p``. We can ask for the mathematical properties
        of ``p``::

            sage: p.categories() # not tested
            [Category of graded bases of Symmetric Functions over Rational Field,
             Category of filtered bases of Symmetric Functions over Rational Field,
             Category of bases of Symmetric Functions over Rational Field,
             Category of graded hopf algebras with basis over Rational Field,
             ...]

        To start with, ``p`` is a graded algebra, the grading being induced
        by the size of the partitions.  As a graded algebra, it is therefore
        automatically a filtered algebra.  The remaining structures will be
        explained in the sections below.

    .. rubric:: Concrete symmetric functions

    Our above abstract symmetric functions represent (possibly very large)
    concrete multivariate polynomials that are invariant under any permutation
    of their variables. Simple examples include

    .. MATH::

        p_k(x_1,x_2,\ldots, x_n) = x_1^k + x_2^k + \cdots + x_n^k
              \quad (\text{for any } k > 0), \quad \text{ or }

    .. MATH::

        e_n(x_1,x_2,\ldots, x_n) = x_1 x_2 \cdots x_n.

    To expand a symmetric function into a concrete polynomial in the set of
    variables :math:`x_0, x_1, \dots, x_{n-1}`, one proceeds as follows::

        sage: p[3].expand(3)
        x0^3 + x1^3 + x2^3

    ::

        sage: h[3].expand(3)
        x0^3 + x0^2*x1 + x0*x1^2 + x1^3 + x0^2*x2 + x0*x1*x2 + x1^2*x2 + x0*x2^2 + x1*x2^2 + x2^3

    ::

        sage: e[3].expand(3)
        x0*x1*x2

    ::

        sage: s[3,1,1].expand(4)
        x0^3*x1*x2 + x0^2*x1^2*x2 + x0*x1^3*x2 + x0^2*x1*x2^2 + x0*x1^2*x2^2 + x0*x1*x2^3 + x0^3*x1*x3 + x0^2*x1^2*x3 + x0*x1^3*x3 + x0^3*x2*x3 + 3*x0^2*x1*x2*x3 + 3*x0*x1^2*x2*x3 + x1^3*x2*x3 + x0^2*x2^2*x3 + 3*x0*x1*x2^2*x3 + x1^2*x2^2*x3 + x0*x2^3*x3 + x1*x2^3*x3 + x0^2*x1*x3^2 + x0*x1^2*x3^2 + x0^2*x2*x3^2 + 3*x0*x1*x2*x3^2 + x1^2*x2*x3^2 + x0*x2^2*x3^2 + x1*x2^2*x3^2 + x0*x1*x3^3 + x0*x2*x3^3 + x1*x2*x3^3

    ::

        sage: m[3,1,1].expand(4)
        x0^3*x1*x2 + x0*x1^3*x2 + x0*x1*x2^3 + x0^3*x1*x3 + x0*x1^3*x3 + x0^3*x2*x3 + x1^3*x2*x3 + x0*x2^3*x3 + x1*x2^3*x3 + x0*x1*x3^3 + x0*x2*x3^3 + x1*x2*x3^3

    ::

        sage: f[3,1,1].expand(4)
        3*x0^5 + 2*x0^4*x1 + x0^3*x1^2 + x0^2*x1^3 + 2*x0*x1^4 + 3*x1^5 + 2*x0^4*x2 + x0^3*x1*x2 + x0*x1^3*x2 + 2*x1^4*x2 + x0^3*x2^2 + x1^3*x2^2 + x0^2*x2^3 + x0*x1*x2^3 + x1^2*x2^3 + 2*x0*x2^4 + 2*x1*x2^4 + 3*x2^5 + 2*x0^4*x3 + x0^3*x1*x3 + x0*x1^3*x3 + 2*x1^4*x3 + x0^3*x2*x3 + x1^3*x2*x3 + x0*x2^3*x3 + x1*x2^3*x3 + 2*x2^4*x3 + x0^3*x3^2 + x1^3*x3^2 + x2^3*x3^2 + x0^2*x3^3 + x0*x1*x3^3 + x1^2*x3^3 + x0*x2*x3^3 + x1*x2*x3^3 + x2^2*x3^3 + 2*x0*x3^4 + 2*x1*x3^4 + 2*x2*x3^4 + 3*x3^5

    One may also use any other set of variables via the optional "alphabet" keyword::

        sage: g = s[2,1]
        sage: g.expand(3, alphabet=['x','y','z'])
        x^2*y + x*y^2 + x^2*z + 2*x*y*z + y^2*z + x*z^2 + y*z^2

    .. TOPIC:: Exercise

        Let :math:`e_k(n) = e_k(x_0,x_1, \dots , x_{n-1})` and similarly for
        the complete homogeneous functions.
        Then we have the following recursion relations for :math:`n \geq 1`:

        .. MATH::

            e_k(n) = e_k(n-1) + x_n e_{k-1}(n-1), \\
            h_k(n) = h_k(n-1) + x_n h_{k-1}(n), \\
            e_k(0) = h_k(0) = \delta_{k,0},

        where :math:`\delta_{k,0}` is the Kronecker delta.

        Check these relations for :math:`k=3` and :math:`2 \leq n \leq 5`.

        ::

            sage: #write your solution here

    .. TOPIC:: Solution

        ::

            sage: k = 3
            sage: R = PolynomialRing(QQ, 'x', 5)
            sage: R.inject_variables()
            Defining x0, x1, x2, x3, x4
            sage: l = list(R.gens())
            sage: for xn, n in zip(l[1:], range(2,6)):
            ....:     f1 = e([k]).expand(n)
            ....:     print(f1)
            ....:     f2 = e([k]).expand(n-1,l[:n-1])+xn*(e([k-1]).expand(n-1,l[:n-1]))
            ....:     print(f2)
            ....:     g1 = h([k]).expand(n)
            ....:     print(g1)
            ....:     g2 = h([k]).expand(n-1,l[:n-1])+xn*(h([k-1]).expand(n,l[:n]))
            ....:     print(g2)
            0
            0
            x0^3 + x0^2*x1 + x0*x1^2 + x1^3
            x0^3 + x0^2*x1 + x0*x1^2 + x1^3
            x0*x1*x2
            x0*x1*x2
            x0^3 + x0^2*x1 + x0*x1^2 + x1^3 + x0^2*x2 + x0*x1*x2 + x1^2*x2 + x0*x2^2 + x1*x2^2 + x2^3
            x0^3 + x0^2*x1 + x0*x1^2 + x1^3 + x0^2*x2 + x0*x1*x2 + x1^2*x2 + x0*x2^2 + x1*x2^2 + x2^3
            x0*x1*x2 + x0*x1*x3 + x0*x2*x3 + x1*x2*x3
            x0*x1*x2 + x0*x1*x3 + x0*x2*x3 + x1*x2*x3
            x0^3 + x0^2*x1 + x0*x1^2 + x1^3 + x0^2*x2 + x0*x1*x2 + x1^2*x2 + x0*x2^2 + x1*x2^2 + x2^3 + x0^2*x3 + x0*x1*x3 + x1^2*x3 + x0*x2*x3 + x1*x2*x3 + x2^2*x3 + x0*x3^2 + x1*x3^2 + x2*x3^2 + x3^3
            x0^3 + x0^2*x1 + x0*x1^2 + x1^3 + x0^2*x2 + x0*x1*x2 + x1^2*x2 + x0*x2^2 + x1*x2^2 + x2^3 + x0^2*x3 + x0*x1*x3 + x1^2*x3 + x0*x2*x3 + x1*x2*x3 + x2^2*x3 + x0*x3^2 + x1*x3^2 + x2*x3^2 + x3^3
            x0*x1*x2 + x0*x1*x3 + x0*x2*x3 + x1*x2*x3 + x0*x1*x4 + x0*x2*x4 + x1*x2*x4 + x0*x3*x4 + x1*x3*x4 + x2*x3*x4
            x0*x1*x2 + x0*x1*x3 + x0*x2*x3 + x1*x2*x3 + x0*x1*x4 + x0*x2*x4 + x1*x2*x4 + x0*x3*x4 + x1*x3*x4 + x2*x3*x4
            x0^3 + x0^2*x1 + x0*x1^2 + x1^3 + x0^2*x2 + x0*x1*x2 + x1^2*x2 + x0*x2^2 + x1*x2^2 + x2^3 + x0^2*x3 + x0*x1*x3 + x1^2*x3 + x0*x2*x3 + x1*x2*x3 + x2^2*x3 + x0*x3^2 + x1*x3^2 + x2*x3^2 + x3^3 + x0^2*x4 + x0*x1*x4 + x1^2*x4 + x0*x2*x4 + x1*x2*x4 + x2^2*x4 + x0*x3*x4 + x1*x3*x4 + x2*x3*x4 + x3^2*x4 + x0*x4^2 + x1*x4^2 + x2*x4^2 + x3*x4^2 + x4^3
            x0^3 + x0^2*x1 + x0*x1^2 + x1^3 + x0^2*x2 + x0*x1*x2 + x1^2*x2 + x0*x2^2 + x1*x2^2 + x2^3 + x0^2*x3 + x0*x1*x3 + x1^2*x3 + x0*x2*x3 + x1*x2*x3 + x2^2*x3 + x0*x3^2 + x1*x3^2 + x2*x3^2 + x3^3 + x0^2*x4 + x0*x1*x4 + x1^2*x4 + x0*x2*x4 + x1*x2*x4 + x2^2*x4 + x0*x3*x4 + x1*x3*x4 + x2*x3*x4 + x3^2*x4 + x0*x4^2 + x1*x4^2 + x2*x4^2 + x3*x4^2 + x4^3

    .. rubric:: Convert a concrete symmetric polynomial into an abstract symmetric function

    Conversely, a "concrete" symmetric polynomial in finitely many
    variables may be written as a formal symmetric function in any
    chosen basis, using the ``from_polynomial()`` method on the
    basis. For example, let us expand :math:`p_2 + e_{2,1}`
    in three variables and recover it back::

        sage: pol1 = (p([2])+e([2,1])).expand(3)
        sage: pol1
        x0^2*x1 + x0*x1^2 + x0^2*x2 + 3*x0*x1*x2 + x1^2*x2 + x0*x2^2 + x1*x2^2 + x0^2 + x1^2 + x2^2
        sage: e.from_polynomial(pol1) == p([2])+e([2,1])
        True

    In general, expansion in finitely many variables is a lossy
    (i.e., non-injective) operation, so this process will not
    always recover the original symmetric function. For example::

        sage: pol1 = p[1, 1, 1].expand(2)
        sage: p.from_polynomial(pol1)
        3*p[2, 1] - 2*p[3]

    Let us write the symmetric polynomial
    :math:`\prod_{0 \leq j < k \leq 2} (x_k - x_j)^2`
    (the discriminant) in three variables `x_0, x_1, x_2` as a
    symmetric function in the elementary basis::

        sage: n = 3
        sage: R = PolynomialRing(FractionField(QQ['q','t']),'x',n)
        sage: X = R.gens()
        sage: R.inject_variables()
        Defining x0, x1, x2
        sage: Discr=mul(mul((X[k]-X[j])^2 for j in range(k)) for k in range(1,n))
        sage: Discr
        x0^4*x1^2 + (-2)*x0^3*x1^3 + x0^2*x1^4 + (-2)*x0^4*x1*x2 + 2*x0^3*x1^2*x2 + 2*x0^2*x1^3*x2 + (-2)*x0*x1^4*x2 + x0^4*x2^2 + 2*x0^3*x1*x2^2 + (-6)*x0^2*x1^2*x2^2 + 2*x0*x1^3*x2^2 + x1^4*x2^2 + (-2)*x0^3*x2^3 + 2*x0^2*x1*x2^3 + 2*x0*x1^2*x2^3 + (-2)*x1^3*x2^3 + x0^2*x2^4 + (-2)*x0*x1*x2^4 + x1^2*x2^4
        sage: e.from_polynomial(Discr)
        e[2, 2, 1, 1] - 4*e[2, 2, 2] - 4*e[3, 1, 1, 1] + 18*e[3, 2, 1] - 27*e[3, 3] - 8*e[4, 1, 1] + 24*e[4, 2]

    The ``pol`` input of the method ``from_polynomial()`` is assumed to
    lie in a polynomial ring over the same base ring as that used for the symmetric
    functions, which therefore has to be declared beforehand.

    ::

        sage: n = 3
        sage: R = PolynomialRing(FractionField(QQ['q','t']),'y',n)
        sage: R.inject_variables()
        Defining y0, y1, y2
        
    Here, we will work with three variables (:math:`y_0, y_1, y_2`).
    Finally, we can declare our polynomial and convert it into a symmetric function
    in the monomial basis for example. ::

        sage: pol2 = y0^2*y1 + y0*y1^2 + y0^2*y2 + 2*y0*y1*y2 + y1^2*y2 + y0*y2^2 + y1*y2^2
        sage: m.from_polynomial(pol2)
        2*m[1, 1, 1] + m[2, 1]

    In the preceding example, the base ring of polynomials is the same as the base
    ring of symmetric polynomials considered, as checked by the following::

        sage: print(s.base_ring())
        Fraction Field of Multivariate Polynomial Ring in q, t over Rational Field
        sage: print(pol2.base_ring())
        Fraction Field of Multivariate Polynomial Ring in q, t over Rational Field

    Thus a concrete symmetric polynomial over :math:`\QQ (q,t)` may be
    transformed into an abstract symmetric function in any basis::

        sage: R = PolynomialRing(QQ['q','t'],'y',3)
        sage: R.inject_variables()
        Defining y0, y1, y2
        sage: pol2 = 1+(y0*y1+y0*y2+y1*y2)*(q+t)+(y0*y1*y2)*(q*t)
        sage: s.from_polynomial(pol2)
        s[] + (q+t)*s[1, 1] + q*t*s[1, 1, 1]

    .. rubric:: Changes of bases

    Many calculations on symmetric functions involve
    a change of (linear) basis. As mentioned above, Sage
    will make such changes automatically if the user
    enters an expression that involves elements of
    different bases. For example, let us compute
    :math:`p_{22}+m_{11}s_{21}` in the elementary basis::

        sage: x = p([2,2]) + m([1,1]) * s([2,1]); x
        1/6*p[1, 1, 1, 1, 1] - 1/6*p[2, 1, 1, 1] + p[2, 2] - 1/6*p[3, 1, 1] + 1/6*p[3, 2]

    A change of basis can also be forced explicitly:
    If ``f`` is a symmetric function, and ``b`` is a basis
    of :math:`\operatorname{Sym}`, then ``b(f)`` will compute
    the expansion of ``f`` in the basis ``b``::

        sage: e(x)
        e[1, 1, 1, 1] - 4*e[2, 1, 1] + 4*e[2, 2] + e[2, 2, 1] - e[3, 2]
        sage: s(p[2,1])
        -s[1, 1, 1] + s[3]
        sage: m(p[3])
        m[3]
        sage: m(p[3,2])
        m[3, 2] + m[5]

    .. TOPIC:: Exercise

        Print all the Schur functions on partitions of size 5 and convert them into the elementary basis.

    .. TOPIC:: Solution

    ::

        sage: for mu in Partitions(5):
        ....:     print(s(mu))
        ....:     print(e(s(mu)))
        s[5]
        e[1, 1, 1, 1, 1] - 4*e[2, 1, 1, 1] + 3*e[2, 2, 1] + 3*e[3, 1, 1] - 2*e[3, 2] - 2*e[4, 1] + e[5]
        s[4, 1]
        e[2, 1, 1, 1] - 2*e[2, 2, 1] - e[3, 1, 1] + 2*e[3, 2] + e[4, 1] - e[5]
        s[3, 2]
        e[2, 2, 1] - e[3, 1, 1] - e[3, 2] + e[4, 1]
        s[3, 1, 1]
        e[3, 1, 1] - e[3, 2] - e[4, 1] + e[5]
        s[2, 2, 1]
        e[3, 2] - e[4, 1]
        s[2, 1, 1, 1]
        e[4, 1] - e[5]
        s[1, 1, 1, 1, 1]
        e[5]


    .. TOPIC:: Exercise

        Compute the sum of the complete homogeneous functions on partitions of size 4 in the power sum basis.

    .. TOPIC:: Solution

    ::

        sage: p(sum(h(mu) for mu in Partitions(4)))
        47/24*p[1, 1, 1, 1] + 7/4*p[2, 1, 1] + 3/8*p[2, 2] + 2/3*p[3, 1] + 1/4*p[4]



    .. TOPIC:: Exercise

        It is well-known that :math:`h_n(X) = \sum \limits_{\mu \vdash n} \dfrac{p_{\mu}(x)}{z_{\mu}}` for every integer :math:`n \geq 0`.
        Verify this result for :math:`n \in \{1,2,3,4\}`.

        Note that there exists a function ``zee()`` which takes a partition
        :math:`\mu` and returns the value of :math:`z_{\mu}`. To use this
        function, you should import it from ``sage.combinat.sf.sfa``::

            sage: from sage.combinat.sf.sfa import zee
            sage: zee([4,4,2,1])
            64

    .. TOPIC:: Solution

    ::

        sage: for n in range(1,5):
        ....:     print(p(h([n])) == sum(p(mu)/zee(mu) for mu in Partitions(n)))
        True
        True
        True
        True

    ::

    .. rubric:: Other well-known bases

    Other bases of symmetric functions are implemented in Sage:

    - The Hall-Littlewood basis
    - The Jack basis
    - The orthogonal basis
    - The symplectic basis
    - The Witt basis
    - The zonal basis
    - The Hecke character basis

    Various bases related to the well-known Macdonald symmetric functions
    are also implemented in Sage.  For more details, you can consult the
    following Sage reference:
    http://doc.sagemath.org/html/en/reference/combinat/sage/combinat/sf/macdonald.html

    For more information, see the documentation of the individual bases.

    We briefly demonstrate how to access these bases. For more information, see
    the documentation of the individual bases.

    The *Jack polynomials* can be obtained as::

        sage: Sym = SymmetricFunctions(FractionField(QQ['t']))
        sage: Jack = Sym.jack()
        sage: P = Jack.P(); J = Jack.J(); Q = Jack.Q()
        sage: J(P[2,1])
        (1/(t+2))*JackJ[2, 1]

    The parameter `t` can be specialized as follows::

        sage: Sym = SymmetricFunctions(QQ)
        sage: Jack = Sym.jack(t = 1)
        sage: P = Jack.P(); J = Jack.J(); Q = Jack.Q()
        sage: J(P[2,1])
        1/3*JackJ[2, 1]

    Similarly one can access the Hall-Littlewood polynomials.
    The Macdonald polynomial bases can be accessed via
    ``Sym.macdonald().H()``, ``Sym.macdonald().Ht()``,
    ``Sym.macdonald().P()``, ``Sym.macdonald().J()``,
    ``Sym.macdonald().Q()``, ``Sym.macdonald().S()``::

        sage: Sym = SymmetricFunctions(FractionField(QQ['q','t']))
        sage: Mcd = Sym.macdonald()
        sage: P = Mcd.P(); J = Mcd.J(); Q = Mcd.Q()
        sage: J(P[2,1])
        (1/(-q*t^4+2*q*t^3-q*t^2+t^2-2*t+1))*McdJ[2, 1]

    Non-generic values of :math:`q` and :math:`t` can be
    specified as parameters in ``Sym.macdonald()``::

        sage: Sym = SymmetricFunctions(QQ)
        sage: Mcd12 = Sym.macdonald(q=1, t=2)
        sage: P = Mcd12.P(); s = Sym.s()
        sage: s(P[2])
        s[1, 1] + s[2]

    Here are some examples involving the "combinatorial"
    Macdonald symmetric functions.  These are eigenfunctions
    of the operator :math:`\nabla` discussed further below.

    ::

        sage: H = Symqt.macdonald().Ht(); s = Symqt.s()

    ::

        sage: s(H([2,1]))
        q*t*s[1, 1, 1] + (q+t)*s[2, 1] + s[3]
        sage: H(s[2,1])
        ((-q)/(-q*t^2+t^3+q^2-q*t))*McdHt[1, 1, 1] + ((q^2+q*t+t^2)/(-q^2*t^2+q^3+t^3-q*t))*McdHt[2, 1] + (t/(-q^3+q^2*t+q*t-t^2))*McdHt[3]


    ::

        sage: [H(mu).nabla() for mu in Partitions(4)]
        [q^6*McdHt[4],
         q^3*t*McdHt[3, 1],
         q^2*t^2*McdHt[2, 2],
         q*t^3*McdHt[2, 1, 1],
         t^6*McdHt[1, 1, 1, 1]]

    We can also construct the `\bar{q}` basis that can be used
    to determine character tables for Hecke algebras (with quadratic
    relation `T_i^2 = (1-q) T_i + q`)::

        sage: Sym = SymmetricFunctions(ZZ['q'].fraction_field())
        sage: qbar = Sym.hecke_character()
        sage: s = Sym.s()
        sage: s(qbar[2,1])
        -s[1, 1, 1] + (q-1)*s[2, 1] + q*s[3]


    .. rubric:: More basic commands on symmetric functions

    We can see that the terms of an expression are always listed in a specific order on the partitions. This order can be changed.

    First, the method  ``get_print_style()`` applied to a basis hands the order used on partitions for this basis. Then, with  ``set_print_style()`` we can set another printing order. The possible orders are:

    -  ``lex``   : lexicographic order.
    -  ``length``   : by partition length, and then by lexicographic order for partitions of same length.
    -  ``maximal_part`` :  by value of the biggest part of the partition, then by length, then by lexicographic order.

    ::

        sage: s.get_print_style()
        'lex'

    ::

        sage: s.set_print_style('lex')
        sage: s(p[4,1,1])
        -s[1, 1, 1, 1, 1, 1] - s[2, 1, 1, 1, 1] + s[2, 2, 1, 1] + s[2, 2, 2] - s[3, 3] - s[4, 2] + s[5, 1] + s[6]


    ::

        sage: s.set_print_style('length')
        sage: s(p[4,1,1])
        s[6] - s[3, 3] - s[4, 2] + s[5, 1] + s[2, 2, 2] + s[2, 2, 1, 1] - s[2, 1, 1, 1, 1] - s[1, 1, 1, 1, 1, 1]

    ::

        sage: s.get_print_style()
        'length'

    ::

        sage: s.set_print_style('maximal_part')
        sage: s(p[4,1,1])
        -s[1, 1, 1, 1, 1, 1] + s[2, 2, 2] + s[2, 2, 1, 1] - s[2, 1, 1, 1, 1] - s[3, 3] - s[4, 2] + s[5, 1] + s[6]
        sage: s.set_print_style('lex')

    The method ``coefficient()`` returns the coefficient
    associated to a given partition. ::

        sage: f = s[5,2,2,1]
        sage: e(f)
        e[4, 3, 1, 1, 1] - 2*e[4, 3, 2, 1] + e[4, 3, 3] - e[4, 4, 1, 1] + e[4, 4, 2] - e[5, 2, 1, 1, 1] + 2*e[5, 2, 2, 1] - e[5, 3, 2] + e[5, 4, 1] + e[6, 2, 1, 1] - e[6, 2, 2] - e[6, 4] - e[7, 2, 1] + e[8, 2]
        sage: e(f).coefficient([4,3,2,1])
        -2

    The method ``degree()`` gives the degree of a symmetric
    function. ::

        sage: f.degree()
        10

    Finally, the method ``support()`` returns a list of
    all partitions that appear with nonzero coefficient
    in a given symmetric function.
    The result will depend on the basis of expression of
    the function. In the following example, we also use
    the method ``sorted()`` to get an ordered list. ::

        sage: f.support()
        [[5, 2, 2, 1]]
        sage: sorted(h(f).support())
        [[5, 2, 2, 1], [5, 3, 1, 1], [5, 3, 2], [5, 4, 1], [6, 2, 1, 1], [6, 3, 1], [6, 4], [7, 1, 1, 1], [7, 2, 1], [8, 1, 1], [8, 2]]

    .. rubric:: Transformations of symmetric functions

    There are many methods in Sage which make it easy to manipulate symmetric
    functions.  For example, if we have some function that acts on partitions
    (say, conjugation), it is a simple matter to apply it to the support of a
    symmetric function.  Here is an example::

        sage: conj = lambda mu: mu.conjugate()
        sage: f = h[4] + 2*h[3,1]
        sage: f.map_support(conj)
        h[1, 1, 1, 1] + 2*h[2, 1, 1]

    We can also easily modify the coefficients::

        sage: def foo(mu, coeff): return mu.conjugate(), -coeff
        sage: f.map_item(foo)
        -h[1, 1, 1, 1] - 2*h[2, 1, 1]

    See also ``map_coefficients``.

    There are also methods for building symmetric functions
    directly as sums or linear combinations of basis elements::

        sage: s.sum_of_monomials(mu for mu in Partitions(3))
        s[1, 1, 1] + s[2, 1] + s[3]
        sage: s.sum_of_monomials(Partitions(3))
        s[1, 1, 1] + s[2, 1] + s[3]
        sage: s.sum_of_terms( (mu, mu[0]) for mu in Partitions(3))
        s[1, 1, 1] + 2*s[2, 1] + 3*s[3]

    These are the preferred ways to build elements within a program;
    the result will usually be faster than using :func:`sum`. It also
    guarantees that empty sums yield the zero of ``s`` (see also
    ``s.sum``).

    Note also that it is a good idea to use::

        sage: s.one()
        s[]
        sage: s.zero()
        0

    instead of ``s(1)`` and ``s(0)`` within programs where speed
    is important (in order to prevent unnecessary coercions) or
    over exotic base rings whose `0` and `1` are not the
    corresponding integers.

    .. rubric:: Different base rings

    Depending on the base ring, the different realizations of the
    symmetric function algebra may not span the same space::

        sage: SZ = SymmetricFunctions(ZZ)
        sage: p = SZ.power(); s = SZ.schur()
        sage: p(s[1,1,1])
        Traceback (most recent call last):
        ...
        TypeError: no conversion of this rational to integer

    Because of this, some functions may not behave as expected
    when working over the integers, even though they make
    mathematical sense::

        sage: s[1,1,1].plethysm(s[1,1,1])
        Traceback (most recent call last):
        ...
        TypeError: no conversion of this rational to integer

    It is possible to work over different base rings simultaneously::

        sage: s = SymmetricFunctions(QQ).schur()
        sage: p = SymmetricFunctions(QQ).power()
        sage: sz = SymmetricFunctions(ZZ).schur(); sz._prefix = 'sz'
        sage: pz = SymmetricFunctions(ZZ).power(); pz._prefix = 'pz'
        sage: p(sz[1,1,1])
        1/6*p[1, 1, 1] - 1/2*p[2, 1] + 1/3*p[3]
        sage: sz( 1/6*p[1, 1, 1] - 1/2*p[2, 1] + 1/3*p[3] )
        sz[1, 1, 1]

    As shown in this example, if you are working over multiple base rings
    simultaneously, it is a good idea to change the prefix in some cases, so that
    you can tell from the output which realization your result is in.

    Let us change the notation back for the remainder of this tutorial::

        sage: sz._prefix = 's'
        sage: pz._prefix = 'p'

    One can also use the Sage standard renaming idiom to get shorter outputs::

        sage: Sym = SymmetricFunctions(QQ)
        sage: Sym.rename("Sym")
        sage: Sym
        Sym
        sage: Sym.rename() # rename it back
        sage: Sym
        Symmetric Functions over Rational Field

    .. rubric:: The omega involution

    The :math:`\omega` involution is the linear transformation
    :math:`\operatorname{Sym} \to \operatorname{Sym}` that sends each
    :math:`e_\lambda` onto :math:`h_{\lambda}`.

    :: 

        sage: f = s[2]^2; f
        s[2, 2] + s[3, 1] + s[4]
        sage: h(f)
        h[2, 2]
        sage: e(f.omega())
        e[2, 2]
        sage: [(s(mu),s(mu).omega()) for mu in Partitions(5)]
        [(s[5], s[1, 1, 1, 1, 1]),
         (s[4, 1], s[2, 1, 1, 1]),
         (s[3, 2], s[2, 2, 1]),
         (s[3, 1, 1], s[3, 1, 1]),
         (s[2, 2, 1], s[3, 2]),
         (s[2, 1, 1, 1], s[4, 1]),
         (s[1, 1, 1, 1, 1], s[5])]

    For more methods than discussed here, create a symmetric function as
    above, and use ``f.<tab>``.

    .. rubric:: Scalar Products

    The Hall scalar product is the standard scalar product on the algebra of
    symmetric functions. It makes the Schur functions into an orthonormal basis. ::

        sage: s[3,1].scalar(s[3,1])
        1
        sage: s[3,1].scalar(s[4])
        0
        sage: f.scalar(f)
        3

    The Hall scalar product of :math:`p_{\mu}` and :math:`p_{\lambda}`
    is :math:`z_{\mu}` if :math:`\mu = \lambda` and zero otherwise.
    In formula,

    .. MATH:: \langle p_\mu,p_\lambda\rangle = z_\mu\,\delta_{\mu,\lambda}

    Equivalently, in terms of matrices, this says that

    .. MATH:: \left(\langle p_\mu,p_\lambda/z_\lambda\rangle\right)_{\mu,\lambda} = \operatorname{Id}

    Thus, we get ::

            sage: p([2,2,1]).scalar(p([2,2,1]))
            8
            sage: Matrix([[p(mu).scalar(p(nu)/zee(mu)) for nu in Partitions(5)] for mu in Partitions(5)])
            [1 0 0 0 0 0 0]
            [0 1 0 0 0 0 0]
            [0 0 1 0 0 0 0]
            [0 0 0 1 0 0 0]
            [0 0 0 0 1 0 0]
            [0 0 0 0 0 1 0]
            [0 0 0 0 0 0 1]


    .. rubric :: Other scalar products, such as the :math:`q,t`-scalar product

    One may specify, as an optional argument named ``zee``
    to the ``scalar`` method, a function on partitions giving
    the value for the scalar product between :math:`p_{\mu}`
    and :math:`p_{\mu}`. Power sums remain orthogonal for the
    resulting scalar product. By default, this value is
    :math:`z_{\mu}`, but other interesting cases include:

    .. MATH:: \langle p_{\mu},p_{\mu}\rangle_{q,t} = z_\mu\,\prod_i\frac{1-q^{\mu_i}}{1-t^{\mu_i}}.

    This is already refined as ``scalar_qt()``::

        sage: Matrix([[p(mu).scalar_qt(p(nu)/zee(mu)) for nu in Partitions(3)] for mu in Partitions(3)])
        [                            (-q^3 + 1)/(-t^3 + 1)                                                 0                                                 0]
        [                                                0           (q^3 - q^2 - q + 1)/(t^3 - t^2 - t + 1)                                                 0]
        [                                                0                                                 0 (-q^3 + 3*q^2 - 3*q + 1)/(-t^3 + 3*t^2 - 3*t + 1)]


     .. rubric:: Schur Positivity


    Sometimes one wants to check whether a given symmetric
    function is Schur-positive or not. In our current
    setup, this means that when the function is expanded
    in the Schur basis :math:`(s_\lambda)`, the
    coefficients are polynomials in :math:`\NN [q,t]`.
    If we worked over :math:`\QQ` instead, it would simply
    mean that the coefficients are nonnegative integers.
    The following function returns ``True`` if the given
    symmetric function is Schur-positive and ``False`` if
    not. ::

        sage: f = s([4,1])+s([3,2])
        sage: print(f.is_schur_positive())
        True
        sage: g = s([4,1])-s([3,2])
        sage: print(g.is_schur_positive())
        False

    For example, we can verify the well-known Schur
    positivity of products of Schur functions::

        sage: for mu in Partitions(2):
        ....:     for nu in Partitions(3):
        ....:         if (s(mu)*s(nu)).is_schur_positive():
        ....:             print('The product of ', s(mu),' and ',s(nu),' is Schur positive.')
        ....:         else:
        ....:             print('The product of ', s(mu),' and ',s(nu),' is not Schur positive.')
        The product of  s[2]  and  s[3]  is Schur positive.
        The product of  s[2]  and  s[2, 1]  is Schur positive.
        The product of  s[2]  and  s[1, 1, 1]  is Schur positive.
        The product of  s[1, 1]  and  s[3]  is Schur positive.
        The product of  s[1, 1]  and  s[2, 1]  is Schur positive.
        The product of  s[1, 1]  and  s[1, 1, 1]  is Schur positive.


    .. TOPIC:: Exercise

        Its representation-theoretic interpretation implies that
        :math:`\nabla (e_n)` is Schur-positive. Check this for
        :math:`1 \leq n \leq 6`.

    .. TOPIC:: Solution

    ::

        sage: for n in range(1,7):
        ....:     print(e([n]).nabla().is_schur_positive())
        True
        True
        True
        True
        True
        True


    Schur positivity is a rare phenomenon in general, but symmetric functions that come from representation theory are Schur-positive, as are those that
    represent effective cohomology classes on a Grassmannian. One can show that the probability of a degree-:math:`n` monomial positive symmetric function
    being Schur positive is

    .. MATH:: \prod_{\mu\vdash n}\frac{1}{k_\mu},\qquad {\rm where}\qquad k_\mu:=\sum_{\nu\vdash n} K_{\mu,\nu},

    with :math:`K_{\mu,\nu}` the **Kostka numbers**. Recall that these occur also in the monomial expansion of Schur functions:

    .. MATH:: s_\mu=\sum_\nu K_{\mu,\nu}\, m_\nu.

    For instance, we have ::

        sage: m(s[3,2])
        5*m[1, 1, 1, 1, 1] + 3*m[2, 1, 1, 1] + 2*m[2, 2, 1] + m[3, 1, 1] + m[3, 2]

    Hence, defining ::

        sage: def K(mu,nu):
        ....:     return s(mu).scalar(h(nu))

    the above expression may indeed be seen as ::

        sage: add(K([3,2],nu)*m(nu) for nu in Partitions(5))
        5*m[1, 1, 1, 1, 1] + 3*m[2, 1, 1, 1] + 2*m[2, 2, 1] + m[3, 1, 1] + m[3, 2]

    Now, we set ::

        sage: def k(mu):
        ....:     n=add(j for j in mu)
        ....:     return add(K(mu,nu) for nu in Partitions(n))

    so that the above probability is computed by the function ::

        sage: def prob_Schur_positive(n): return 1/mul(k(mu) for mu in Partitions(n))

    One can then illustrate how very rare Schur-positivity
    is, as a function of the degree::

        sage: [prob_Schur_positive(n) for n in range(1,8)]
        [1, 1/2, 1/9, 1/560, 1/480480, 1/1027458432000, 1/2465474364698304960000]

    .. rubric:: Plethysm

    **Plethysm** of symmetric functions is a binary operation (defined
    when the base ring is a :math:`\QQ`-algebra, and in some
    other situations) on :math:`\operatorname{Sym}` that is
    characterized by the properties
    - :math:`(f_1+f_2)\circ g = (f_1\circ g) + (f_2\circ g)`,
    - :math:`(f_1\cdot f_2)\circ g = (f_1\circ g) \cdot (f_2\circ g)`,
    - :math:`p_k\circ (g_1+g_2) = (p_k\circ g_1) + (p_k\circ g_2)`,
    - :math:`p_k\circ (g_1\cdot g_2) = (p_k\circ g_1) \cdot (p_k\circ g_2)`,
    - :math:`p_k\circ p_n = p_{kn}`,
    - :math:`p_k\circ x = x^k`, if :math:`x` is a **variable**,
    - :math:`p_k\circ c = c`, if :math:`c` is a **constant**.

    It corresponds to composition of representations of the general
    linear group.  See [EnumComb2]_ Chapter 7, Appendix 2 for details.

    As its name strongly suggests, the ``plethysm()`` function computes
    the plethysm :math:`f\circ g` of two symmetric functions
    :math:`f` and :math:`g`. ::

        sage: s[2].plethysm(s[2])
        s[2, 2] + s[4]

    One may specify a list of Sage-variables to be treated as **variables**
    in a plethysm, using the option ``include=[v1,v2,...,vk]``, and/or a list
    of Sage-variables to be considered as **constants**, using the option
    ``exclude=[c1,c2,...,ck]``. Here are some examples. ::

        sage: F = QQ['q','t'].fraction_field()
        sage: F.inject_variables()
        Defining q, t
        sage: s = SymmetricFunctions(F).s()
        sage: p = SymmetricFunctions(F).p()
        sage: p([3,2]).plethysm(h([3,1]))
        1/36*p[3, 3, 3, 3, 2, 2, 2, 2] + 1/12*p[4, 3, 3, 3, 3, 2, 2] + 1/12*p[6, 3, 3, 2, 2, 2, 2] + 1/18*p[6, 3, 3, 3, 3, 2] + 1/4*p[6, 4, 3, 3, 2, 2] + 1/6*p[6, 6, 3, 3, 2] + 1/18*p[9, 3, 2, 2, 2, 2] + 1/6*p[9, 4, 3, 2, 2] + 1/9*p[9, 6, 3, 2]
        sage: g = p([1]) + t*s([2,1])
        sage: p([2]).plethysm(g,include=[t])
        p[2] + 1/3*t^2*p[2, 2, 2] + (-1/3*t^2)*p[6]
        sage: p([2]).plethysm(g,exclude=[t])
        p[2] + 1/3*t*p[2, 2, 2] + (-1/3*t)*p[6]

    If the arguments ``include`` and ``exclude`` are not provided,
    then Sage automatically considers the designated generators of
    the base ring (e.g., the variables in a polynomial ring) to be
    the variables::

        sage: R.<t> = QQ[]; s = SymmetricFunctions(R).schur()
        sage: s[2]( (1-t)*s[1] )
        (t^2-t)*s[1, 1] + (-t+1)*s[2]

    It is customary to also write :math:`f[g]` for :math:`f\circ g` in
    mathematical texts; likewise, Sage uses the shorthand notation
    ``f(g)`` for this (with parentheses, not square brackets). For
    instance, the plethysm :math:`s_4\circ s_2` may also be computed
    as follows::

        sage: s[4](s[2])
        s[2, 2, 2, 2] + s[4, 2, 2] + s[4, 4] + s[6, 2] + s[8]

    To have nice expressions for plethystic substitutions, one may set aliases
    for the symmetric function on the empty partition
    (i.e. :math:`s_0, m_0, \dots`, all equal to the constant 1), and the
    symmetric function (unique up to a scalar) of degree 1.

    ::

        sage: R = Frac(PolynomialRing(QQ, 'q'))
        sage: q = R.gen()
        sage: s = SymmetricFunctions(R).schur()
        sage: One = s.one()
        sage: X = s[1]

    ::

        sage: s[3](s[4](One*(1+q)))
        (q^12+q^11+2*q^10+3*q^9+4*q^8+4*q^7+5*q^6+4*q^5+4*q^4+3*q^3+2*q^2+q+1)*s[]


    One should compare this with

    ::

        sage: q_binomial(7,3)
        q^12 + q^11 + 2*q^10 + 3*q^9 + 4*q^8 + 4*q^7 + 5*q^6 + 4*q^5 + 4*q^4 + 3*q^3 + 2*q^2 + q + 1


    ::

        sage: s[4](X*(1+q))
        q^2*s[2, 2] + (q^3+q^2+q)*s[3, 1] + (q^4+q^3+q^2+q+1)*s[4]


    ::

        sage: s[4](X/(1-q)).map_coefficients(factor)
        ((q-1)^-4*(q+1)^-2*q^6*(q^2+1)^-1*(q^2+q+1)^-1)*s[1, 1, 1, 1] + ((q-1)^-4*(q+1)^-2*q^3*(q^2+1)^-1)*s[2, 1, 1] + ((q-1)^-4*(q+1)^-2*q^2*(q^2+q+1)^-1)*s[2, 2] + ((q-1)^-4*(q+1)^-2*q*(q^2+1)^-1)*s[3, 1] + ((q-1)^-4*(q+1)^-2*(q^2+1)^-1*(q^2+q+1)^-1)*s[4]

    ::

        sage: s[3](s[4])-s[2](s[6])
        s[4, 4, 4] + s[6, 4, 2] + s[7, 4, 1] + s[8, 2, 2] + s[9, 3]


    These calculations suggest that we have the following
    polynomial with positive coefficients::

        sage: q_binomial(7,3)-q_binomial(8,2)
        q^9 + q^8 + q^7 + q^6 + q^5 + q^4 + q^3

    .. rubric:: Some interesting operators on symmetric functions

    Operators on symmetric functions may be found in Sage. Among these,
    the **nabla operator** is characterized as having the combinatorial
    Macdonald symmetric functions :math:`H_{\mu}=H_{\mu}(\mathbf{x};q,t)`
    as eigenfunctions:

    .. MATH:: \nabla H_{\mu} = t^{n(\mu)} q^{n(\mu')} H_{\mu},

    where :math:`\mu` is a partition, :math:`\mu'` its conjugate, and
    :math:`n(\mu)` is defined as :math:`\sum_i (i-1)\mu_i`.
    This operator :math:`\nabla` is thus defined on the symmetric
    functions with coefficients in the fraction field of
    :math:`\QQ [q,t]`, as declared above.

    It has been shown by Mark Haiman that :math:`\nabla(e_n)` is the Frobenius transform
    of the bigraded character of the :math:`S_n`-module of diagonal harmonic
    polynomials. Recall that the Frobenius transform maps irreducible modules onto Schur
    functions. ::

        sage: F = QQ['q','t'].fraction_field()
        sage: F.inject_variables()
        Defining q, t
        sage: Symqt = SymmetricFunctions(F)
        sage: Symqt.inject_shorthands(verbose=False)
        sage: s(e[3].nabla())
        (q^3+q^2*t+q*t^2+t^3+q*t)*s[1, 1, 1] + (q^2+q*t+t^2+q+t)*s[2, 1] + s[3]


    The dimension of this module is :math:`(n+1)^{n-1}`, and the dimension of its alternating component (see exercise below) is the Catalan number :math:`C_n=\frac{1}{n+1}\binom{2n}{n}`. The bigraded version has many other interesting properties. ::

        sage: Hilb_qt=s(e[3].nabla()).scalar(p[1]^3); Hilb_qt
        q^3 + q^2*t + q*t^2 + t^3 + 2*q^2 + 3*q*t + 2*t^2 + 2*q + 2*t + 1
        sage: Hilb_qt.substitute({q:1,t:1})
        16


    There are also interesting conjectures on the effect of :math:`\nabla` on Schur functions.

    ::

        sage: (-s([2,2,1])).nabla()
            (q^6*t^3+q^5*t^4+q^4*t^5+q^3*t^6)*s[1, 1, 1, 1, 1] + (q^6*t^2+2*q^5*t^3+2*q^4*t^4+2*q^3*t^5+q^2*t^6+q^4*t^3+q^3*t^4)*s[2, 1, 1, 1] + (q^5*t^2+2*q^4*t^3+2*q^3*t^4+q^2*t^5)*s[2, 2, 1] + (q^5*t^2+q^4*t^3+q^3*t^4+q^2*t^5+q^4*t^2+2*q^3*t^3+q^2*t^4)*s[3, 1, 1] + (q^4*t^2+q^3*t^3+q^2*t^4)*s[3, 2] + (q^3*t^2+q^2*t^3)*s[4, 1]

    .. TOPIC:: Exercise

        We have the following relation between :math:`\nabla (e_n)` and the q,t-Catalan numbers:

        .. MATH:: C_n(q,t) = \langle \nabla e_n , e_n \rangle.

        Check this relation for :math:`1 \leq n \leq 5`

        Note that the :math:`n`-th q,t-Catalan number can be computed
        using the command ``qt_catalan_number(n)``, which has to be
        imported from ``sage.combinat.q_analogues`` if not done yet.

        ::

            sage: from sage.combinat.q_analogues import *
            sage: for n in range (1,6):
            ....:     print((n,qt_catalan_number(n)))
            (1, 1)
            (2, q + t)
            (3, q^3 + q^2*t + q*t^2 + t^3 + q*t)
            (4, q^6 + q^5*t + q^4*t^2 + q^3*t^3 + q^2*t^4 + q*t^5 + t^6 + q^4*t + q^3*t^2 + q^2*t^3 + q*t^4 + q^3*t + q^2*t^2 + q*t^3)
            (5, q^10 + q^9*t + q^8*t^2 + q^7*t^3 + q^6*t^4 + q^5*t^5 + q^4*t^6 + q^3*t^7 + q^2*t^8 + q*t^9 + t^10 + q^8*t + q^7*t^2 + q^6*t^3 + q^5*t^4 + q^4*t^5 + q^3*t^6 + q^2*t^7 + q*t^8 + q^7*t + 2*q^6*t^2 + 2*q^5*t^3 + 2*q^4*t^4 + 2*q^3*t^5 + 2*q^2*t^6 + q*t^7 + q^6*t + q^5*t^2 + 2*q^4*t^3 + 2*q^3*t^4 + q^2*t^5 + q*t^6 + q^4*t^2 + q^3*t^3 + q^2*t^4)
            sage: for n in range (1,6):
            ....:     print((n,e([n]).nabla().scalar(e([n])).substitute({q:1,t:1})))
            (1, 1)
            (2, 2)
            (3, 5)
            (4, 14)
            (5, 42)

        ::

            sage: for n in range (1,6):
            ....:     print((n,factor(e([n]).nabla().scalar(e([n])).substitute({t:1/q}))))
            (1, 1)
            (2, q^-1 * (q^2 + 1))
            (3, q^-3 * (q^2 - q + 1) * (q^4 + q^3 + q^2 + q + 1))
            (4, q^-6 * (q^2 - q + 1) * (q^4 + 1) * (q^6 + q^5 + q^4 + q^3 + q^2 + q + 1))
            (5, q^-10 * (q^4 + 1) * (q^4 - q^3 + q^2 - q + 1) * (q^6 + q^3 + 1) * (q^6 + q^5 + q^4 + q^3 + q^2 + q + 1))


    .. TOPIC:: Solution

    ::

        sage: for n in range (1,6):
        ....:     print(e([n]).nabla().scalar(e([n])) == qt_catalan_number(n))
        True
        True
        True
        True
        True

    .. rubric:: `k`-Schur functions

    The `k`-Schur functions live in the `k`-bounded subspace of the ring of
    symmetric functions. It is possible to compute in the `k`-bounded subspace
    directly::

        sage: Sym = SymmetricFunctions(QQ)
        sage: ks = Sym.kschur(3,1)
        sage: f = ks[2,1]*ks[2,1]; f
        ks3[2, 2, 1, 1] + ks3[2, 2, 2] + ks3[3, 1, 1, 1]

    or to lift to the ring of symmetric functions::

        sage: f.lift()
        s[2, 2, 1, 1] + s[2, 2, 2] + s[3, 1, 1, 1] + 2*s[3, 2, 1] + s[3, 3] + s[4, 1, 1] + s[4, 2]

    However, it is not always possible to convert a symmetric function to the `k`-bounded subspace::

        sage: s = Sym.schur()
        sage: ks(s[2,1,1])
        Traceback (most recent call last):
        ...
        ValueError: s[2, 1, 1] is not in the image

    The `k`-Schur functions are more generally defined with a parameter `t` and they are
    a basis of the subspace spanned by the Hall-Littlewood ``Qp`` symmetric functions
    indexed by partitions whose first part is less than or equal to `k`::

        sage: Sym = SymmetricFunctions(QQ['t'].fraction_field())
        sage: SymS3 = Sym.kBoundedSubspace(3) # default t='t'
        sage: ks = SymS3.kschur()
        sage: Qp = Sym.hall_littlewood().Qp()
        sage: ks(Qp[2,1,1,1])
        ks3[2, 1, 1, 1] + (t^2+t)*ks3[2, 2, 1] + (t^3+t^2)*ks3[3, 1, 1] + t^4*ks3[3, 2]

    The subspace spanned by the `k`-Schur functions with a parameter `t` is not known
    to form a natural algebra.  However it is known that the product of a `k`-Schur
    function and an `\ell`-Schur function is in the linear span of the `k+\ell`-Schur
    functions::

        sage: ks(ks[2,1]*ks[1,1])
        Traceback (most recent call last):
        ...
        ValueError: s[2, 1, 1, 1] + s[2, 2, 1] + s[3, 1, 1] + s[3, 2] is not in the image
        sage: ks[2,1]*ks[1,1]
        s[2, 1, 1, 1] + s[2, 2, 1] + s[3, 1, 1] + s[3, 2]
        sage: ks6 = Sym.kBoundedSubspace(6).kschur()
        sage: ks6(ks[3,1,1]*ks[3])
        ks6[3, 3, 1, 1] + ks6[4, 2, 1, 1] + (t+1)*ks6[4, 3, 1] + t*ks6[4, 4]
        + ks6[5, 1, 1, 1] + ks6[5, 2, 1] + t*ks6[5, 3] + ks6[6, 1, 1]

    The `k`-split basis is a second basis of the ring spanned by the `k`-Schur
    functions with a parameter `t`.  The `k`-split basis has the property that
    `Q'_\lambda[X;t]` expands positively in the `k`-split basis and the
    `k`-split basis conjecturally expands positively in the `k`-Schur functions.
    The definition can be found in [LLMSSZ]_ p. 81.::

        sage: ksp3 = SymS3.ksplit()
        sage: ksp3(Qp[2,1,1,1])
        ksp3[2, 1, 1, 1] + t^2*ksp3[2, 2, 1] + (t^3+t^2)*ksp3[3, 1, 1] + t^4*ksp3[3, 2]
        sage: [ks(ksp3(la)) for la in sorted(ksp3(Qp[2,1,1,1]).support())]
        [ks3[2, 1, 1, 1] + t*ks3[2, 2, 1], ks3[2, 2, 1], ks3[3, 1, 1], ks3[3, 2]]

    .. rubric:: dual `k`-Schur functions

    The dual space to the subspace spanned by the `k`-Schur functions is most naturally
    realized as a quotient of the ring of symmetric functions by an ideal.  When `t=1`
    the ideal is generated by the monomial symmetric functions indexed by partitions
    whose first part is greater than `k`.::

        sage: Sym = SymmetricFunctions(QQ)
        sage: SymQ3 = Sym.kBoundedQuotient(3,t=1)
        sage: km = SymQ3.kmonomial()
        sage: km[2,1]*km[2,1]
        4*m3[2, 2, 1, 1] + 6*m3[2, 2, 2] + 2*m3[3, 2, 1] + 2*m3[3, 3]
        sage: F = SymQ3.affineSchur()
        sage: F[2,1]*F[2,1]
        2*F3[1, 1, 1, 1, 1, 1] + 4*F3[2, 1, 1, 1, 1] + 4*F3[2, 2, 1, 1] + 4*F3[2, 2, 2]
        + 2*F3[3, 1, 1, 1] + 4*F3[3, 2, 1] + 2*F3[3, 3]

    When `t` is not equal to `1`, the subspace spanned by the `k`-Schur functions is
    realized as a quotient of the ring of symmetric functions by the ideal generated by
    the Hall-Littlewood symmetric functions in the `P` basis indexed by partitions with
    first part greater than `k`.::

        sage: Sym = SymmetricFunctions(FractionField(QQ['t']))
        sage: SymQ3 = Sym.kBoundedQuotient(3)
        sage: kHLP = SymQ3.kHallLittlewoodP()
        sage: kHLP[2,1]*kHLP[2,1]
        (t^2+2*t+1)*HLP3[2, 2, 1, 1] + (t^3+2*t^2+2*t+1)*HLP3[2, 2, 2]
        + (-t^4-t^3+t+1)*HLP3[3, 1, 1, 1] + (-t^2+t+2)*HLP3[3, 2, 1] + (t+1)*HLP3[3, 3]
        sage: HLP = Sym.hall_littlewood().P()
        sage: kHLP(HLP[3,1])
        HLP3[3, 1]
        sage: kHLP(HLP[4])
        0

    In this space, the basis which is dual to the `k`-Schur functions conjecturally
    expands positively in the `k`-bounded Hall-Littlewood functions and has positive
    structure coefficients.::

        sage: dks = SymQ3.dual_k_Schur()
        sage: kHLP(dks[2,2])
        (t^4+t^2)*HLP3[1, 1, 1, 1] + t*HLP3[2, 1, 1] + HLP3[2, 2]
        sage: dks[2,1]*dks[1,1]
        (t^2+t)*dks3[1, 1, 1, 1, 1] + (t+1)*dks3[2, 1, 1, 1] + (t+1)*dks3[2, 2, 1]
        + dks3[3, 1, 1] + dks3[3, 2]

    At `t=1` the `k`-bounded Hall-Littlewood basis is equal to the `k`-bounded monomial
    basis and the dual `k`-Schur elements are equal to the affine Schur basis.  The
    `k`-bounded monomial basis and affine Schur functions are faster and should be used
    instead of the `k`-bounded Hall-Littlewood P basis and dual `k`-Schur functions when
    `t=1`. ::

        sage: SymQ3 = Sym.kBoundedQuotient(3,t=1)
        sage: dks = SymQ3.dual_k_Schur()
        sage: F = SymQ3.affineSchur()
        sage: F[3,1]==dks[3,1]
        True


    .. _`Representation theory of the symmetric group`:

    .. rubric:: Representation theory of the symmetric group

    The Schur functions `s_\lambda` can be interpreted as irreducible characters of the symmetric
    group `S_n`, where `n` is the size of the partition `\lambda`. Since the Schur functions of
    degree `n` form a basis of the symmetric functions of degree `n`, it
    follows that an arbitrary symmetric function (homogeneous of degree
    `n`) may be interpreted as a function on the symmetric group. In this
    interpretation, the power sum symmetric function `p_\lambda` is the characteristic
    function of the conjugacy class with shape `\lambda`, multiplied by the order of
    the centralizer of an element. Hence the irreducible characters can be computed
    as follows::

        sage: Sym = SymmetricFunctions(QQ)
        sage: s = Sym.schur()
        sage: p = Sym.power()
        sage: P = Partitions(5).list()
        sage: P = [P[i] for i in range(len(P)-1,-1,-1)]
        sage: M = matrix([[s[P[i]].scalar(p[P[j]]) for j in range(len(P))] for i in range(len(P))])
        sage: M
        [ 1 -1  1  1 -1 -1  1]
        [ 4 -2  0  1  1  0 -1]
        [ 5 -1  1 -1 -1  1  0]
        [ 6  0 -2  0  0  0  1]
        [ 5  1  1 -1  1 -1  0]
        [ 4  2  0  1 -1  0 -1]
        [ 1  1  1  1  1  1  1]

    We can indeed check that this agrees with the character table of `S_5`::

        sage: SymmetricGroup(5).character_table() == M
        True


    .. rubric:: Inner plethysm

    The operation of inner plethysm (``f.inner_plethysm(g)``) models the
    composition of the `S_n`-representation represented by `g` with the
    `GL_m`-representation whose character is `f`.  See the documentation of
    ``inner_plethysm``, [ST94]_ or [EnumComb2]_, exercise 7.74 solutions for more
    information::

        sage: s = SymmetricFunctions(QQ).schur()
        sage: f = s[2]^2
        sage: f.inner_plethysm(s[2])
        s[2]


    .. rubric:: More specific applications

    The first part of this tutorial was meant to present general use of symmetric functions in Sage.
    Here are now more specific applications.

    Let us explore the other operations of ``p`` (the power-sum
    realization of :math:`\operatorname{Sym}`). We can ask for
    the mathematical properties of ``p``.

    ::

        sage: p.categories() #not tested
        [Category of graded bases of Symmetric Functions over Fraction Field of Multivariate Polynomial Ring in q, t over Rational Field,
         Category of filtered bases of Symmetric Functions over Fraction Field of Multivariate Polynomial Ring in q, t over Rational Field,
         Category of bases of Symmetric Functions over Fraction Field of Multivariate Polynomial Ring in q, t over Rational Field,
         Category of graded hopf algebras with basis over Fraction Field of Multivariate Polynomial Ring in q, t over Rational Field,
         Category of filtered hopf algebras with basis over Fraction Field of Multivariate Polynomial Ring in q, t over Rational Field,
         Category of hopf algebras with basis over Fraction Field of Multivariate Polynomial Ring in q, t over Rational Field,
         Category of realizations of hopf algebras over Fraction Field of Multivariate Polynomial Ring in q, t over Rational Field,
         Category of hopf algebras over Fraction Field of Multivariate Polynomial Ring in q, t over Rational Field,
         Category of graded algebras with basis over Fraction Field of Multivariate Polynomial Ring in q, t over Rational Field,
         Category of filtered algebras with basis over Fraction Field of Multivariate Polynomial Ring in q, t over Rational Field,
         Category of bialgebras with basis over Fraction Field of Multivariate Polynomial Ring in q, t over Rational Field,
         Category of algebras with basis over Fraction Field of Multivariate Polynomial Ring in q, t over Rational Field,
         Category of graded algebras over Fraction Field of Multivariate Polynomial Ring in q, t over Rational Field,
         Category of commutative algebras over Fraction Field of Multivariate Polynomial Ring in q, t over Rational Field,
         Category of filtered algebras over Fraction Field of Multivariate Polynomial Ring in q, t over Rational Field,
         Category of bialgebras over Fraction Field of Multivariate Polynomial Ring in q, t over Rational Field,
         Category of algebras over Fraction Field of Multivariate Polynomial Ring in q, t over Rational Field,
         Category of commutative rings,
         Category of rings,
         Category of associative algebras over Fraction Field of Multivariate Polynomial Ring in q, t over Rational Field,
         Category of rngs,
         Category of semirings,
         Category of associative additive commutative additive associative additive unital distributive magmas and additive magmas,
         Category of unital algebras with basis over Fraction Field of Multivariate Polynomial Ring in q, t over Rational Field,
         Category of magmatic algebras with basis over Fraction Field of Multivariate Polynomial Ring in q, t over Rational Field,
         Category of unital algebras over Fraction Field of Multivariate Polynomial Ring in q, t over Rational Field,
         Category of magmatic algebras over Fraction Field of Multivariate Polynomial Ring in q, t over Rational Field,
         Category of additive commutative additive associative additive unital distributive magmas and additive magmas,
         Category of additive commutative additive associative distributive magmas and additive magmas,
         Category of additive associative distributive magmas and additive magmas,
         Category of distributive magmas and additive magmas,
         Category of magmas and additive magmas,
         Category of commutative monoids,
         Category of monoids,
         Category of semigroups,
         Category of realizations of unital magmas,
         Category of realizations of magmas,
         Category of commutative magmas,
         Category of unital magmas,
         Category of magmas,
         Category of graded modules with basis over Fraction Field of Multivariate Polynomial Ring in q, t over Rational Field,
         Category of filtered modules with basis over Fraction Field of Multivariate Polynomial Ring in q, t over Rational Field,
         Category of coalgebras with basis over Fraction Field of Multivariate Polynomial Ring in q, t over Rational Field,
         Category of vector spaces with basis over Fraction Field of Multivariate Polynomial Ring in q, t over Rational Field,
         Category of modules with basis over Fraction Field of Multivariate Polynomial Ring in q, t over Rational Field,
         Category of graded modules over Fraction Field of Multivariate Polynomial Ring in q, t over Rational Field,
         Category of realizations of coalgebras over Fraction Field of Multivariate Polynomial Ring in q, t over Rational Field,
         Category of filtered modules over Fraction Field of Multivariate Polynomial Ring in q, t over Rational Field,
         Category of coalgebras over Fraction Field of Multivariate Polynomial Ring in q, t over Rational Field,
         Category of vector spaces over Fraction Field of Multivariate Polynomial Ring in q, t over Rational Field,
         Category of modules over Fraction Field of Multivariate Polynomial Ring in q, t over Rational Field,
         Category of bimodules over Fraction Field of Multivariate Polynomial Ring in q, t over Rational Field on the left and Fraction Field of Multivariate Polynomial Ring in q, t over Rational Field on the right,
         Category of right modules over Fraction Field of Multivariate Polynomial Ring in q, t over Rational Field,
         Category of left modules over Fraction Field of Multivariate Polynomial Ring in q, t over Rational Field,
         Category of commutative additive groups,
         Category of additive groups,
         Category of additive inverse additive unital additive magmas,
         Category of commutative additive monoids,
         Category of additive monoids,
         Category of additive unital additive magmas,
         Category of commutative additive semigroups,
         Category of additive commutative additive magmas,
         Category of additive semigroups,
         Category of additive magmas,
         Category of realizations of Symmetric Functions over Fraction Field of Multivariate Polynomial Ring in q, t over Rational Field,
         Category of realizations of sets,
         Category of sets,
         Category of sets with partial maps,
         Category of objects]


    .. rubric:: Hopf structure and important identities

    Many important identities between symmetric functions can be linked to "the"
    Hopf algebra structure on the ring of symmetric functions.
    In part, this means that we have a **coproduct** on the
    ring of symmetric functions. It is an algebra homomorphism
    :math:`\Delta : \operatorname{Sym} \to \operatorname{Sym} \otimes \operatorname{Sym}`
    that can be defined by

    .. MATH::

        \Delta(g) = \sum_{\mu, \nu} a_{\mu,\nu}\, s_\mu\otimes s_\nu

    where

    .. MATH::

        g(\mathbf{x},\mathbf{y})
        = \sum_{\mu, \nu} a_{\mu,\nu}\, s_\mu(\mathbf{x}) s_\nu(\mathbf{y}),

    where :math:`\mathbf{x} = (x_1, x_2, x_3, \ldots)` and
    :math:`\mathbf{y} = (y_1, y_2, y_3, \ldots)` are two
    infinite sequences of indeterminates and where
    :math:`g(\mathbf{x},\mathbf{y})` stands for
    :math:`g(x_1, x_2, x_3, \ldots, y_1, y_2, y_3, \ldots)`.

    For instance, we have ::

        sage: Sym = SymmetricFunctions(QQ)
        sage: Sym.inject_shorthands(verbose=false)
        sage: One=s[0]
        sage: X=s[1]
        sage: Y=tensor([X,One])
        sage: Z=tensor([One,X])

    ::

        sage: s[3].coproduct()
        s[] # s[3] + s[1] # s[2] + s[2] # s[1] + s[3] # s[]
        sage: s[3](Y+Z)
        s[] # s[3] + s[1] # s[2] + s[2] # s[1] + s[3] # s[]
        sage: s[3,2,1].coproduct()
        s[] # s[3, 2, 1] + s[1] # s[2, 2, 1] + s[1] # s[3, 1, 1] + s[1] # s[3, 2] + s[1, 1] # s[2, 1, 1] + s[1, 1] # s[2, 2] + s[1, 1] # s[3, 1] + s[1, 1, 1] # s[2, 1] + s[2] # s[2, 1, 1] + s[2] # s[2, 2] + s[2] # s[3, 1] + s[2, 1] # s[1, 1, 1] + 2*s[2, 1] # s[2, 1] + s[2, 1] # s[3] + s[2, 1, 1] # s[1, 1] + s[2, 1, 1] # s[2] + s[2, 2] # s[1, 1] + s[2, 2] # s[2] + s[2, 2, 1] # s[1] + s[3] # s[2, 1] + s[3, 1] # s[1, 1] + s[3, 1] # s[2] + s[3, 1, 1] # s[1] + s[3, 2] # s[1] + s[3, 2, 1] # s[]
        sage: s[3,2,1](Y+Z)
        s[] # s[3, 2, 1] + s[1] # s[2, 2, 1] + s[1] # s[3, 1, 1] + s[1] # s[3, 2] + s[1, 1] # s[2, 1, 1] + s[1, 1] # s[2, 2] + s[1, 1] # s[3, 1] + s[1, 1, 1] # s[2, 1] + s[2] # s[2, 1, 1] + s[2] # s[2, 2] + s[2] # s[3, 1] + s[2, 1] # s[1, 1, 1] + 2*s[2, 1] # s[2, 1] + s[2, 1] # s[3] + s[2, 1, 1] # s[1, 1] + s[2, 1, 1] # s[2] + s[2, 2] # s[1, 1] + s[2, 2] # s[2] + s[2, 2, 1] # s[1] + s[3] # s[2, 1] + s[3, 1] # s[1, 1] + s[3, 1] # s[2] + s[3, 1, 1] # s[1] + s[3, 2] # s[1] + s[3, 2, 1] # s[]

    In particular, we get

    .. MATH:: \Delta(h_n) = \sum_{k+j=n} h_k\otimes h_j

    for every :math:`n \geq 0`::

        sage: h[4].coproduct()
        h[] # h[4] + h[1] # h[3] + h[2] # h[2] + h[3] # h[1] + h[4] # h[]
        sage: h[4](Y+Z)
        h[] # h[4] + h[1] # h[3] + h[2] # h[2] + h[3] # h[1] + h[4] # h[]
        sage: tensor([h,e])(h[4](Y-Z))
        h[] # e[4] - h[1] # e[3] + h[2] # e[2] - h[3] # e[1] + h[4] # e[]
        sage: s[3,1](Y-Z)
        s[] # s[2, 1, 1] - s[1] # s[1, 1, 1] - s[1] # s[2, 1] + s[1, 1] # s[1, 1] + s[2] # s[1, 1] + s[2] # s[2] - s[2, 1] # s[1] - s[3] # s[1] + s[3, 1] # s[]

    The coproduct is an algebra morphism, and is therefore
    determined by its values on the generators; the power sum generators
    are primitive::

        sage: p[1].coproduct()
        p[] # p[1] + p[1] # p[]
        sage: p[2].coproduct()
        p[] # p[2] + p[2] # p[]

    The coproduct, being cocommutative on the generators,
    is cocommutative everywhere::

        sage: p[2, 1].coproduct()
        p[] # p[2, 1] + p[1] # p[2] + p[2] # p[1] + p[2, 1] # p[]

    This coproduct, along with the counit which maps every symmetric function
    to its 0-th homogeneous component, makes the ring of symmetric functions
    into a graded connected bialgebra. It is known that every graded connected
    bialgebra has an antipode. For the ring of symmetric functions, the antipode
    can be characterized explicitly: The antipode is an anti-algebra morphism
    (thus an algebra morphism, since our algebra is commutative) which sends
    :math:`p_{\lambda}` to :math:`(-1)^{\operatorname{length}(\lambda)} p_{\lambda}` for every
    partition :math:`\lambda`. Thus, in particular, it sends the generators on the
    :math:`p` basis to their opposites::

        sage: p[3].antipode()
        -p[3]
        sage: p[3](-X)
        -p[3]
        sage: s[3,1,1,1,1].antipode()
        -s[5, 1, 1]
        sage: s[3,1,1,1,1](-X)
        -s[5, 1, 1]

    If the base ring is a :math:`\QQ`-algebra, then
    the graded connected bialgebra :math:`\operatorname{Sym}`
    has a rather simple structure: It is (isomorphic to) the
    symmetric algebra of its space of primitives (which is spanned by the
    power-sum symmetric functions).

    Here are further examples::

        sage: g = s[2]^2
        sage: g.antipode()
        s[1, 1, 1, 1] + s[2, 1, 1] + s[2, 2]
        sage: g.coproduct()
        s[] # s[2, 2] + s[] # s[3, 1] + s[] # s[4] + 2*s[1] # s[2, 1]
         + 2*s[1] # s[3] + s[1, 1] # s[1, 1] + s[1, 1] # s[2]
         + s[2] # s[1, 1] + 3*s[2] # s[2] + 2*s[2, 1] # s[1]
         + s[2, 2] # s[] + 2*s[3] # s[1] + s[3, 1] # s[] + s[4] # s[]
        sage: g.coproduct().apply_multilinear_morphism( lambda x,y: x*y.antipode() )
        0

    In the interpretation of symmetric functions as characters on the symmetric group,
    the multiplication and comultiplication are interpreted as induction
    (from :math:`S_n\times S_m` to :math:`S_{n+m}`) and restriction, respectively.
    The Schur functions can also be interpreted as characters of :math:`GL_n`, 
    see `Partitions and Schur functions`__.

    .. rubric:: Skew Schur functions and skewing

    Skew Schur functions arise when one considers the effect
    of the coproduct on Schur functions themselves:

    .. MATH:: \Delta(s_\lambda) = \sum_{\mu\subseteq \lambda} s_{\lambda/\mu}\otimes s_\mu.

    Skew Schur functions are also implemented in Sage.
    For instance, we have the skew Schur function
    :math:`s_{321/2}`::

        sage: s([[3,2,1],[2]])
        s[2, 1, 1] + s[2, 2] + s[3, 1]
        sage: s([[2,1],[1]])
        s[1, 1] + s[2]

    Given a symmetric function :math:`f`, we can define a linear
    map :math:`f^\perp : \operatorname{Sym} \to \operatorname{Sym}` by
    requiring :math:`f^\perp(g)` to depend linearly on each of
    :math:`f` and :math:`g` and to satisfy

    .. MATH:: `s_\mu^\perp s_\lambda = s_{\lambda/\mu}`

    for all partitions :math:`\lambda` and :math:`\mu`, where
    :math:`s_{\lambda/\mu}` is defined to be `0` if
    :math:`\mu \not\subseteq \lambda`. This linear map
    :math:`f^\perp` is called *skewing by* :math:`f`, and is
    implemented as ``skew_by(f)`` in Sage.  For example::

        sage: s[3,2,1].skew_by(s[2]) # same as s([[3,2,1],[2]])
        s[2, 1, 1] + s[2, 2] + s[3, 1]
        sage: s[3,2,1].skew_by(s[4])
        0

    ::

        sage: add(tensor([s[3,2,1].skew_by(s(mu)),s(mu)]) for k in range(7) for mu in Partitions(k))
        s[] # s[3, 2, 1] + s[1] # s[2, 2, 1] + s[1] # s[3, 1, 1] + s[1] # s[3, 2] + s[1, 1] # s[2, 1, 1] + s[1, 1] # s[2, 2] + s[1, 1] # s[3, 1] + s[1, 1, 1] # s[2, 1] + s[2] # s[2, 1, 1] + s[2] # s[2, 2] + s[2] # s[3, 1] + s[2, 1] # s[1, 1, 1] + 2*s[2, 1] # s[2, 1] + s[2, 1] # s[3] + s[2, 1, 1] # s[1, 1] + s[2, 1, 1] # s[2] + s[2, 2] # s[1, 1] + s[2, 2] # s[2] + s[2, 2, 1] # s[1] + s[3] # s[2, 1] + s[3, 1] # s[1, 1] + s[3, 1] # s[2] + s[3, 1, 1] # s[1] + s[3, 2] # s[1] + s[3, 2, 1] # s[]

    .. rubric:: Cauchy kernel formula

    Consider the family

    .. MATH:: \mathbf{x}\mathbf{y} := (x_i y_j)_{i \geq 1,\ j \geq 1}

    of indeterminates. The Cauchy kernel is the expression

    .. MATH:: \sum_{n\geq 0} h_n(\mathbf{x}\mathbf{y})
        = \prod_{i,j}\frac{1}{1-x_iy_j}

    written here using plethystic notation.
    Its degree-:math:`n` homogeneous component plays a crucial role in the description 
    of "dual bases" with respect to the scalar product.
    We can embed the tensor product
    :math:`\operatorname{Sym} \otimes \operatorname{Sym}` into
    the ring of formal power series in the indeterminates
    :math:`x_1, x_2, x_3, \ldots, y_1, y_2, y_3, \ldots` by
    identifying each pure tensor :math:`f \otimes g` with
    the power series :math:`f(x_1, x_2, x_3, \ldots, y_1, y_2, y_3, \ldots)`;
    then, for a given :math:`n \geq 0`, we have

    .. MATH:: h_n(\mathbf{x}\mathbf{y})
        = \sum_{\mu\vdash n} F_\mu\otimes G_\mu
        \qquad \text{iff} \qquad
        \left\langle F_\mu,G_\lambda\right\rangle
        = \delta_{\mu, \lambda}
        \qquad \text{for all partitions } \mu \text{ and }
        \lambda \text{ of } n ,

    where :math:`\delta_{\mu, \lambda}` is the Kronecker
    delta. One says that :math:`\{F_\mu\}_\mu` and
    :math:`\{G_\lambda\}_\lambda` are **dual bases**.
    Schur functions are self-dual; the dual of the
    :math:`h_{\mu}` are the :math:`m_\mu`; that of the
    :math:`p_\mu` are the :math:`p_{\mu}/z_{\mu}`.
    The "forgotten" symmetric functions :math:`f_{\mu}`
    appear as the dual of the :math:`e_{\mu}`.

    ::

        sage: h4xy=add(tensor([s(mu),s(mu)]) for mu in Partitions(4)); h4xy
        s[1, 1, 1, 1] # s[1, 1, 1, 1] + s[2, 1, 1] # s[2, 1, 1] + s[2, 2] # s[2, 2] + s[3, 1] # s[3, 1] + s[4] # s[4]
        sage: s[4](Y*Z)
        s[1, 1, 1, 1] # s[1, 1, 1, 1] + s[2, 1, 1] # s[2, 1, 1] + s[2, 2] # s[2, 2] + s[3, 1] # s[3, 1] + s[4] # s[4]
        sage: tensor([h,m])(h4xy)
        h[1, 1, 1, 1] # m[1, 1, 1, 1] + h[2, 1, 1] # m[2, 1, 1] + h[2, 2] # m[2, 2] + h[3, 1] # m[3, 1] + h[4] # m[4]
        sage: tensor([e,h])(h4xy)
        e[1, 1, 1, 1] # h[4] + e[2, 1, 1] # h[3, 1] - 4*e[2, 1, 1] # h[4] + e[2, 2] # h[2, 2] - 2*e[2, 2] # h[3, 1] + 2*e[2, 2] # h[4] + e[3, 1] # h[2, 1, 1] - 2*e[3, 1] # h[2, 2] - e[3, 1] # h[3, 1] + 4*e[3, 1] # h[4] + e[4] # h[1, 1, 1, 1] - 4*e[4] # h[2, 1, 1] + 2*e[4] # h[2, 2] + 4*e[4] # h[3, 1] - 4*e[4] # h[4]
        sage: tensor([p,p])(h4xy)
        1/24*p[1, 1, 1, 1] # p[1, 1, 1, 1] + 1/4*p[2, 1, 1] # p[2, 1, 1] + 1/8*p[2, 2] # p[2, 2] + 1/3*p[3, 1] # p[3, 1] + 1/4*p[4] # p[4]


    __ ../../../../../thematic_tutorials/lie/lie_basics.html#partitions-and-schur-polynomials

    .. rubric:: The Kronecker product

    As in the section about the **Representation theory of the symmetric group**,
    a symmetric function may be considered as a class function on the symmetric
    group where the elements :math:`p_\mu/z_\mu` are the indicators of a permutation
    having cycle structure :math:`\mu`.  The Kronecker product of two symmetric
    functions corresponds to the pointwise product of these class functions.

    Since the Schur functions are the irreducible characters
    of the symmetric group under this identification, the Kronecker
    product of two Schur functions corresponds to the internal
    tensor product of two irreducible symmetric group representations.

    Under this identification, the Kronecker
    product of :math:`p_\mu/z_\mu` and :math:`p_\nu/z_\nu` is :math:`p_\mu/z_\mu`
    if :math:`\mu=\nu`, and is :math:`0` otherwise.

    ``internal_product``, ``kronecker_product``, ``inner_tensor`` and
    ``itensor`` are different names for the same function.

    ::

        sage: g
        s[2, 2] + s[3, 1] + s[4]
        sage: g.kronecker_product(g)
        s[1, 1, 1, 1] + 3*s[2, 1, 1] + 4*s[2, 2] + 5*s[3, 1] + 3*s[4]
        sage: g.kronecker_product(s[4])
        s[2, 2] + s[3, 1] + s[4]
        sage: g.kronecker_product(e[4])
        s[1, 1, 1, 1] + s[2, 1, 1] + s[2, 2]
        sage: g.omega()
        s[1, 1, 1, 1] + s[2, 1, 1] + s[2, 2]
        sage: Matrix([[p(mu).kronecker_product(p(nu)/zee(nu)) for nu in Partitions(5)] for mu in Partitions(5)])
        [            p[5]                0                0                0                0                0                0]
        [               0          p[4, 1]                0                0                0                0                0]
        [               0                0          p[3, 2]                0                0                0                0]
        [               0                0                0       p[3, 1, 1]                0                0                0]
        [               0                0                0                0       p[2, 2, 1]                0                0]
        [               0                0                0                0                0    p[2, 1, 1, 1]                0]
        [               0                0                0                0                0                0 p[1, 1, 1, 1, 1]]



    .. rubric:: Implementing new bases

    In order to implement a new symmetric function basis, Sage will need
    to know at a minimum how to change back and forth between it and
    at least one other basis (although they do not necessarily have to
    be the same basis).
    All of the standard functions associated with the basis will then
    automatically have a default implementation by way of these
    conversions (although a more specific implementation may be more
    efficient).

    To present an idea of how this is done, we will create
    here the example of how to implement the basis :math:`s_\mu[X(1-t)]`.

    To begin, we import the class
    :class:`sage.combinat.sf.sfa.SymmetricFunctionAlgebra_generic()`.  Our
    new basis will inherit all of the default methods from this class::

        sage: from sage.combinat.sf.sfa import SymmetricFunctionAlgebra_generic as SFA_generic

    Now the basis we are creating has a parameter :math:`t` which one can
    specialize. In this example we will convert to and from the Schur
    basis.  For this we implement methods ``_self_to_s`` and ``_s_to_self``.
    By registering these two functions as coercions, Sage then knows
    automatically how to change between any two bases for
    which there is a path of changes of basis.

    ::

        sage: from sage.categories.morphism import SetMorphism
        sage: class SFA_st(SFA_generic):
        ....:     def __init__(self, Sym, t):
        ....:         SFA_generic.__init__(self, Sym, basis_name=
        ....:           "Schur functions with a plethystic substitution of X -> X(1-t)",
        ....:           prefix='st')
        ....:         self._s = Sym.s()
        ....:         self.t = Sym.base_ring()(t)
        ....:         cat = HopfAlgebras(Sym.base_ring()).WithBasis()
        ....:         self.register_coercion(
        ....:           SetMorphism(Hom(self._s, self, cat), self._s_to_self))
        ....:         self._s.register_coercion(
        ....:           SetMorphism(Hom(self, self._s, cat), self._self_to_s))
        ....:     def _s_to_self(self, f):
        ....:         # f is a Schur function and the output is in the st basis
        ....:         return self._from_dict(f.theta_qt(0,self.t)._monomial_coefficients)
        ....:     def _self_to_s(self, f):
        ....:         # f is in the st basis and the output is in the Schur basis
        ....:         return self._s.sum(cmu*self._s(mu).theta_qt(self.t,0) for mu,cmu in f)
        ....:     class Element(SFA_generic.Element):
        ....:         pass

    An instance of this basis is created by calling it with a symmetric
    function ring ``Sym`` and a parameter ``t`` which is in the base ring
    of ``Sym``.  The ``Element`` class inherits all of the methods from
    :class:`sage.combinat.sf.sfa.SymmetricFunctionAlgebra_generic_Element`.

    In Macdonald's work ([Mac1995]_ page 354), this basis is denoted by
    :math:`S_\lambda(x;t)` and the change of basis coefficients of the
    Macdonald ``J`` basis are the coefficients :math:`K_{\lambda\mu}(q,t)`.
    Here is an example of its use::

        sage: QQqt = QQ['q','t'].fraction_field()
        sage: (q,t) = QQqt.gens()
        sage: st = SFA_st(SymmetricFunctions(QQqt),t)
        sage: st
        Symmetric Functions over Fraction Field of Multivariate Polynomial
         Ring in q, t over Rational Field in the Schur functions with a
         plethystic substitution of X -> X(1-t) basis
        sage: st[2,1] * st[1]
        st[2, 1, 1] + st[2, 2] + st[3, 1]
        sage: st([2]).coproduct()
         st[] # st[2] + st[1] # st[1] + st[2] # st[]
        sage: J = st.symmetric_function_ring().macdonald().J()
        sage: st(J[2,1])
        q*st[1, 1, 1] + (q*t+1)*st[2, 1] + t*st[3]


    .. rubric:: Acknowledgements

    The design is heavily inspired from the implementation of
    symmetric functions in MuPAD-Combinat (see [HT04]_ and [FD06]_).


    .. rubric:: Further tests

    TESTS::

        sage: Sym = SymmetricFunctions(QQ)
        sage: Sym
        Symmetric Functions over Rational Field
        sage: h = Sym.h(); e = Sym.e(); s = Sym.s(); m = Sym.m(); p = Sym.p()
        sage: ( ( h[2,1] * ( 1 + 3 * h[2,1]) ) + s[2]. antipode()) . coproduct()
        h[] # h[1, 1] - h[] # h[2] + h[] # h[2, 1] + 3*h[] # h[2, 2, 1, 1] + h[1] # h[1] + h[1] # h[1, 1]
        + h[1] # h[2] + 6*h[1] # h[2, 1, 1, 1] + 6*h[1] # h[2, 2, 1] + h[1, 1] # h[] + h[1, 1] # h[1]
        + 3*h[1, 1] # h[1, 1, 1, 1] + 12*h[1, 1] # h[2, 1, 1] + 3*h[1, 1] # h[2, 2] + 6*h[1, 1, 1] # h[1, 1, 1]
        + 6*h[1, 1, 1] # h[2, 1] + 3*h[1, 1, 1, 1] # h[1, 1] - h[2] # h[] + h[2] # h[1] + 6*h[2] # h[2, 1, 1]
        + h[2, 1] # h[] + 6*h[2, 1] # h[1, 1, 1] + 12*h[2, 1] # h[2, 1] + 12*h[2, 1, 1] # h[1, 1]
        + 6*h[2, 1, 1] # h[2] + 6*h[2, 1, 1, 1] # h[1] + 3*h[2, 2] # h[1, 1] + 6*h[2, 2, 1] # h[1] + 3*h[2, 2, 1, 1] # h[]

    .. TODO::

        - Introduce fields with degree 1 elements as in
          MuPAD-Combinat, to get proper plethysm.
        - Use UniqueRepresentation to get rid of all the manual cache
          handling for the bases
        - Devise a mechanism so that pickling bases of symmetric
          functions pickles the coercions which have a cache.
    """

    def __init__(self, R):
        r"""
        Initialization of ``self``.

        INPUT:

        - ``R`` -- a ring

        EXAMPLES::

            sage: Sym = SymmetricFunctions(QQ)

        TESTS::

            sage: Sym1 = SymmetricFunctions(FiniteField(23))
            sage: Sym2 = SymmetricFunctions(Integers(23))
            sage: TestSuite(Sym).run()

        """
        # change the line below to assert(R in Rings()) once MRO issues from #15536, #15475 are resolved
        assert(R in Fields() or R in Rings()) # side effect of this statement assures MRO exists for R
        self._base = R # Won't be needed when CategoryObject won't override anymore base_ring
        Parent.__init__(self, category = GradedHopfAlgebras(R).WithRealizations())

    def a_realization(self):
        r"""
        Return a particular realization of ``self`` (the Schur basis).

        EXAMPLES::

            sage: Sym = SymmetricFunctions(QQ)
            sage: Sym.a_realization()
            Symmetric Functions over Rational Field in the Schur basis
        """
        return self.schur()

    def _repr_(self): # could be taken care of by the category
        r"""
        Representation of ``self``

        TESTS::

            sage: SymmetricFunctions(RR) # indirect doctest
            Symmetric Functions over Real Field with 53 bits of precision
        """
        return "Symmetric Functions over %s"%self.base_ring()

    def schur(self):
        r"""
        The Schur basis of the Symmetric Functions

        EXAMPLES::

            sage: SymmetricFunctions(QQ).schur()
            Symmetric Functions over Rational Field in the Schur basis
        """
        return schur.SymmetricFunctionAlgebra_schur(self)
    s = schur
    Schur = schur # Currently needed by SymmetricFunctions.__init_extra__
                  # and sfa.GradedSymmetricFunctionsBases.corresponding_basis_over

    def powersum(self):
        r"""
        The power sum basis of the Symmetric Functions

        EXAMPLES::

            sage: SymmetricFunctions(QQ).powersum()
            Symmetric Functions over Rational Field in the powersum basis
        """
        return powersum.SymmetricFunctionAlgebra_power(self)
    p = powersum
    power = powersum # Todo: get rid of this one when it won't be needed anymore

    def complete(self):
        r"""
        The complete basis of the Symmetric Functions

        EXAMPLES::

            sage: SymmetricFunctions(QQ).complete()
            Symmetric Functions over Rational Field in the homogeneous basis
        """
        return homogeneous.SymmetricFunctionAlgebra_homogeneous(self)
    h = complete
    homogeneous = complete

    def elementary(self):
        r"""
        The elementary basis of the Symmetric Functions

        EXAMPLES::

            sage: SymmetricFunctions(QQ).elementary()
            Symmetric Functions over Rational Field in the elementary basis
        """
        return elementary.SymmetricFunctionAlgebra_elementary(self)
    e = elementary

    def monomial(self):
        r"""
        The monomial basis of the Symmetric Functions

        EXAMPLES::

            sage: SymmetricFunctions(QQ).monomial()
            Symmetric Functions over Rational Field in the monomial basis
        """
        return monomial.SymmetricFunctionAlgebra_monomial(self)
    m = monomial

    def witt(self, coerce_h=True, coerce_e=False, coerce_p=False):
        r"""
        The Witt basis of the symmetric functions.

        EXAMPLES::

            sage: SymmetricFunctions(QQ).witt()
            Symmetric Functions over Rational Field in the Witt basis
            sage: SymmetricFunctions(QQ).witt(coerce_p=True)
            Symmetric Functions over Rational Field in the Witt basis
            sage: SymmetricFunctions(QQ).witt(coerce_h=False, coerce_e=True, coerce_p=True)
            Symmetric Functions over Rational Field in the Witt basis
        """
        from . import witt
        return witt.SymmetricFunctionAlgebra_witt(self, coerce_h=coerce_h, coerce_e=coerce_e, coerce_p=coerce_p)
    w = witt
    # Currently needed by sfa.GradedSymmetricFunctionsBases.corresponding_basis_over
    Witt = witt

    def irreducible_symmetric_group_character(self):
        r"""
        The irreducible `S_n` character basis of the Symmetric Functions.

        This basis has the property that if the element indexed by the
        partition `\lambda` is evaluated at the roots of a permutation of
        cycle structure `\rho` then the value is the irreducible character
        `\chi^{(|\rho|-|\lambda|,\lambda)}(\rho)`.

        In terms of methods that are implemented in Sage, if ``n`` is
        a sufficiently large integer, then
        ``st(lam).character_to_frobenius_image(n)`` is equal the Schur function
        indexed by ``[n-sum(lam)]+lam``.

        This basis is introduced in [OZ2015]_.

        .. SEEALSO::

            :meth:`~sage.combinat.sf.sfa.SymmetricFunctionAlgebra_generic_Element.character_to_frobenius_image`,
            :meth:`~sage.combinat.sf.sfa.SymmetricFunctionAlgebra_generic_Element.eval_at_permutation_roots`

        EXAMPLES::

            sage: SymmetricFunctions(QQ).irreducible_symmetric_group_character()
            Symmetric Functions over Rational Field in the irreducible symmetric group character basis
            sage: st = SymmetricFunctions(QQ).st()
            sage: s = SymmetricFunctions(QQ).s()
            sage: s(st([3,2]).character_to_frobenius_image(9))
            s[4, 3, 2]
            sage: s(st([3,2]).character_to_frobenius_image(7))
            0
            sage: s(st([3,2]).character_to_frobenius_image(6))
            -s[2, 2, 2]
            sage: list(SymmetricGroup(5).character_table()[-2])
            [4, 2, 0, 1, -1, 0, -1]
            sage: list(reversed([st([1]).eval_at_permutation_roots(rho)
            ....:   for rho in Partitions(5)]))
            [4, 2, 0, 1, -1, 0, -1]
        """
        from .character import irreducible_character_basis
        return irreducible_character_basis(self, 'st')

    st = irreducible_symmetric_group_character

    def induced_trivial_character(self):
        r"""
        The induced trivial character basis of the Symmetric Functions.

        The trivial character of

        .. MATH::

            S_{n-|\lambda|} \times S_{\lambda_1} \times S_{\lambda_2} \times
            \cdots \times S_{\lambda_\ell(\lambda)}

        induced to the group `S_{n}` is a symmetric function in the
        eigenvalues of a permutation matrix.  This basis is that character.

        It has the property that if the element indexed by the
        partition `\lambda` is evaluated at the roots of a permutation of
        cycle structure `\rho` then the value is the coefficient
        `\left< h_{(n-|\lambda|,\lambda)}, p_\rho \right>`.

        In terms of methods that are implemented in Sage, if ``n`` is
        a sufficiently large integer, then
        ``ht(lam).character_to_frobenius_image(n)`` is equal the complete
        function indexed by ``[n-sum(lam)]+lam``.

        This basis is introduced in [OZ2015]_.

        .. SEEALSO::

            :meth:`~sage.combinat.sf.sfa.SymmetricFunctionAlgebra_generic_Element.character_to_frobenius_image`,
            :meth:`~sage.combinat.sf.sfa.SymmetricFunctionAlgebra_generic_Element.eval_at_permutation_roots`

        EXAMPLES::

            sage: SymmetricFunctions(QQ).induced_trivial_character()
            Symmetric Functions over Rational Field in the induced trivial symmetric group character basis
            sage: ht = SymmetricFunctions(QQ).ht()
            sage: h = SymmetricFunctions(QQ).h()
            sage: h(ht([3,2]).character_to_frobenius_image(9))
            h[4, 3, 2]
            sage: h(ht([3,2]).character_to_frobenius_image(7))
            h[3, 2, 2]
            sage: h(ht([3,2]).character_to_frobenius_image(5))
            h[3, 2]
            sage: h(ht([3,2]).character_to_frobenius_image(4))
            0
            sage: p = SymmetricFunctions(QQ).p()
            sage: [h([4,1]).scalar(p(rho)) for rho in Partitions(5)]
            [0, 1, 0, 2, 1, 3, 5]
            sage: [ht([1]).eval_at_permutation_roots(rho) for rho in Partitions(5)]
            [0, 1, 0, 2, 1, 3, 5]
        """
        from .character import induced_trivial_character_basis
        return induced_trivial_character_basis(self, 'ht')

    ht = induced_trivial_character


    def forgotten(self):
        r"""
        The forgotten basis of the Symmetric Functions (or the basis dual to
        the elementary basis with respect to the Hall scalar product).

        EXAMPLES::

            sage: SymmetricFunctions(QQ).forgotten()
            Symmetric Functions over Rational Field in the forgotten basis

        TESTS:

        Over the rationals::

            sage: Sym = SymmetricFunctions(QQ)
            sage: e = Sym.e()
            sage: f = Sym.f()
            sage: h = Sym.h()
            sage: p = Sym.p()
            sage: s = Sym.s()
            sage: m = Sym.m()
            sage: e(f([2,1]))
            -2*e[1, 1, 1] + 5*e[2, 1] - 3*e[3]
            sage: f(e([2,1]))
            3*f[1, 1, 1] + 2*f[2, 1] + f[3]
            sage: h(f([2,1]))
            h[2, 1] - 3*h[3]
            sage: f(h([2,1]))
            3*f[1, 1, 1] + f[2, 1]
            sage: p(f([2,1]))
            -p[2, 1] - p[3]
            sage: f(p([2,1]))
            -f[2, 1] - f[3]
            sage: s(f([2,1]))
            s[2, 1] - 2*s[3]
            sage: f(s([2,1]))
            2*f[1, 1, 1] + f[2, 1]
            sage: m(f([2,1]))
            -m[2, 1] - 2*m[3]
            sage: f(m([2,1]))
            -f[2, 1] - 2*f[3]

        Over the integers::

            sage: Sym = SymmetricFunctions(ZZ)
            sage: e = Sym.e()
            sage: f = Sym.f()
            sage: h = Sym.h()
            sage: p = Sym.p()
            sage: s = Sym.s()
            sage: m = Sym.m()
            sage: e(f([2,1]))
            -2*e[1, 1, 1] + 5*e[2, 1] - 3*e[3]
            sage: f(e([2,1]))
            3*f[1, 1, 1] + 2*f[2, 1] + f[3]
            sage: h(f([2,1]))
            h[2, 1] - 3*h[3]
            sage: f(h([2,1]))
            3*f[1, 1, 1] + f[2, 1]
            sage: f(p([2,1]))
            -f[2, 1] - f[3]
            sage: s(f([2,1]))
            s[2, 1] - 2*s[3]
            sage: f(s([2,1]))
            2*f[1, 1, 1] + f[2, 1]
            sage: m(f([2,1]))
            -m[2, 1] - 2*m[3]
            sage: f(m([2,1]))
            -f[2, 1] - 2*f[3]

        Conversion from the forgotten basis to the power-sum basis over the
        integers is not well-defined in general, even if the result happens
        to have integral coefficients::

            sage: p(f([2,1]))
            Traceback (most recent call last):
            ...
            TypeError: no conversion of this rational to integer

        Fun exercise: prove that `p(f_{\lambda})` and `p(m_{\lambda})` have
        integral coefficients whenever `\lambda` is a strict partition.
        """
        return self.elementary().dual_basis()
    f = forgotten

    def symplectic(self):
        """
        The symplectic basis of the symmetric functions.

        .. SEEALSO:: :class:`~sage.combinat.sf.symplectic.SymmetricFunctionAlgebra_symplectic`

        EXAMPLES::

            sage: SymmetricFunctions(QQ).symplectic()
            Symmetric Functions over Rational Field in the symplectic basis
        """
        from . import symplectic
        return symplectic.SymmetricFunctionAlgebra_symplectic(self)
    sp = symplectic

    def orthogonal(self):
        """
        The orthogonal basis of the symmetric functions.

        .. SEEALSO:: :class:`~sage.combinat.sf.orthogonal.SymmetricFunctionAlgebra_orthogonal`

        EXAMPLES::

            sage: SymmetricFunctions(QQ).orthogonal()
            Symmetric Functions over Rational Field in the orthogonal basis
        """
        from . import orthogonal
        return orthogonal.SymmetricFunctionAlgebra_orthogonal(self)
    o = orthogonal

    def hecke_character(self, q='q'):
        """
        The basis of symmetric functions that determines the character
        tables for Hecke algebras.

        EXAMPLES::

            sage: SymmetricFunctions(ZZ['q'].fraction_field()).hecke_character()
            Symmetric Functions over
             Fraction Field of Univariate Polynomial Ring in q over Integer Ring
             in the Hecke character with q=q basis
            sage: SymmetricFunctions(QQ).hecke_character(1/2)
            Symmetric Functions over Rational Field in the Hecke character with q=1/2 basis
        """
        from sage.combinat.sf.hecke import HeckeCharacter
        return HeckeCharacter(self, q)
    qbar = hecke_character

    def macdonald(self, q='q', t='t'):
        r"""
        Returns the entry point for the various Macdonald bases.

        INPUT:

        - ``q``, ``t`` -- parameters

        Macdonald symmetric functions including bases `P`, `Q`, `J`, `H`, `Ht`.
        This also contains the `S` basis which is dual to the Schur basis with
        respect to the `q,t` scalar product.

        The parameters `q` and `t` must be in the base_ring of parent.

        EXAMPLES::

            sage: Sym = SymmetricFunctions(FractionField(QQ['q','t']))
            sage: P = Sym.macdonald().P(); P
            Symmetric Functions over Fraction Field of Multivariate Polynomial Ring in q, t over Rational Field in the Macdonald P basis
            sage: P[2]
            McdP[2]
            sage: Q = Sym.macdonald().Q(); Q
            Symmetric Functions over Fraction Field of Multivariate Polynomial Ring in q, t over Rational Field in the Macdonald Q basis
            sage: S = Sym.macdonald().S()
            sage: s = Sym.schur()
            sage: matrix([[S(la).scalar_qt(s(mu)) for la in Partitions(3)] for mu in Partitions(3)])
            [1 0 0]
            [0 1 0]
            [0 0 1]
            sage: H = Sym.macdonald().H()
            sage: s(H[2,2])
            q^2*s[1, 1, 1, 1] + (q^2*t+q*t+q)*s[2, 1, 1] + (q^2*t^2+1)*s[2, 2] + (q*t^2+q*t+t)*s[3, 1] + t^2*s[4]

            sage: Sym = SymmetricFunctions(QQ['z','q'].fraction_field())
            sage: (z,q) = Sym.base_ring().gens()
            sage: Hzq = Sym.macdonald(q=z,t=q).H()
            sage: H1z = Sym.macdonald(q=1,t=z).H()
            sage: s = Sym.schur()
            sage: s(H1z([2,2]))
            s[1, 1, 1, 1] + (2*z+1)*s[2, 1, 1] + (z^2+1)*s[2, 2] + (z^2+2*z)*s[3, 1] + z^2*s[4]
            sage: s(Hzq[2,2])
            z^2*s[1, 1, 1, 1] + (z^2*q+z*q+z)*s[2, 1, 1] + (z^2*q^2+1)*s[2, 2] + (z*q^2+z*q+q)*s[3, 1] + q^2*s[4]
            sage: s(H1z(Hzq[2,2]))
            z^2*s[1, 1, 1, 1] + (z^2*q+z*q+z)*s[2, 1, 1] + (z^2*q^2+1)*s[2, 2] + (z*q^2+z*q+q)*s[3, 1] + q^2*s[4]
        """
        return macdonald.Macdonald(self, q=q, t=t)

    def hall_littlewood(self, t='t'):
        """
        Returns the entry point for the various Hall-Littlewood bases.

        INPUT:

        - ``t`` -- parameter

        Hall-Littlewood symmetric functions including bases `P`, `Q`, `Qp`.
        The Hall-Littlewood `P` and `Q` functions at `t=-1` are the
        Schur-P and Schur-Q functions when indexed by strict partitions.

        The parameter `t` must be in the base ring of parent.

        EXAMPLES::

            sage: Sym = SymmetricFunctions(FractionField(QQ['t']))
            sage: P = Sym.hall_littlewood().P(); P
            Symmetric Functions over Fraction Field of Univariate Polynomial Ring in t over Rational Field in the Hall-Littlewood P basis
            sage: P[2]
            HLP[2]
            sage: Q = Sym.hall_littlewood().Q(); Q
            Symmetric Functions over Fraction Field of Univariate Polynomial Ring in t over Rational Field in the Hall-Littlewood Q basis
            sage: Q[2]
            HLQ[2]
            sage: Qp = Sym.hall_littlewood().Qp(); Qp
            Symmetric Functions over Fraction Field of Univariate Polynomial Ring in t over Rational Field in the Hall-Littlewood Qp basis
            sage: Qp[2]
            HLQp[2]
        """
        return hall_littlewood.HallLittlewood(self, t=t)

    def jack(self, t='t'):
        """
        Returns the entry point for the various Jack bases.

        INPUT:

        - ``t`` -- parameter

        Jack symmetric functions including bases `P`, `Q`, `Qp`.

        The parameter `t` must be in the base ring of parent.

        EXAMPLES::

            sage: Sym = SymmetricFunctions(FractionField(QQ['t']))
            sage: JP = Sym.jack().P(); JP
            Symmetric Functions over Fraction Field of Univariate Polynomial Ring in t over Rational Field in the Jack P basis
            sage: JQ = Sym.jack().Q(); JQ
            Symmetric Functions over Fraction Field of Univariate Polynomial Ring in t over Rational Field in the Jack Q basis
            sage: JJ = Sym.jack().J(); JJ
            Symmetric Functions over Fraction Field of Univariate Polynomial Ring in t over Rational Field in the Jack J basis
            sage: JQp = Sym.jack().Qp(); JQp
            Symmetric Functions over Fraction Field of Univariate Polynomial Ring in t over Rational Field in the Jack Qp basis
        """
        return jack.Jack( self, t=t )

    def zonal(self):
        """
        The zonal basis of the Symmetric Functions

        EXAMPLES::

            sage: SymmetricFunctions(QQ).zonal()
            Symmetric Functions over Rational Field in the zonal basis
        """
        return jack.SymmetricFunctionAlgebra_zonal( self )

    def llt(self, k, t='t'):
        """
        The LLT symmetric functions.

        INPUT:

        - ``k`` -- a positive integer indicating the level
        - ``t`` -- a parameter (default: `t`)

        LLT polynomials in `hspin` and `hcospin` bases.

        EXAMPLES::

            sage: llt3 = SymmetricFunctions(QQ['t'].fraction_field()).llt(3); llt3
            level 3 LLT polynomials over Fraction Field of Univariate Polynomial Ring in t over Rational Field
            sage: llt3.hspin()
            Symmetric Functions over Fraction Field of Univariate Polynomial Ring in t over Rational Field in the level 3 LLT spin basis
            sage: llt3.hcospin()
            Symmetric Functions over Fraction Field of Univariate Polynomial Ring in t over Rational Field in the level 3 LLT cospin basis
            sage: llt3.hcospin()
            Symmetric Functions over Fraction Field of Univariate Polynomial Ring in t over Rational Field in the level 3 LLT cospin basis
        """
        return llt.LLT_class( self, k, t=t )

    def from_polynomial(self, f):
        """
        Converts a symmetric polynomial ``f`` to a symmetric function.

        INPUT:

        - ``f`` -- a symmetric polynomial

        This function converts a symmetric polynomial `f` in a polynomial ring in finitely
        many variables to a symmetric function in the monomial
        basis of the ring of symmetric functions over the same base ring.

        EXAMPLES::

            sage: P = PolynomialRing(QQ, 'x', 3)
            sage: x = P.gens()
            sage: f = x[0] + x[1] + x[2]
            sage: S = SymmetricFunctions(QQ)
            sage: S.from_polynomial(f)
            m[1]

            sage: f = x[0] + 2*x[1] + x[2]
            sage: S.from_polynomial(f)
            Traceback (most recent call last):
            ...
            ValueError: x0 + 2*x1 + x2 is not a symmetric polynomial
        """
        return self.m().from_polynomial(f)

    def register_isomorphism(self, morphism, only_conversion=False):
        r"""
        Register an isomorphism between two bases of ``self``, as a canonical coercion
        (unless the optional keyword ``only_conversion`` is set to ``True``,
        in which case the isomorphism is registered as conversion only).

        EXAMPLES:

        We override the canonical coercion from the Schur basis to the
        powersum basis by a (stupid!) map `s_\lambda\mapsto 2p_\lambda`.
        ::

            sage: Sym = SymmetricFunctions(QQ['zorglub']) # make sure we are not going to screw up later tests
            sage: s = Sym.s(); p = Sym.p().dual_basis()
            sage: phi = s.module_morphism(diagonal = lambda t: 2, codomain = p)
            sage: phi(s[2, 1])
            2*d_p[2, 1]
            sage: Sym.register_isomorphism(phi)
            sage: p(s[2,1])
            2*d_p[2, 1]

        The map is supposed to implement the canonical isomorphism
        between the two bases. Otherwise, the results will be
        mathematically wrong, as above. Use with care!
        """
        if only_conversion:
            morphism.codomain().register_conversion(morphism)
        else:
            morphism.codomain().register_coercion(morphism)

    _shorthands = ['e', 'f', 'h', 'm', 'p', 's']
    _shorthands_all = sorted(_shorthands + ['ht', 'o', 'sp', 'st', 'w'])

    def __init_extra__(self):
        """
        Sets up the coercions between the different bases

        EXAMPLES::

            sage: Sym = SymmetricFunctions(QQ) # indirect doctest
            sage: s = Sym.s(); p = Sym.p()
            sage: f = s.coerce_map_from(p); f
            Generic morphism:
              From: Symmetric Functions over Rational Field in the powersum basis
              To:   Symmetric Functions over Rational Field in the Schur basis
            sage: p.an_element()
            2*p[] + 2*p[1] + 3*p[2]
            sage: f(p.an_element())
            2*s[] + 2*s[1] - 3*s[1, 1] + 3*s[2]
            sage: f(p.an_element()) == p.an_element()
            True

        """
        #powersum   = self.powersum  ()
        #complete   = self.complete  ()
        #elementary = self.elementary()
        #schur      = self.schur     ()
        #monomial   = self.monomial  ()

        iso = self.register_isomorphism

        from sage.combinat.sf.classical import conversion_functions

        for (basis1_name, basis2_name) in conversion_functions:
            basis1 = getattr(self, basis1_name)()
            basis2 = getattr(self, basis2_name)()
            on_basis = SymmetricaConversionOnBasis(t = conversion_functions[basis1_name,basis2_name], domain = basis1, codomain = basis2)
            from sage.rings.rational_field import RationalField
            if basis2_name != "powersum" or self._base.has_coerce_map_from(RationalField()):
                iso(basis1._module_morphism(on_basis, codomain = basis2))
            else:
                # Don't register conversions to powersums as coercions,
                # unless the base ring is a `\QQ`-algebra
                # (otherwise the coercion graph loses commutativity).
                iso(basis1._module_morphism(on_basis, codomain = basis2), only_conversion = True)

        # Todo: fill in with other conversion functions on the classical bases

    def kBoundedSubspace(self, k, t='t'):
        r"""
        Return the `k`-bounded subspace of the ring of symmetric functions.

        INPUT:

        - ``k`` - a positive integer
        - ``t`` a formal parameter; `t=1` yields a subring

        The subspace of the ring of symmetric functions spanned by
        `\{ s_{\lambda}[X/(1-t)] \}_{\lambda_1\le k} = \{ s_{\lambda}^{(k)}[X,t]\}_{\lambda_1 \le k}`
        over the base ring `\QQ[t]`. When `t=1`, this space is in fact a subalgebra of
        the ring of symmetric functions generated by the complete homogeneous symmetric functions
        `h_i` for `1\le i \le k`.

        .. SEEALSO:: :meth:`sage.combinat.sf.new_kschur.KBoundedSubspace`

        EXAMPLES::

            sage: Sym = SymmetricFunctions(QQ)
            sage: KB = Sym.kBoundedSubspace(3,1); KB
            3-bounded Symmetric Functions over Rational Field with t=1

            sage: Sym = SymmetricFunctions(QQ['t'])
            sage: Sym.kBoundedSubspace(3)
            3-bounded Symmetric Functions over Univariate Polynomial Ring in t over Rational Field

            sage: Sym = SymmetricFunctions(QQ['z'])
            sage: z = Sym.base_ring().gens()[0]
            sage: Sym.kBoundedSubspace(3,t=z)
            3-bounded Symmetric Functions over Univariate Polynomial Ring in z over Rational Field with t=z
        """
        from sage.combinat.sf.new_kschur import KBoundedSubspace
        return KBoundedSubspace(self, k, t=t)

    def kschur(self, k, t='t'):
        r"""
        Returns the `k`-Schur functions.

        EXAMPLES::

            sage: Sym = SymmetricFunctions(QQ)
            sage: ks = Sym.kschur(3,1)
            sage: ks[2]*ks[2]
            ks3[2, 2] + ks3[3, 1]
            sage: ks[2,1,1].lift()
            s[2, 1, 1] + s[3, 1]

            sage: Sym = SymmetricFunctions(QQ['t'])
            sage: ks = Sym.kschur(3)
            sage: ks[2,2,1].lift()
            s[2, 2, 1] + t*s[3, 2]
        """
        return self.kBoundedSubspace(k, t=t).kschur()

    def ksplit(self, k, t='t'):
        r"""
        Return the `k`-split basis of the `k`-bounded subspace.

        EXAMPLES::

            sage: Sym = SymmetricFunctions(QQ)
            sage: ksp = Sym.ksplit(3,1)
            sage: ksp[2]*ksp[2]
            ksp3[2, 2] + ksp3[3, 1]
            sage: ksp[2,1,1].lift()
            s[2, 1, 1] + s[2, 2] + s[3, 1]

            sage: Sym = SymmetricFunctions(QQ['t'])
            sage: ksp = Sym.ksplit(3)
            sage: ksp[2,1,1].lift()
            s[2, 1, 1] + t*s[2, 2] + t*s[3, 1]
        """
        return self.kBoundedSubspace(k, t=t).ksplit()

    def khomogeneous(self, k):
        r"""
        Returns the homogeneous symmetric functions in the `k`-bounded subspace.

        EXAMPLES::

            sage: Sym = SymmetricFunctions(QQ)
            sage: kh = Sym.khomogeneous(4)
            sage: kh[3]*kh[4]
            h4[4, 3]
            sage: kh[4].lift()
            h[4]
        """
        return self.kBoundedSubspace(k, t=1).khomogeneous()

    def kBoundedQuotient(self, k, t='t'):
        r"""
        Returns the `k`-bounded quotient space of the ring of symmetric functions.

        INPUT:

        - ``k`` - a positive integer

        The quotient of the ring of symmetric functions ...

        .. SEEALSO:: :meth:`sage.combinat.sf.k_dual.KBoundedQuotient`

        EXAMPLES::

            sage: Sym = SymmetricFunctions(QQ)
            sage: KQ = Sym.kBoundedQuotient(3); KQ
            Traceback (most recent call last):
            ...
            TypeError: unable to convert 't' to a rational
            sage: KQ = Sym.kBoundedQuotient(3,t=1); KQ
            3-Bounded Quotient of Symmetric Functions over Rational Field with t=1
            sage: Sym = SymmetricFunctions(QQ['t'].fraction_field())
            sage: KQ = Sym.kBoundedQuotient(3); KQ
            3-Bounded Quotient of Symmetric Functions over Fraction Field of Univariate Polynomial Ring in t over Rational Field
        """
        from sage.combinat.sf.k_dual import KBoundedQuotient
        return KBoundedQuotient(self, k, t)

class SymmetricaConversionOnBasis:
    def __init__(self, t, domain, codomain):
        """
        Initialization of ``self``.

        INPUT:

        - ``t`` -- a function taking a monomial in CombinatorialFreeModule(QQ, Partitions()),
           and returning a (partition, coefficient) list.

        - ``domain``, ``codomain`` -- parents

        Construct a function mapping a partition to an element of ``codomain``.

        This is a temporary quick hack to wrap around the existing
        symmetrica conversions, without changing their specs.

        EXAMPLES::

            sage: Sym = SymmetricFunctions(QQ['x'])
            sage: p = Sym.p(); s = Sym.s()
            sage: def t(x): [(p,c)] = x; return [ (p,2*c), (p.conjugate(), c) ]
            sage: f = sage.combinat.sf.sf.SymmetricaConversionOnBasis(t, p, s)
            sage: f(Partition([3,1]))
            s[2, 1, 1] + 2*s[3, 1]
        """
        self._domain = domain
        self.fake_sym = CombinatorialFreeModule(QQ, Partitions())
        self._codomain = codomain
        self._t = t

    def __call__(self, partition):
        """
            sage: Sym = SymmetricFunctions(QQ['x'])
            sage: p = Sym.p(); s = Sym.s()
            sage: p[1] + s[1]                           # indirect doctest
            2*p[1]
        """
        # TODO: use self._codomain.sum_of_monomials, when the later
        # will have an optional optimization for the case when there
        # is no repetition in the support
        return self._codomain._from_dict(dict(self._t(self.fake_sym.monomial(partition))), coerce = True)

