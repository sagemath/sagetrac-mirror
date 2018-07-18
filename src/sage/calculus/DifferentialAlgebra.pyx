r"""
Differential elimination

The DifferentialAlgebra module implements elimination methods for systems
of polynomial differential equations. Mathematically, it relies on the
differential algebra theory of Joseph Fels Ritt and Ellis Robert Kolchin.



AUTHORS:

- Nicolas M. Thiery (2011-12): initial version

- François Boulier (2012): the only other developer, so far

- Charles Bouillaguet (2012): helped with SAGE integration

- Brent Baccala (2018): updated to Sage 8.2


REFERENCES:

.. [Ritt50] Joseph Fels Ritt. Differential Algebra. Dover Publications Inc. New York, 1950.
.. [Kolchin73] Ellis Robert Kolchin. Differential Algebra and Algebraic Groups. Academic Press. New York, 1973.

Chemical Reaction Systems
~~~~~~~~~~~~~~~~~~~~~~~~~

EXAMPLES:

    This example shows how to build the Henri Michaelis Menten formula
    by differential elimination. One considers a chemical reaction
    system describing the enzymatic reaction::

                   k(1)
        E + S  -----------> ES
                   k(-1)
        ES     -----------> E + S
                   k(2)
        ES     -----------> E + P

    A substrate S is transformed into a product P, in the presence
    of an enzyme E. An intermediate complex ES is formed::

        sage: from sage.calculus.DifferentialAlgebra import DifferentialRing, RegularDifferentialChain, BaseFieldExtension
        sage: leader,order,rank = var ('leader,order,rank')
        sage: derivative = function ('derivative')
        sage: t = var('t')
        sage: k,F_1,E,S,ES,P = function('k,F_1,E,S,ES,P')
        sage: params = [k(-1),k(1),k(2)]
        sage: params
        [k(-1), k(1), k(2)]

    The main assumption is that `k(1), k(-1) >> k(2)` i.e. that the
    revertible reaction is much faster than the last one. One performs
    a quasi-steady state approximation by considering the following
    differential-algebraic system (it comes from the mass-action law
    kinetics, replacing the contribution of the fast reactions by an
    unknown function `F_1(t)`, on the algebraic variety where the fast
    reaction would equilibrate if they were alone)::

        sage: syst = [diff(E(t),t)   == - F_1(t) + k(2)*ES(t),
        ....:         diff(S(t),t)   == - F_1(t),
        ....:         diff (ES(t),t) == - k(2)*ES(t) + F_1(t),
        ....:         diff (P(t),t)  == k(2)*ES(t),
        ....:         k(1)*E(t)*S(t) - k(-1)*ES(t) == 0]
        sage: syst
        [diff(E(t), t) == ES(t)*k(2) - F_1(t),
         diff(S(t), t) == -F_1(t),
         diff(ES(t), t) == -ES(t)*k(2) + F_1(t),
         diff(P(t), t) == ES(t)*k(2),
         E(t)*S(t)*k(1) - ES(t)*k(-1) == 0]

    Differential elimination permits to simplify this DAE.
    To avoid discussing the possible vanishing of ``params``, one
    moves them to the base field of the equations::

        sage: Field = BaseFieldExtension (generators = params)
        sage: Field
        differential_field
        sage: R = DifferentialRing (derivations = [t], blocks = [F_1, [E,ES,P,S], params], parameters = params)
        sage: R
        Differential Ring over Rational Field
            with dependents F_1(t), E(t), ES(t), P(t), S(t),
                 independent t, and parameters k(-1), k(1), k(2)

    The Rosenfeld-Groebner algorihtm considers three cases. The two last ones
    are degenerate cases::

        sage: ideal = R.RosenfeldGroebner (syst, basefield = Field)
        sage: ideal
        [Regular differential chain (E(t)*S(t)*k(1) - ES(t)*k(-1),
             ES(t)*S(t)^2*k(2)*k(1) + ES(t)*S(t)*k(2)*k(-1) + S(t)^2*k(1)*diff(S(t), t)
                 + ES(t)*k(-1)*diff(S(t), t) + S(t)*k(-1)*diff(S(t), t),
            -ES(t)*k(2) + diff(P(t), t),
             ES(t)^2*k(2)*k(-1) + S(t)^2*k(1)*diff(ES(t), t) + ES(t)*k(-1)*diff(ES(t), t)
                 + S(t)*k(-1)*diff(ES(t), t),
            -ES(t)*S(t)^2*k(2)*k(1) + F_1(t)*S(t)^2*k(1) - ES(t)*S(t)*k(2)*k(-1)
                 + ES(t)*F_1(t)*k(-1) + F_1(t)*S(t)*k(-1)),
         Regular differential chain (S(t)*k(1) + k(-1), ES(t), E(t), diff(P(t), t), F_1(t)),
         Regular differential chain (S(t), ES(t), diff(P(t), t), diff(E(t), t), F_1(t))]
        sage: [ C.equations (solved = true) for C in ideal ]
        [[E(t) == ES(t)*k(-1)/(S(t)*k(1)),
          diff(S(t), t) == -(ES(t)*S(t)^2*k(2)*k(1) + ES(t)*S(t)*k(2)*k(-1))/(S(t)^2*k(1) + ES(t)*k(-1) + S(t)*k(-1)),
          diff(P(t), t) == ES(t)*k(2),
          diff(ES(t), t) == -ES(t)^2*k(2)*k(-1)/(S(t)^2*k(1) + ES(t)*k(-1) + S(t)*k(-1)),
          F_1(t) == (ES(t)*S(t)^2*k(2)*k(1) + ES(t)*S(t)*k(2)*k(-1))/(S(t)^2*k(1) + ES(t)*k(-1) + S(t)*k(-1))],
         [S(t) == -k(-1)/k(1), ES(t) == 0, E(t) == 0, diff(P(t), t) == 0, F_1(t) == 0],
         [S(t) == 0, ES(t) == 0, diff(P(t), t) == 0, diff(E(t), t) == 0, F_1(t) == 0]]

    The sought equation, below, is not yet the Henri-Michaelis-Menten
    formula. This is expected, since some minor hypotheses have not yet
    been taken into account::

        sage: ideal [0].equations (solved = true, selection = leader == derivative (S(t)))
        [diff(S(t), t) == -(ES(t)*S(t)^2*k(2)*k(1) + ES(t)*S(t)*k(2)*k(-1))/(S(t)^2*k(1) + ES(t)*k(-1) + S(t)*k(-1))]

    Let us take them into account. First create two new constants.
    Put them among ``params``, together with initial values::

        sage: K,V_max = var ('K,V_max')
        sage: params = [k(-1),k(1),k(2),E(0),ES(0),P(0),S(0),K,V_max]
        sage: params
        [k(-1), k(1), k(2), E(0), ES(0), P(0), S(0), K, V_max]

        sage: R = DifferentialRing (blocks = [F_1, [ES,E,P,S], params], parameters = params, derivations = [t])
        sage: R
        Differential Ring over Rational Field
            with dependents F_1(t), ES(t), E(t), P(t), S(t),
                 independent t,
             and parameters k(-1), k(1), k(2), E(0), ES(0), P(0), S(0), K, V_max

    There are relations among the parameters: initial values supposed
    to be zero, and equations meant to rename constants::

        sage: relations_among_params = RegularDifferentialChain ([P(0) == 0, ES(0) == 0, K == k(-1)/k(1), V_max == k(2)*E(0)], R)
        sage: relations_among_params
        Regular differential chain (P(0), ES(0), -E(0)*k(2) + V_max, K*k(1) - k(-1))

    Coming computations will be performed over a base field defined
    by generators and relations::

        sage: Field = BaseFieldExtension (generators = params, relations = relations_among_params)
        sage: Field
        differential_field

    Extend the DAE with linear conservation laws. They could have
    been computed from the stoichimetry matrix of the chemical system::

        sage: newsyst = syst + [E(t) + ES(t) == E(0) + ES(0), S(t) + ES(t) + P(t) == S(0) + ES(0) + P(0)]
        sage: newsyst
        [diff(E(t), t) == ES(t)*k(2) - F_1(t),
         diff(S(t), t) == -F_1(t),
         diff(ES(t), t) == -ES(t)*k(2) + F_1(t),
         diff(P(t), t) == ES(t)*k(2),
         E(t)*S(t)*k(1) - ES(t)*k(-1) == 0,
         E(t) + ES(t) == E(0) + ES(0),
         ES(t) + P(t) + S(t) == ES(0) + P(0) + S(0)]

    Simplify again. Only one case is left::

        sage: ideal = R.RosenfeldGroebner (newsyst, basefield = Field)
        sage: ideal
        [Regular differential chain (P(0), ES(0), E(0)*k(2) - V_max, -K*k(1) + k(-1),
             K*P(t) - K*S(0) + K*S(t) + E(0)*S(t) + P(t)*S(t) - S(0)*S(t) + S(t)^2,
            -K*E(0) + K*E(t) + E(t)*S(t), K*ES(t) - E(0)*S(t) + ES(t)*S(t),
             K*V_max*S(t) + V_max*S(t)^2 + K^2*diff(S(t), t) + K*E(0)*diff(S(t), t)
                 + 2*K*S(t)*diff(S(t), t) + S(t)^2*diff(S(t), t),
             K^2*F_1(t) + K*E(0)*F_1(t) - K*V_max*S(t) + 2*K*F_1(t)*S(t) - V_max*S(t)^2 + F_1(t)*S(t)^2)]

    To get the traditional Henri-Michaelis-Menten formula, one still needs
    to neglect the term ``K*E(0)``::

        sage: ideal[0].equations (solved = true, selection = leader == derivative (S(t)))
        [diff(S(t), t) == -(K*V_max*S(t) + V_max*S(t)^2)/(K^2 + K*E(0) + 2*K*S(t) + S(t)^2)]

    One can also get it by computing the right hand side of the equation
    which gives the evolution of the product `P`::

        sage: ideal[0].normal_form (diff(P(t),t))
        V_max*S(t)/(K + S(t))

REFERENCES:

..  [BLLM11] François Boulier, Marc Lefranc, François Lemaire and
    Pierre-Emmanuel Morant. Model Reduction of Chemical Reaction
    Systems using Elimination. Mathematics in Computer Science vol. 5,
    pp. 289-301. 2011. http://hal.archives-ouvertes.fr/hal-00184558

I/O Relations
~~~~~~~~~~~~~

EXAMPLE:

    This example shows how to compute an input/output relation, having
    a nice form, for a compartmental model. The degradation of the
    studied product, from compartment 1, is supposed to follow a
    Henri-Michaelis-Menten formula. The exchanges between the two
    compartments are supposed to be linear::

                          k(e),V(e)
        compartment 1  --------------> outside the model

                            k(12)
        compartment 1 ---------------> compartment 2

                            k(21)
        compartment 2 ---------------> compartment 1

    Functions `x_1(t)` and `x_2(t)` are associated to both compartments::

        sage: from sage.calculus.DifferentialAlgebra import DifferentialRing, RegularDifferentialChain, BaseFieldExtension
        sage: leader,order,rank = var ('leader,order,rank')
        sage: derivative = function ('derivative')

        sage: t = var ('t')
        sage: k,V,x1,x2 = function ('k,V,x1,x2')
        sage: params = [k(e),V(e),k(12),k(21)]
        sage: params
        [k(e), V(e), k(12), k(21)]

    The dynamical system associated to the above compartmental model
    already is a regular differential chain w.r.t. any orderly ranking,
    such as the one of ``R``::

        sage: R = DifferentialRing (derivations = [t], blocks = [[x1,x2], params], parameters = params)
        sage: R
        Differential Ring over Rational Field with dependents x1(t), x2(t),
            independent t, and parameters k(e), V(e), k(12), k(21)

    Here are the equations which give the dynamics of the system::

        sage: edoA = diff (x1(t),t) == -k(12)*x1(t) + k(21)*x2(t) - (V(e)*x1(t))/(k(e) + x1(t))
        sage: edoA
        diff(x1(t), t) == -k(12)*x1(t) + k(21)*x2(t) - V(e)*x1(t)/(k(e) + x1(t))

        sage: edoB = diff (x2(t),t) == k(12)*x1(t) - k(21)*x2(t)
        sage: edoB
        diff(x2(t), t) == k(12)*x1(t) - k(21)*x2(t)

    The parameters are moved to the base field of the equations, in
    order to avoid discussing their possible vanishing::

        sage: F = BaseFieldExtension (generators = params)
        sage: F
        differential_field

    The following computation does not do anything. The two equations
    are just bundled into a regular differential chain::

        sage: ideal = R.RosenfeldGroebner ([edoA, edoB], basefield = F)
        sage: ideal
        [Regular differential chain (-k(12)*x1(t) + k(21)*x2(t) + diff(x2(t), t),
             k(12)*k(e)*x1(t) + k(12)*x1(t)^2 - k(21)*k(e)*x2(t) - k(21)*x1(t)*x2(t)
                 + V(e)*x1(t) + k(e)*diff(x1(t), t) + x1(t)*diff(x1(t), t))]
        sage: C = ideal[0]
        sage: C.equations (solved = true)
        [diff(x2(t), t) == k(12)*x1(t) - k(21)*x2(t), diff(x1(t), t) == -(k(12)*k(e)*x1(t) + k(12)*x1(t)^2 - k(21)*k(e)*x2(t) - k(21)*x1(t)*x2(t) + V(e)*x1(t))/(k(e) + x1(t))]

    Let us assume now that compartment 1 is observed while compartment 2
    is not. In other words, let us assume that the output `y = x_1(t)`.
    The idea is to eliminate the non observed variable `x_2(t)` and
    compute a relation (the input/output equation), which is a consequence
    of the dynamical system, but does not involve `x_2(t)` and its
    derivatives. For this purpose, one defines a new differential ring,
    which is mathematically equivalent to `R`, but with a different ranking::

        sage: IO_R = DifferentialRing (derivations = [t], blocks = [x2,x1,params], parameters = params)
        sage: IO_R
        Differential Ring over Rational Field with dependents x2(t), x1(t),
            independent t, and parameters k(e), V(e), k(12), k(21)

    One just has to perform a change of ranking over `C`::

        sage: IO_C = C.change_ranking (IO_R)
        sage: IO_C
        Regular differential chain (V(e)*k(21)*k(e)*x1(t) + V(e)*k(21)*x1(t)^2
                + k(21)*k(e)^2*diff(x1(t), t) + k(12)*k(e)^2*diff(x1(t), t)
                + 2*k(21)*k(e)*x1(t)*diff(x1(t), t) + 2*k(12)*k(e)*x1(t)*diff(x1(t), t)
                + k(21)*x1(t)^2*diff(x1(t), t) + k(12)*x1(t)^2*diff(x1(t), t)
                + V(e)*k(e)*diff(x1(t), t) + k(e)^2*diff(x1(t), t, t)
                + 2*k(e)*x1(t)*diff(x1(t), t, t) + x1(t)^2*diff(x1(t), t, t),
           -k(12)*k(e)*x1(t) - k(12)*x1(t)^2 + k(21)*k(e)*x2(t) + k(21)*x1(t)*x2(t)
                - V(e)*x1(t) - k(e)*diff(x1(t), t) - x1(t)*diff(x1(t), t))

    Here is the input/output equation. However, it looks quite complicated::

        sage: IO_rel = IO_C.equations (selection = leader == derivative(x1(t)))[0]
        sage: IO_rel
        V(e)*k(21)*k(e)*x1(t) + V(e)*k(21)*x1(t)^2 + k(21)*k(e)^2*diff(x1(t), t) + k(12)*k(e)^2*diff(x1(t), t) + 2*k(21)*k(e)*x1(t)*diff(x1(t), t) + 2*k(12)*k(e)*x1(t)*diff(x1(t), t) + k(21)*x1(t)^2*diff(x1(t), t) + k(12)*x1(t)^2*diff(x1(t), t) + V(e)*k(e)*diff(x1(t), t) + k(e)^2*diff(x1(t), t, t) + 2*k(e)*x1(t)*diff(x1(t), t, t) + x1(t)^2*diff(x1(t), t, t)

    One way to simplify it consists in integrating it. For this purpose,
    one divides it by its initial and then one and then one integrates it::

        sage: IO_rel = IO_rel / IO_R.initial (IO_rel)
        sage: IO_rel
        (V(e)*k(21)*k(e)*x1(t) + V(e)*k(21)*x1(t)^2 + k(21)*k(e)^2*diff(x1(t), t) + k(12)*k(e)^2*diff(x1(t), t) + 2*k(21)*k(e)*x1(t)*diff(x1(t), t) + 2*k(12)*k(e)*x1(t)*diff(x1(t), t) + k(21)*x1(t)^2*diff(x1(t), t) + k(12)*x1(t)^2*diff(x1(t), t) + V(e)*k(e)*diff(x1(t), t) + k(e)^2*diff(x1(t), t, t) + 2*k(e)*x1(t)*diff(x1(t), t, t) + x1(t)^2*diff(x1(t), t, t))/(k(e)^2 + 2*k(e)*x1(t) + x1(t)^2)

    Then, one integrates it. Simpler, isn't it::

        sage: L = IO_R.integrate (IO_rel, t)
        sage: L
        [V(e)*k(21)*x1(t)/(k(e) + x1(t)),
         -(k(21)*k(e)^2 + k(12)*k(e)^2 - k(21)*x1(t)^2 - k(12)*x1(t)^2 + V(e)*k(e))/(k(e) + x1(t)),
         x1(t)]

    The precise relationship between the input/output equation and the
    list `L` is given by the next formula::

        sage: zero = add (IO_R.differentiate (L[i], t^i) for i in range (len (L))) - IO_rel

    The variable ``zero`` contains `0` but sage does not recognize it.
    Let us force the simplification::

        sage: IO_R.normal_form (zero)
        0
"""

#*****************************************************************************
#  Copyright (C) 2011-2012 Nicolas M. Thiery <nthiery at users.sf.net>
#                2011-2012 François Boulier <Francois.Boulier at univ-lille1.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

# To avoid eval ('1/7') giving 0, use eval (preparse (...))
from sage.repl.preparse import preparse
from sage.calculus.var import function
from sage.all import Integer, flatten, derivative, latex
from sage.symbolic.expression import Expression
from sage.symbolic.function_factory import function as symbolic_function

# __all__ = (DifferentialRing, RegularDifferentialChain, BaseFieldExtension)

import sys
import string

from cpython cimport bool

# blad_c.pxd is automatically built from blad.h
# bmi_c.pxd is automatically built from bmi_balsa.h
# bmi_strings.pyx is automatically built from bmi_indices-h.in

cimport sage.libs.blad_c  as blad_c
cimport sage.libs.bmi_c as bmi_c

include "sage/libs/bmi_strings.pyx"

import re
import itertools

# DifferentialRing
cdef class DifferentialRing:
    r"""
    The class DifferentialRing implements differential polynomial rings

        Differential rings can be endowed with zero, one or many different
        derivations.

        Each derivation is denoted using a symbol and is associated to
        an independent variable (for each derivation `x`, it is possible
        to differentiate w.r.t. `x`, and to manipulate differential polynomials
        depending on `x`).

        Differential rings are built over an alphabet of dependent variables
        (called differential indeterminates by Ritt and Kolchin), which
        are presented in a list of blocks.

        Differential polynomials are plain polynomials, built over the
        infinite alphabet of the derivatives of the dependent
        variables, often called *derivatives* in this module, plus the
        finite set of the independent variables.

        While building differential rings, it is possible to customize
        the dependencies of the dependent variables w.r.t. the independent
        variables. By default, each dependent variable depends on all
        the independent variables, following the order of the list of
        derivations. It is possible to customize the order of the
        independent variables. It is also possible to specify the list
        of the independent variables the dependent variables depend on.
        This feature is especially useful to define parameters, i.e.
        dependent variables which do not depend on any independent variables.

        The list of blocks, together with the list of derivations, define
        a so-called *ranking*. Observe that rankings apply to parameters also
        and that parameters do not need to lie at the bottom of the rankings.

        A ranking is any total ordering over the set of the derivatives (in
        this package, rankings are extended to the independent variables),
        which satisfies the two axioms of rankings:

        -- Every derivative `u` is less than any of its proper derivatives.

        -- If `u` and `v` are derivatives such that `u` is less than `v` and `x` is
           any independent variable, then the derivative of `u` with respect
           to `x` is less than the derivative of `v` with respect to `x`.

        In this package, rankings are defined by the list ``x[1] > ... > x[p]``
        of the independent variables plus a list ``b[1] >> ... >> b[n]`` of blocks.
        Each block b is defined by a list ``u[1] > ... > u[m]`` of dependent
        variables. Any dependent variable must appear in exactly one block.

        The ``>>`` operator between blocks indicates a block elimination
        ranking: if ``b[i] >> b[j]`` are two blocks, ``v[i]`` is any derivative of
        any dependent variable occuring in ``b[i]``, and, ``v[j]`` is any derivative
        of any dependent variable occuring in ``b[j]``, then ``v[i] > v[j]``.

        Within a given block ``b = u[1] > ... > u[m]``, in the ordinary
        differential case (only one derivation), the derivatives are ordered
        by the unique orderly ranking such that ``u[1] > ... > u[m]``:
        the derivatives of higher order are ranked above the derivatives
        if lower order ; two derivatives having the same order are
        ranked using the ordering defined by the list.

    NOTES:

        Rankings for partial differential equations are not yet fully
        supported, due to limitations of Sage functions.

    EXAMPLES:

    The examples focus on the concept of rankings and parameters.
    See the help page of :mod:`~sage.calculus.DifferentialAlgebra` for relevant examples::

        sage: from sage.calculus.DifferentialAlgebra import DifferentialRing
        sage: x,a,b = var ('x,a,b')
        sage: y,z = function ('y,z')

    Mathematically, ``R`` is a differential ring, with a single derivation.
    Among the four dependent variables listed in the block list, two
    of them do not depend on any independent variable. They thus should
    be viewed as parameters, in the usual sense. Following Ritt and Kolchin
    notation, `R = Q[x,a,b]\{z,y\}` with derivation `d/dx`::

        sage: R = DifferentialRing(derivations = [x], blocks = [z,y,a,b], parameters = [a,b])
        sage: R
        Differential Ring over Rational Field with dependents z(x), y(x),
            independent x, and parameters a, b

    Here is a differential polynomial which belongs to `R`::

        sage: poly = a*diff(y(x),x) + b + x*y(x) + diff(z(x),x)^2 + z(x)

    The ranking is a pure elimination ranking::

        sage: R.sort( R.indets(poly, selection = 'all'), 'descending' )
        [diff(z(x), x), z(x), diff(y(x), x), y(x), a, b, x]

    Mathematically, `S` (defined below) and `R` define the same ring. The ranking of `S` is
    different from that of `R`. Observe that parameters do not lie at the
    bottom of the ranking::

        sage: S = DifferentialRing( derivations = [x], blocks = [b,a,y,z], parameters = [a,b] )
        sage: S
        Differential Ring over Rational Field with dependents y(x), z(x),
            independent x, and parameters b, a
        sage: S.sort( R.indets(poly, selection = 'all'), 'descending' )
        [b, a, diff(y(x), x), y(x), diff(z(x), x), z(x), x]

    Mathematically, `T` (defined below), `S` and `R` define the same ring. The derivatives of
    z and y are ranked orderly::

        sage: T = DifferentialRing( derivations = [x], blocks = [[z,y],a,b], parameters = [a,b] )
        sage: T
        Differential Ring over Rational Field with dependents z(x), y(x),
            independent x, and parameters a, b
        sage: T.sort( R.indets(poly, selection = 'all'), 'descending' )
        [diff(z(x), x), diff(y(x), x), z(x), y(x), a, b, x]

    Last, here is an example of a PDE ring. The ranking is orderly.
    The dependent variable v does not depend on x::

        sage: x,y = var('x,y')
        sage: u,v = function('u,v')
        sage: R = DifferentialRing(derivations = [x,y], blocks = [[u,v]], parameters = [v(y)])
        sage: poly = u(x,y) + v(y)
        sage: R.sort( R.indets(poly, selection = 'all'), 'descending' )
        [u(x, y), v(y)]
    """

    cdef bmi_c.ALGEB dring
    """ Pointer to the BALSA internal data structure """

    cdef dict _names
    """ Dictionary mapping BMI names to Sage objects, for parsing """

# __init__
    def __init__ (
                self,
                list derivations = [],
                list blocks = [],
                list parameters = [],
                char* notation = BMI_IX_D):
        r"""
        Constructor of DifferentialRing

        INPUT:

        - ``derivations`` -- a list of variables
        - ``blocks``      -- a list of blocks, each block being either
                             a dependent variable or, a list of
                             dependent variables
        - ``parameters``  -- a list of dependent variables, with
                             special dependencies
        - ``notation``   -- (default: 'undefined') a string

        OUTPUT:

        Create a new differential ring with derivations (independent
        variables) ``derivations``, differential indeterminates (dependent
        variables) ``blocks``, parameters ``parameters`` and default
        notation ``notation``.

        The ranking is mostly defined by the order of the blocks
        in ``blocks``. The leftmost blocks are ranked higher than
        the rightmost ones. Dependent variables which belong to a
        same block are ranked orderly.

        Parameters are dependent variables which do not need to
        depend on all independent variables, or have special dependencies.
        This definition generalizes parameters, in the usual sense:
        dependent variables which do not depend on any independent
        variables. Parameters do not need to lie at the bottom of
        the ranking.

        See the help page of :class:`~sage.calculus.DifferentialAlgebra.DifferentialRing`
        for more details.

        NOTES:

        It would be more consistent to define parameters as function.
        They need sometimes to be defined as var because of limitations
        of function.

        So far, there is a single notation available, because of
        limitations of function.

        EXAMPLES::

            sage: from sage.calculus.DifferentialAlgebra import DifferentialRing
            sage: x,a,b = var ('x,a,b')
            sage: y,z = function ('y,z')
            sage: R = DifferentialRing (derivations = [x], blocks = [z,y,a,b], parameters = [a,b])
            sage: R
            Differential Ring over Rational Field with dependents z(x), y(x),
                independent x, and parameters a, b
            sage: latex(R)
            \Bold{Q}[x]\{z\left(x\right),y\left(x\right),a,b\}

        """
#
        cdef bmi_c.ALGEB A
        cdef bytes mesgerr, strders, strblks, strpars
        strders = bytes (derivations)
        strblks = bytes (blocks)
        strpars = bytes (parameters)
# First build the BMI differential ring and exit on error
        A = bmi_c.bmi_sage_differential_ring (
                strders, strblks, strpars, notation)
        if bmi_c.bmi_sage_is_error (A):
            mesgerr = bmi_c.bmi_sage_mesgerr (A)
            bmi_c.bmi_balsa_clear_ALGEB (A)
            raise RuntimeError, mesgerr
# Build a dictionary to be used for future parsing operations
        self._names = {'diff' : derivative, 'Integer' : Integer, 'true' : True, 'false' : False}

        for p in flatten([parameters, blocks, derivations]):
            if isinstance(p, Expression) and p.operator() != None:
                self._names[bytes(p.operator())] = p.operator()
                for o in p.operands():
                    if not isinstance(o, Integer):
                        self._names[bytes(o)] = o
            else:
                self._names[bytes(p)] = p
# The dring field with its reference counter set to 1
        self.dring = A
        bmi_c.bmi_balsa_increment_nbref (self.dring)

    def __dealloc__ (self):
        """ The destructor """
        bmi_c.bmi_balsa_decrement_nbref (self.dring)
        bmi_c.bmi_balsa_clear_ALGEB (self.dring)

    def _translate_str(self, x):
        r"""
        Translates a string from Sage format to BMI's format.

        Used internally by the DifferentialAlgebra package.

        TESTS::

            sage: from sage.calculus.DifferentialAlgebra import DifferentialRing
            sage: R = DifferentialRing()
            sage: R._translate_str('diff(x(t), t)')
            'D[0](x)(t)'
            sage: R._translate_str('diff(x(t), t), diff(y(t), t)')
            'D[0](x)(t), D[0](y)(t)'
            sage: R._translate_str('diff(x(t), t, t)')
            'D[0,0](x)(t)'
            sage: R._translate_str('diff(f(x,y), x)')
            'D[0](f)(x,y)'
            sage: R._translate_str('diff(f(x,y), y)')
            'D[1](f)(x,y)'

        """

        def callable(m):
           func = m.group(1)
           funcvar1 = m.group(2)
           funcvars = [funcvar1] + re.split(r',\s*',m.group(3))[1:]
           diffvar1 = m.group(5)
           diffvars = [diffvar1] + re.split(r',\s*',m.group(6))[1:]

           return 'D[' + ','.join([str(funcvars.index(v)) for v in diffvars]) + "](" + func + ")(" + ','.join(funcvars) + ")"

        return re.sub(r'diff\((\w+)\((\w)+((,\s*\w+)*)\),\s*(\w+)((,\s*\w+)*)\)', callable, x)

    def _eval_sage (self, expr, extra_names = None):
        r"""
        Translates a string from BMI to Sage by evaluating the string
        as a Sage expression, using the ``_names`` dictionary to
        convert variables to Sage objects.

        Used internally by the DifferentialAlgebra package.

        TESTS::

            sage: from sage.calculus.DifferentialAlgebra import DifferentialRing
            sage: x,y = var('x,y')
            sage: f = function('f')
            sage: R = DifferentialRing(derivations = [x,y], blocks = [f])
            sage: R._eval_sage('x + y + diff(f(x), x)')
            x + y + diff(f(x), x)

        """
        if extra_names is None:
            return eval (preparse (expr), self._names)
        else:
            names = dict(self._names, **extra_names)
            return eval (preparse (expr), names)

    def __repr__ (self):
        """ The external representation """

        independents = self.indets(selection = 'independent')
        dependents = self.indets(selection = 'dependent')
        params = self.indets(selection = 'parameters')

        for p in params:
            dependents.remove(p)

        def convert_to_string(name, list):
            if len(list) == 0:
                return ''
            elif len(list) == 1:
                return name + ' ' + repr(list.pop())
            else:
                return name + 's ' + ', '.join([repr(i) for i in list])

        if len(params) > 0 and len(dependents) > 0 and len(independents) > 0:
            descr = convert_to_string(' with dependent', dependents) \
                    + convert_to_string(', independent', independents) \
                    + convert_to_string(', and parameter', params)
        elif len(dependents) > 0:
            descr = convert_to_string(' with dependent', dependents) \
                    + convert_to_string(' and independent', independents) \
                    + convert_to_string(' and parameter', params)
        elif len(independents) > 0:
            descr = convert_to_string(' with independent', independents) \
                    + convert_to_string(' and parameter', params)
        else:
            descr = convert_to_string(' with parameter', params)

        return 'Differential Ring over Rational Field' + descr

    def _latex_ (self):
        """ The LaTeX representation """

        independents = self.indets(selection = 'independent')
        dependents = self.indets(selection = 'dependent')

        descr = '\Bold{Q}'
        if len(independents) > 0:
            descr = descr + '[' + ','.join([latex(i) for i in independents]) + ']'
        if len(dependents) > 0:
            descr = descr + '\{' + ','.join([latex(d) for d in dependents]) + '\}'

        return descr

# coeffs
    def coeffs (
            self,
            equation,
            variable = None,
            BaseFieldExtension basefield = None,
            char* notation = BMI_IX_undefined):
        r"""
        The coefficients and monomials of a rational differential fraction

        INPUT:

        - ``equation``  -- a rational fraction
        - ``variable``  -- (default: None) a derivative
        - ``basefield`` -- (default: None) a base field
        - ``notation``  -- (default: 'undefined') a string

        OUTPUT:

        A sequence of two lists: the list ``C`` of the coefficients and
        the list ``M`` of the monomials. Both lists have the same length
        and the following relation holds:

        ``equation`` == add (C[i]*M[i] for i in range (len (C)))``

        The monomials ``M[i]`` are possibly rational fractions.

        If ``variable`` and ``basefield`` are both ``None``, then
        the coefficients ``C[i]`` are numerical.

        If ``variable`` is not ``None``, then ``basefield`` must be ``None``,
        the ``C[i]`` only depend on variables strictly less than
        ``variable``, and the ``M[i]``'s only depend on variables greater
        than or equal to ``variable``, w.r.t. the ranking.

        If ``basefield`` is not ``None``, then ``variable`` must be ``None``,
        the ``C[i]`` only depend on variables which belong to the field,
        and the ``M[i]`` only depend on variables which do not belong
        to the field.

        Restriction: the denominator of ``equation`` must depend on
        variables which are all, either strictly less than, or greater
        than ``variable``.

        NOTES:

        So far, there is a single notation available, because of
        limitations of function.

        EXAMPLES::

            sage: from sage.calculus.DifferentialAlgebra import DifferentialRing, RegularDifferentialChain, BaseFieldExtension
            sage: x,y = var('x,y')
            sage: f,u = function('f,u')

        The differential ring is equipped with an elimination ranking such
        that ``u >> f``::

            sage: R = DifferentialRing(derivations = [x,y], blocks = [u,f], parameters = [u(x)])

        Let us define a differential field ``F``::

            sage: eqns = [diff(f(x,y),x) == f(x,y)^2 + 1, diff(f(x,y),y) == 2*y*(f(x,y)^2 + 1)]
            sage: C = RegularDifferentialChain(eqns, R, pretend=false)
            sage: F = BaseFieldExtension(relations = C)

        Let us start with a differential polynomial::

            sage: eqn = x*diff(u(x),x)^2-(f(x,y)^2+1)*u(x) == 0

        The numerical coefficients and the corresponding monomials::

            sage: R.coeffs(eqn)
            ([1, -1, -1], [x*diff(u(x), x)^2, f(x, y)^2*u(x), u(x)])

        The coefficients of eqn w.r.t. ``diff(f(x,y),x)``. The variables ``x``
        and ``f(x,y)``, both less than ``diff(f(x,y),x)``, w.r.t. the ranking,
        arise in the coefficients::

            sage: R.coeffs( eqn, diff(f(x,y),x) )
            ([x, -f(x, y)^2 - 1], [diff(u(x), x)^2, u(x)])

        The coefficients of ``eqn``, viewed as a differential polynomials
        with coefficients in ``F``::

            sage: R.coeffs(eqn, basefield = F)
            ([x, -f(x, y)^2 - 1], [diff(u(x), x)^2, u(x)])

        The coefficients of ``eqn``, viewed as a differential polynomials
        with coefficients in `\QQ(x,y)`::

            sage: R.coeffs (eqn, basefield = BaseFieldExtension ())
            ([x, -1, -1], [diff(u(x), x)^2, f(x, y)^2*u(x), u(x)])

        The :meth:`coeffs` method also applies to rational fractions::

            sage: eqn = x*diff(u(x),x)^2-u(x)/diff(f(x,y),x)

        The numerical coefficients and the corresponding monomials,
        in the extended sense::

            sage: R.coeffs (eqn)
            ([1, -1], [x*diff(u(x), x)^2, u(x)/diff(f(x, y), x)])

        The coefficients of ``eqn``, viewed as a differential polynomials
        with coefficients in ``F``::

            sage: R.coeffs (eqn, basefield = F)
            ([x, -1/diff(f(x, y), x)], [diff(u(x), x)^2, u(x)])

        The coefficients of ``eqn``, viewed as a differential polynomials
        with coefficients in `\QQ(x,y)`::

            sage: R.coeffs (eqn, basefield = BaseFieldExtension ())
            ([x, -1], [diff(u(x), x)^2, u(x)/diff(f(x, y), x)])
        """
#
        cdef result
        cdef bytes streqns, strvar, strgens, strrels, mesgerr
        cdef bmi_c.ALGEB_string A
        if variable != None and basefield != None:
            raise RuntimeError, "the variable and the base field must not be both specified"
        streqns = bytes (equation)
        streqns = self._translate_str(streqns)
        if variable == None:
            if basefield == None:
                strvar = bytes (0)
                strgens = bytes ("")
                strrels = bytes ("")
            else:
                strvar = bytes ("")
                strgens = bytes (basefield.generators ())
                strrels = basefield.__relations ()
                strrels = self._translate_str(strrels)
        else:
            strvar = bytes (variable)
            strvar = self._translate_str(strvar)
            if basefield == None:
                strgens = bytes ("")
                strrels = bytes ("")
            else:
                raise RuntimeError, "the variable and the base field must not be both specified"
        A = bmi_c.bmi_sage_coeffs (
                    streqns, strvar, strgens, strrels, self.dring,
                    BMI_IX_undefined, notation, 0, 0)
        if bmi_c.bmi_sage_is_error (A):
            mesgerr = bmi_c.bmi_sage_mesgerr (A)
            bmi_c.bmi_balsa_clear_ALGEB (A)
            raise RuntimeError, mesgerr
        result = self._eval_sage(A.value)
        bmi_c.bmi_balsa_clear_ALGEB (A)
        return result

# equations
    def equations (
            self,
            equations,
            bool solved = False,
            selection = None,
            char* notation = BMI_IX_undefined):
        r"""
        A selection in a list of polynomials

        INPUT:

        - ``equations`` -- a polynomial or a list of polynomials
        - ``solved``    -- (default: False) a boolean
        - ``selection`` -- (default: None) a relational expression
        - ``notation``  -- (default: 'undefined') a string

        OUTPUT:

        The same list, but processed according to the other parameters.

        If ``solved`` is ``True``, the differential polynomials are
        displayed as equations, solved with respect to their
        leading rank.

        The optional parameter ``selection`` is a relational
        expression combining one of the keywords 'order', 'leader',
        'rank' on the one side, with a value on the other side,
        using one of the relational operators. Only the equations
        which satisfy the condition are selected.

        - If the keyword is ``order``, then the value must be a
          nonnegative integer.

        - If the keyword is ``rank``, then the value must be a
          rank, including the special ranks 0 and 1.

        - If the keyword is ``leader``, and the relational
          operator is ``>=``, ``>``, ``<=`` or ``<``, then the value
          must be a derivative.

        - If the keyword is ``leader``, and the relational
          operator is ``==`` or ``!=``, then the value may also have
          the form: ``modifier`` (derivative), where ``modifier`` is
          one of the keywords ``derivative`` or ``proper``.

        NOTES:

        So far, there is a single notation available, because of
        limitations of function.

        EXAMPLES::

            sage: from sage.calculus.DifferentialAlgebra import DifferentialRing
            sage: x,y = var('x,y')
            sage: w,z,a = function ('w,z,a')
            sage: R = DifferentialRing (derivations = [x,y], blocks = [[w,a],z], parameters = [a(y)])
            sage: L = [0, 18, a(y)^2, diff(a(y),y), w(x,y) - z(x,y)]
            sage: leader,order,rank = var('leader,order,rank')
            sage: derivative,proper = function('derivative,proper')

        Just check that the elements of ``L`` are polynomials of ``R``::

            sage: R.equations (L)
            [0, 18, a(y)^2, diff(a(y), y), w(x, y) - z(x, y)]

        The only polynomial with rank 0 is 0::

            sage: R.equations (L, selection = rank == 0)
            [0]

        The polynomials with rank 1 are the nonzero numbers::

            sage: R.equations (L, selection = rank == 1)
            [18]

        The polynomials which are not numbers, solved with respect
        to their leaders::

            sage: R.equations (L, solved = True, selection = rank > 1)
            [a(y)^2 == 0, diff(a(y), y) == 0, w(x, y) == z(x, y)]

        The polynomials whose leaders are derivatives of `a(y)`::

            sage: R.equations (L, selection = leader == derivative (a(y)))
            [a(y)^2, diff(a(y), y)]

        The polynomials whose leaders are derivatives of `w(x,y)`, solved
        with respect to their leaders::

            sage: R.equations (L, solved = True, selection = leader == derivative (w(x,y)))
            [w(x, y) == z(x, y)]

        Solving a single equation w.r.t. its leading derivative::

            sage: R.equations (w(x,y) - a(y), solved=True)
            w(x, y) == a(y)
        """
#
        cdef result
        cdef bytes streqns, strsel, mesgerr
        cdef bmi_c.ALGEB_string A
        if isinstance (equations, list):
            streqns = bytes (equations)
        else:
            streqns = bytes ([equations])
        streqns = self._translate_str(streqns)
        if selection == None:
            strsel = bytes ('rank >= 0')
        else:
            strsel = bytes (selection)
        A = bmi_c.bmi_sage_equations_with_criterion_DR (
                    streqns, self.dring, int (False), int (solved), strsel,
                    BMI_IX_undefined, notation, 0, 0)
        if bmi_c.bmi_sage_is_error (A):
            mesgerr = bmi_c.bmi_sage_mesgerr (A)
            bmi_c.bmi_balsa_clear_ALGEB (A)
            raise RuntimeError, mesgerr
        result = self._eval_sage(A.value)
        if not isinstance (equations, list):
            result = result [0]
        bmi_c.bmi_balsa_clear_ALGEB (A)
        return result

# factor_derivative
    def factor_derivative (
            self,
            derivative,
            char* notation = BMI_IX_undefined):
        r"""
        Splits a derivative into a derivation operator and a dependent variable

        INPUT:

        - ``derivative`` -- a derivative
        - ``notation``   -- (default: 'undefined') a string

        OUTPUT:

        A pair (``theta``, ``symb``) such that ``derivative`` is
        the derivative of ``symb`` w.r.t. ``theta``. The derivation
        operator ``theta`` is expressed as a monomial on the
        independent variables.

        NOTES:

        So far, there is a single notation available, because of
        limitations of function.

        EXAMPLES::

            sage: from sage.calculus.DifferentialAlgebra import DifferentialRing
            sage: a,x,y = var ('a,x,y')
            sage: u,v = function ('u,v')
            sage: R = DifferentialRing (derivations = [x,y], blocks = [u,v,a], parameters = [a,v(y)])

        A basic example::

            sage: R.factor_derivative (diff (u(x,y),x))
            (x, u(x, y))

        A parameter is considered as a dependent variable which does not
        depend on any independent variable. The monomial ``theta`` is `1`,
        meaning that `a` is not differentiated::

            sage: R.factor_derivative (a)
            (1, a)

        The following examples illustrates the relationship between
        :meth:`factor_derivative` and :meth:`differentiate`::

            sage: derv = diff (u(x,y),x,x,y)
            sage: derv
            diff(u(x, y), x, x, y)
            sage: theta, symb = R.factor_derivative (derv)
            sage: theta
            x^2*y
            sage: symb
            u(x, y)
            sage: R.differentiate (symb, theta)
            diff(u(x, y), x, x, y)
        """
        cdef result
        cdef bytes strder, mesgerr
        cdef bmi_c.ALGEB_string A
        strder = bytes (derivative)
        strder = self._translate_str(strder)
        A = bmi_c.bmi_sage_factor_derivative (
                    strder, self.dring,
                    BMI_IX_undefined, notation, 0, 0)
        if bmi_c.bmi_sage_is_error (A):
            mesgerr = bmi_c.bmi_sage_mesgerr (A)
            bmi_c.bmi_balsa_clear_ALGEB (A)
            raise RuntimeError, mesgerr
        result = self._eval_sage(A.value)
        bmi_c.bmi_balsa_clear_ALGEB (A)
        return result

# indets
    def indets (
                self,
                equations = None,
                char* selection = BMI_IX_derivs,
                derivation = None,
                char* notation = BMI_IX_undefined):
        r"""
        The indeterminates of a DifferentialRing or a rational fraction

        INPUT:

        - ``equations``  -- (default: None) a rational fraction
                             or a list of rational fractions
        - ``selection``  -- (default: 'derivatives') a string
        - ``derivation`` -- (default: None) a variable
        - ``notation``   -- (default: 'undefined') a string

        OUTPUT:

        A list of variable and function, selected according to
        ``selection``. See the examples below.

        NOTES:

        So far, there is a single notation available, because of
        limitations of function.

        EXAMPLES::

            sage: from sage.calculus.DifferentialAlgebra import DifferentialRing
            sage: x,y,z,t = var ('x,y,z,t')
            sage: u,v,w = function ('u,v,w')
            sage: R = DifferentialRing (derivations = [x,y,t], blocks = [u,[v,w],z], parameters = [u(t,y,x),z,w(y)])

        With no selection, or with ``derivatives`` or ``dependent``,
        return the list of the dependent variables (or differential
        indeterminates) of `R`::

            sage: R.indets()
            [u(t, y, x), v(x, y, t), w(y), z]
            sage: R.indets (selection = 'dependent')
            [u(t, y, x), v(x, y, t), w(y), z]
            sage: R.indets (selection = 'derivatives')
            [u(t, y, x), v(x, y, t), w(y), z]

        With ``derivations`` or ``independent``, return the list of
        the derivations (or independent variables) of `R`::

            sage: from sage.calculus.DifferentialAlgebra import DifferentialRing
            sage: R.indets (selection = 'derivations')
            [x, y, t]
            sage: R.indets (selection = 'independent')
            [x, y, t]

        With ``parameters``, return the list of the parameters of `R`::

            sage: R.indets (selection = 'parameters')
            [u(t, y, x), w(y), z]

        With ``all``, return the list of the independent and the dependent
        variables of `R`::

            sage: R.indets (selection = 'all')
            [x, y, t, u(t, y, x), v(x, y, t), w(y), z]

        With ``constants``, possibly followed by a derivation, return
        the list of the dependent variables the derivatives (w.r.t.
        the derivation, if any) of which, are zero::

            sage: R.indets (selection = 'constants')
            [z]
            sage: R.indets (selection = 'constants', derivation = y)
            [z]
            sage: R.indets (selection = 'constants', derivation = t)
            [w(y), z]

        The indets of a rational fraction::

            sage: eq = u(t,y,x) + 1/w(y)
            sage: R.indets (eq)
            [u(t, y, x), w(y)]
            sage: R.indets (eq, selection = 'constants', derivation = x)
            [w(y)]
        """
        cdef list result
        cdef bytes streqns, strder, mesgerr
        cdef bmi_c.ALGEB_string L
        if equations == None:
            streqns = bytes ("")
        elif isinstance (equations, list):
            streqns = bytes (equations)
        else:
            streqns = bytes ([equations])
        streqns = self._translate_str(streqns)
        if derivation == None:
            strder = bytes (0)
        else:
            strder = bytes (derivation)
        L = bmi_c.bmi_sage_indets (
                streqns, self.dring, selection, strder,
                BMI_IX_undefined, notation, 0, 0)
        if bmi_c.bmi_sage_is_error (L):
            mesgerr = bmi_c.bmi_sage_mesgerr (L)
            bmi_c.bmi_balsa_clear_ALGEB (L)
            raise RuntimeError, mesgerr
        result = self._eval_sage(L.value)
        bmi_c.bmi_balsa_clear_ALGEB (L)
        return result

# is_constant
    def is_constant (
            self,
            ratfrac,
            derivation = None):
        r"""
        Test if a rational fraction is constant.

        INPUT:

        - ``ratfrac``    -- a rational fraction or a list of rational
                            fractions
        - ``derivation`` -- (default: None) an independent variable

        OUTPUT:

        ``True``, if the rational fraction is a constant w.r.t.
        ``derivation``, else ``False``. If ``derivation`` is omitted,
        true if the rational fraction is a constant w.r.t.
        all the derivations of the ring.

        EXAMPLES::

            sage: from sage.calculus.DifferentialAlgebra import DifferentialRing
            sage: a,x,y = var ('a,x,y')
            sage: u,v = function ('u,v')
            sage: R = DifferentialRing (derivations = [x,y], blocks = [u,v,a], parameters = [a,v(y)])

            sage: R.is_constant (3)
            True

            sage: R.is_constant( [1/a, v(y), diff(u(x,y),x)/(a+1)] )
            [True, False, False]

            sage: R.is_constant( [1/a, v(y), diff(u(x,y),x)/(a+1)], x )
            [True, True, False]
        """
        cdef result
        cdef bytes streqns, strder, mesgerr
        cdef bmi_c.ALGEB_string A
        if isinstance (ratfrac, list):
            streqns = bytes (ratfrac)
        else:
            streqns = bytes ([ratfrac])
        streqns = self._translate_str(streqns)
        if derivation == None:
            strder = bytes (0)
        else:
            strder = bytes (derivation)
        A = bmi_c.bmi_sage_is_constant (
                    streqns, strder, self.dring,
                    BMI_IX_undefined, BMI_IX_undefined, 0, 0)
        if bmi_c.bmi_sage_is_error (A):
            mesgerr = bmi_c.bmi_sage_mesgerr (A)
            bmi_c.bmi_balsa_clear_ALGEB (A)
            raise RuntimeError, mesgerr
        result = self._eval_sage(A.value)
        if not isinstance (ratfrac, list):
            result = result [0]
        bmi_c.bmi_balsa_clear_ALGEB (A)
        return result

# leading_derivative
    def leading_derivative (
            self,
            equations,
            char* notation = BMI_IX_undefined):
        r"""
        The leading derivative of a rational fraction.

        INPUT:

        - ``equations`` -- a rational fraction or a list of
                           rational fractions
        - ``notation``  -- (default: 'undefined') a string

        OUTPUT:

        The list of the leading derivatives of the equations.

        The leading derivative of a rational fraction F is
        defined as the highest derivative u such that the
        partial derivative of F w.r.t. u is nonzero.

        This definition is extended to rational fractions
        which only depend on independent variables.

        The leading derivative of constant rational fractions
        is not defined.

        NOTES:

        So far, there is a single notation available, because of
        limitations of function.

        EXAMPLES::

            sage: from sage.calculus.DifferentialAlgebra import DifferentialRing
            sage: t = var ('t')
            sage: u,v = function ('u,v')
            sage: R = DifferentialRing (derivations = [t], blocks = [u,v])
            sage: L = [t, u(t)**2 - u(t), diff(v(t),t), 1/v(t)]
            sage: R.leading_derivative (L)
            [t, u(t), diff(v(t), t), v(t)]
            sage: R.leading_derivative (1/u(t)**2)
            u(t)
            sage: R.leading_derivative (u(t)**2)
            u(t)
        """
#
        cdef result
        cdef bytes streqns, mesgerr
        cdef bmi_c.ALGEB_string A
        if isinstance (equations, list):
            streqns = bytes (equations)
        else:
            streqns = bytes ([equations])
        streqns = self._translate_str(streqns)
        A = bmi_c.bmi_sage_leading_derivative (
                    streqns, self.dring, int (False),
                    BMI_IX_undefined, notation, 0, 0)
        if bmi_c.bmi_sage_is_error (A):
            mesgerr = bmi_c.bmi_sage_mesgerr (A)
            bmi_c.bmi_balsa_clear_ALGEB (A)
            raise RuntimeError, mesgerr
        result = self._eval_sage(A.value)
        if not isinstance (equations, list):
            result = result [0]
        bmi_c.bmi_balsa_clear_ALGEB (A)
        return result

# leading_rank
    def leading_rank (
            self,
            equations,
            bool listform = False,
            char* notation = BMI_IX_undefined):
        r"""
        The leading rank of a rational fraction.

        INPUT:

        - ``equations`` -- a rational fraction or a list of
                           rational fractions
        - ``listform``  -- (default: False) a bool
        - ``notation``  -- (default: 'undefined') a string

        OUTPUT:

        The list of the leading ranks of the equations.

        If the optional parameter ``listform`` is ``True``, then
        the ranks are returned as lists rather than monomials.
        This is feature is relevant when the rank is a derivative
        raised at the 0-th power.

        The leading rank of a rational fraction `F = P/Q` is defined
        as the leading derivative `u` of `F`, raised at the power
        `deg (P, u) - deg (Q, u)`.

        Nonzero constant rational fractions have leading rank 1.

        The zero rational fraction has leading rank 0.

        NOTES:

        So far, there is a single notation available, because of
        limitations of function.

        EXAMPLES::

            sage: from sage.calculus.DifferentialAlgebra import DifferentialRing
            sage: t = var ('t')
            sage: u,v = function ('u,v')
            sage: R = DifferentialRing (derivations = [t], blocks = [u,v])
            sage: L = [t, u(t)**2 - u(t), diff(v(t),t), 1/v(t)]

        The list of leading ranks of the elements of `L`, in monomial form::

            sage: R.leading_rank (L)
            [t, u(t)^2, diff(v(t), t), 1/v(t)]

        The same list, in list form::

            sage: R.leading_rank (L, listform = True)
            [[t, 1], [u(t), 2], [diff(v(t), t), 1], [v(t), -1]]

        The case of a single rational fraction::

            sage: R.leading_rank (1/u(t)**2)
            u(t)^(-2)
            sage: R.leading_rank (1/u(t)**2, listform = True)
            [u(t), -2]
            sage: R.leading_rank (u(t)**2)
            u(t)^2
            sage: R.leading_rank (u(t)**2, listform = True)
            [u(t), 2]

        The case of a rank, being a derivative raised at the 0-th power::

            sage: R.leading_rank ((u(t)-1)/(u(t)+2))
            1
            sage: R.leading_rank ((u(t)-1)/(u(t)+2), listform = True)
            [u(t), 0]
        """
        cdef result
        cdef bytes streqns, mesgerr
        cdef bmi_c.ALGEB_string A
        if isinstance (equations, list):
            streqns = bytes (equations)
        else:
            streqns = bytes ([equations])
        streqns = self._translate_str(streqns)
        A = bmi_c.bmi_sage_leading_rank (
                    streqns, self.dring, int (False), int (listform),
                    BMI_IX_undefined, notation, 0, 0)
        if bmi_c.bmi_sage_is_error (A):
            mesgerr = bmi_c.bmi_sage_mesgerr (A)
            bmi_c.bmi_balsa_clear_ALGEB (A)
            raise RuntimeError, mesgerr
        result = self._eval_sage(A.value)
        if not isinstance (equations, list):
            result = result [0]
        bmi_c.bmi_balsa_clear_ALGEB (A)
        return result

    def initial (
            self,
            equations,
            char* notation = BMI_IX_undefined):
        r"""
        The initial of a rational fraction.

        INPUT:

        - ``equations`` -- a rational fraction or a list of
                           rational fractions
        - ``notation``  -- (default: 'undefined') a string

        OUTPUT:

        The list of the initials of the equations, i.e. the
        leading coefficients of the equations w.r.t. their
        leading derivatives.

        The initial of a rational fraction `F = P/Q`, with
        leading derivative `u`, is itself a rational fraction:
        the leading coefficient of `P` w.r.t. `u`, divided by
        the leading coefficient of `Q` w.r.t. `u`, recalling
        that the leading coefficient w.r.t. `u`, of a polynomial
        `R` which does not depend on `u`, is defined as `R` itself.

        NOTES:

        The above definition is not yet fully implemented.

        So far, there is a single notation available, because of
        limitations of function.

        EXAMPLES::

            sage: from sage.calculus.DifferentialAlgebra import DifferentialRing
            sage: t = var ('t')
            sage: u,v = function ('u,v')
            sage: R = DifferentialRing (derivations = [t], blocks = [u,v])
            sage: L = [t, v(t)*u(t)**2 - u(t) + diff (v(t),t) + 1, t*u(t)/diff(v(t),t)]
            sage: R.initial (L)
            [1, v(t), t/diff(v(t), t)]
            sage: R.initial (u(t)**2 + 1/v(t))
            1
        """
        cdef result
        cdef bytes streqns, strvar, mesgerr
        cdef bmi_c.ALGEB_string A
        if isinstance (equations, list):
            streqns = bytes (equations)
        else:
            streqns = bytes ([equations])
        strvar = bytes (0)
        streqns = self._translate_str(streqns)
        A = bmi_c.bmi_sage_leading_coefficient (
                    streqns, self.dring, int (False), strvar,
                    BMI_IX_undefined, notation, 0, 0)
        if bmi_c.bmi_sage_is_error (A):
            mesgerr = bmi_c.bmi_sage_mesgerr (A)
            bmi_c.bmi_balsa_clear_ALGEB (A)
            raise RuntimeError, mesgerr
        result = self._eval_sage(A.value)
        if not isinstance (equations, list):
            result = result [0]
        bmi_c.bmi_balsa_clear_ALGEB (A)
        return result

# leading_coefficient
    def leading_coefficient (
            self,
            equations,
            variable = None,
            char* notation = BMI_IX_undefined):
        r"""
        The leading coefficient of a rational fraction.

        INPUT:

        - ``equations`` -- a rational fraction or a list of
                           rational fractions
        - ``variable``  -- (default: None) a variable
        - ``notation``  -- (default: 'undefined') a string

        OUTPUT:

        The list of the leading coefficients of the equations,
        w.r.t. ``variable``. If ``variable`` is ``None``, the
        list of the initials of the equations is returned.

        The leading coefficient of a rational fraction `F = P/Q`,
        w.r.t. some derivative `u`, is a rational fraction:
        it is the leading coefficient of `P` w.r.t. `u`, divided
        by the leading coefficient of `Q` w.r.t. `u`,  recalling
        that the leading coefficient w.r.t. `u`, of a polynomial
        `R` which does not depend on `u`, is defined as `R` itself.

        NOTES:

        The above definition is not yet fully implemented.

        So far, there is a single notation available, because of
        limitations of function.

        EXAMPLES::

            sage: from sage.calculus.DifferentialAlgebra import DifferentialRing
            sage: t = var ('t')
            sage: u,v = function ('u,v')
            sage: R = DifferentialRing (derivations = [t], blocks = [u,v])
            sage: L = [t, v(t)*u(t)**2 - u(t) + diff (v(t),t) + 1, t*u(t)/diff(v(t),t)]

        The initials of the elements of ``L``::

            sage: R.leading_coefficient (L)
            [1, v(t), t/diff(v(t), t)]

        The initial of a single rational fraction::

            sage: R.leading_coefficient (u(t)**2 + diff(v(t),t)/v(t))
            1

        The leading coefficient of the elements of ``L``, w.r.t. t.::

            sage: R.leading_coefficient (L, t)
            [1, u(t)^2*v(t) - u(t) + diff(v(t), t) + 1, u(t)/diff(v(t), t)]

        The leading coefficient of a single rational fraction,
        w.r.t. a derivative.::

            sage: R.leading_coefficient (u(t)**2 + diff(v(t),t)/v(t), diff(v(t),t))
            1/v(t)
        """
        cdef result
        cdef bytes streqns, strvar, mesgerr
        cdef bmi_c.ALGEB_string A
        if isinstance (equations, list):
            streqns = bytes (equations)
        else:
            streqns = bytes ([equations])
        streqns = self._translate_str(streqns)
        if variable == None:
            strvar = bytes (0)
        else:
            strvar = bytes (variable)
        strvar = self._translate_str(strvar)
        A = bmi_c.bmi_sage_leading_coefficient (
                    streqns, self.dring, int (False), strvar,
                    BMI_IX_undefined, notation, 0, 0)
        if bmi_c.bmi_sage_is_error (A):
            mesgerr = bmi_c.bmi_sage_mesgerr (A)
            bmi_c.bmi_balsa_clear_ALGEB (A)
            raise RuntimeError, mesgerr
        result = self._eval_sage(A.value)
        if not isinstance (equations, list):
            result = result [0]
        bmi_c.bmi_balsa_clear_ALGEB (A)
        return result

# tail
    def tail (
            self,
            equations,
            variable = None,
            char* notation = BMI_IX_undefined):
        r"""
        The tail of a rational fraction

        INPUT:

        - ``equations`` -- a rational fraction or a list of
                           rational fractions
        - ``variable``  -- (default: None) a variable
        - ``notation``  -- (default: 'undefined') a string

        OUTPUT:

        The list of tails of the equations, w.r.t. ``variable``.
        If ``variable`` is omitted, the list of tails w.r.t. the
        leading derivatives of the equations is returned.

        The tail of a rational fraction ``F`` is defined as
        the rational fraction ``F - i(F) * rg(F)`` where ``i(F)``
        and ``rg(F)`` denote the initial and the leading rank of ``F``.

        More generally, the tail, w.r.t. some derivative `u`, of
        a rational fraction `F = P/Q` is defined as the rational
        fraction ``F - lc(F,u) * rg(F,u)``, where ``lc(F,u)`` is defined
        as the leading coefficient of ``F`` w.r.t. u and ``rg(F,u)`` is
        the monomial obtained by raising `u` at the degree
        `deg (P,u) - deg (Q,u)`

        NOTES:

        The above definition is not yet fully implemented.

        So far, there is a single notation available, because of
        limitations of function.

        EXAMPLES::

            sage: from sage.calculus.DifferentialAlgebra import DifferentialRing
            sage: t = var ('t')
            sage: u,v = function ('u,v')
            sage: R = DifferentialRing (derivations = [t], blocks = [u,v])
            sage: L = [t + 3, v(t)*u(t)**2 - u(t) + diff (v(t),t) + 1, t*u(t)/diff(v(t),t)]

        The list of tails of the elements of `L`, w.r.t. their leading
        derivatives::

            sage: R.tail (L)
            [3, -u(t) + diff(v(t), t) + 1, 0]

        The tail of a single rational fraction, w.r.t. its leading
        derivative::

            sage: eq = u(t)**2 + diff(v(t),t)/v(t)
            sage: R.tail (eq)
            diff(v(t), t)/v(t)

        Check the definition of tails::

            sage: R.initial (eq) * R.leading_rank (eq) + R.tail (eq) - eq
            0

        The list of tails of the elements of `L`, w.r.t. `t`::

            sage: R.tail (L, t)
            [3, 0, 0]

        The tail of a single rational fraction, w.r.t. a derivative::

            sage: eq = u(t)**2 + diff(v(t),t)/v(t)
            sage: derv = diff (v(t), t)
            sage: R.tail (eq, derv)
            u(t)^2

        Check the definition of tails w.r.t. some derivative::

            sage: R.leading_coefficient (eq, derv) * derv + R.tail (eq, derv) - eq
            0
        """
        cdef result
        cdef bytes streqns, strvar, mesgerr
        cdef bmi_c.ALGEB_string A
        if isinstance (equations, list):
            streqns = bytes (equations)
        else:
            streqns = bytes ([equations])
        streqns = self._translate_str(streqns)
        if variable == None:
            strvar = bytes (0)
        else:
            strvar = bytes (variable)
            strvar = self._translate_str (strvar)
        A = bmi_c.bmi_sage_tail (
                    streqns, self.dring, int (False), strvar,
                    BMI_IX_undefined, notation, 0, 0)
        if bmi_c.bmi_sage_is_error (A):
            mesgerr = bmi_c.bmi_sage_mesgerr (A)
            bmi_c.bmi_balsa_clear_ALGEB (A)
            raise RuntimeError, mesgerr
        result = self._eval_sage(A.value)
        if not isinstance (equations, list):
            result = result [0]
        bmi_c.bmi_balsa_clear_ALGEB (A)
        return result

    def separant (
            self,
            equations,
            variable = None,
            char* notation = BMI_IX_undefined):
        r"""
        The separant of a rational fraction.

        INPUT:

        - ``equations`` -- a rational fraction or a list of
                           rational fractions
        - ``variable``  -- (default: None) a variable
        - ``notation``  -- (default: 'undefined') a string

        OUTPUT:

        The list of separants of the equations, w.r.t. ``variable``.
        If ``variable`` is omitted, the list of separants of the
        equations, w.r.t. their leading derivatives, is returned.

        The separant of a rational fraction `F` w.r.t. some derivative
        `u`, is defined as the partial derivative of `F` w.r.t. u.

        NOTES:

        So far, there is a single notation available, because of
        limitations of function.

        EXAMPLES::

            sage: from sage.calculus.DifferentialAlgebra import DifferentialRing
            sage: t = var ('t')
            sage: u,v = function ('u,v')
            sage: R = DifferentialRing (derivations = [t], blocks = [u,v])
            sage: L = [t + 3, v(t)*u(t)**2 - u(t) + diff (v(t),t) + 1, t*u(t)/diff(v(t),t)]

        The list of separants of the elements of `L`, w.r.t. their
        leading derivatives::

            sage: R.separant (L)
            [1, 2*u(t)*v(t) - 1, t/diff(v(t), t)]

        The separant of a single rational fraction, w.r.t. its
        leading derivative::

            sage: R.separant (u(t)**2 + diff(v(t),t)/v(t))
            2*u(t)

        The list of separants of the elements of `L`, w.r.t. `t`::

            sage: R.separant (L, t)
            [1, 0, u(t)/diff(v(t), t)]

        The separant of a single rational fraction, w.r.t. a derivative::

            sage: R.separant (u(t)**2 + diff(v(t),t)/v(t), diff(v(t),t))
            1/v(t)
        """
        cdef result
        cdef bytes streqns, strvar, mesgerr
        cdef bmi_c.ALGEB_string A
        if isinstance (equations, list):
            streqns = bytes (equations)
        else:
            streqns = bytes ([equations])
        streqns = self._translate_str(streqns)
        if variable == None:
            strvar = bytes (0)
        else:
            strvar = bytes (variable)
            strvar = self._translate_str(strvar)
        A = bmi_c.bmi_sage_separant (
                    streqns, self.dring, int (False), strvar,
                    BMI_IX_undefined, notation, 0, 0)
        if bmi_c.bmi_sage_is_error (A):
            mesgerr = bmi_c.bmi_sage_mesgerr (A)
            bmi_c.bmi_balsa_clear_ALGEB (A)
            raise RuntimeError, mesgerr
        result = self._eval_sage(A.value)
        if not isinstance (equations, list):
            result = result [0]
        bmi_c.bmi_balsa_clear_ALGEB (A)
        return result

# integrate
    def integrate (
            self,
            equations,
            variable,
            bool iterated = True,
            char* notation = BMI_IX_undefined):
        r"""
        Integrate a rational differential fraction.

        INPUT:

        - ``equations`` -- a rational fraction or a list of
                           rational fractions
        - ``variable``  -- a variable
        - ``iterated``  -- (default: True) a bool
        - ``notation``  -- (default: 'undefined') a string

        OUTPUT:

        A list of rational fractions, representing a decomposition
        of the equation or, a list of lists of rational fractions
        if ``equations`` is a list.

        The parameter ``variable`` must be an independent variable.

        Assume ``equations`` is a single rational fraction.
        Denote ``R`` the differential ring, ``F`` the rational
        fraction, ``t`` the independent variable and ``L`` the
        returned list. Then

        ``F == add( R.differentiate( L[i], t^i ) for i in range(len (L)) )``

        Moreover, if ``F`` is the derivative of some other
        rational fraction, w.r.t. ``t``, then ``L[0]`` is zero.

        The returned list is ranking dependent. Orderly rankings
        usually give more interesting results.

        NOTES:

        So far, there is a single notation available, because of
        limitations of function.

        EXAMPLES::

            sage: from sage.calculus.DifferentialAlgebra import DifferentialRing
            sage: t = var ('t')
            sage: u,v = function ('u,v')
            sage: R = DifferentialRing (derivations = [t], blocks = [[u,v]])

        Observe that the ranking is orderly::

            sage: eq = diff (u(t),t)^2 + diff (v(t)^2,t,t) - diff (u(t)*v(t) + t*v(t), t)
            sage: L = R.integrate (eq, t)
            sage: L
            [diff(u(t), t)^2, -t*v(t) - u(t)*v(t), v(t)^2]

        The two expressions are equal, for their difference is `0`::

            sage: eq - add (R.differentiate (L[i], t^i) for i in range (len (L)))
            0

        The case of a rational fraction which is the derivative of
        some other fraction::

            sage: eq = (v(t)^2 + t)/(3*diff(v(t),t))
            sage: L = R.integrate (diff (eq, t), t)
            sage: L
            [0, 1/3*(v(t)^2 + t)/diff(v(t), t)]
            sage: L[1] - eq
            0
        """
        cdef list result
        cdef bytes streqns, strvar, mesgerr
        cdef bmi_c.ALGEB_string A

        if isinstance (equations, list):
            streqns = bytes (equations)
        else:
            streqns = bytes ([equations])
        strvar = bytes (variable)
        streqns = self._translate_str(streqns)
        A = bmi_c.bmi_sage_integrate (
                    self.dring, streqns, strvar, int (iterated),
                    BMI_IX_undefined, notation, 0, 0)
        if bmi_c.bmi_sage_is_error (A):
            mesgerr = bmi_c.bmi_sage_mesgerr (A)
            bmi_c.bmi_balsa_clear_ALGEB (A)
            raise RuntimeError, mesgerr
        result = self._eval_sage(A.value)
        if not isinstance (equations, list):
            result = result [0]
        bmi_c.bmi_balsa_clear_ALGEB (A)
        return result

    def differentiate (
            self,
            equations,
            *args,
            char* notation = BMI_IX_undefined):
        r"""
        Differentiate a rational differential fraction.

        INPUT:

        - ``equations`` -- a rational fraction or a list of
                           rational fractions
        - ``args``      -- a possibly empty sequence of monomials
                           involving independent variables only
        - ``notation``  -- (default: 'undefined') a string

        OUTPUT:

        The list of derivatives of the equations, w.r.t. the
        product of the monomials given in ``args``. The monomials
        must depend on independent variables only.

        NOTES:

        So far, there is a single notation available, because of
        limitations of function.

        EXAMPLES::

            sage: from sage.calculus.DifferentialAlgebra import DifferentialRing
            sage: t = var ('t')
            sage: u,v = function ('u,v')
            sage: R = DifferentialRing (derivations = [t], blocks = [u,v])
            sage: eq = u(t)^2

        Do not differentiate at all::

            sage: R.differentiate (eq)
            u(t)^2
            sage: R.differentiate (eq, 1)
            u(t)^2

        Similar behaviour as ``diff``::

            sage: R.differentiate (eq, t)
            2*u(t)*diff(u(t), t)
            sage: R.differentiate (eq, t, t)
            2*diff(u(t), t)^2 + 2*u(t)*diff(u(t), t, t)

        Another way to differentiate twice w.r.t. `t`::

            sage: R.differentiate (eq, t^2)
            2*diff(u(t), t)^2 + 2*u(t)*diff(u(t), t, t)

        Differentiate three times w.r.t. `t`::

            sage: R.differentiate (eq, t^2, t)
            6*diff(u(t), t)*diff(u(t), t, t) + 2*u(t)*diff(u(t), t, t, t)
        """
        cdef result
        cdef bytes streqns, strders, mesgerr
        cdef bmi_c.ALGEB_string A

        if isinstance (equations, list):
            streqns = bytes (equations)
        else:
            streqns = bytes ([equations])
        streqns = self._translate_str(streqns)
        strders = bytes ([arg for arg in args])
        A = bmi_c.bmi_sage_differentiate (
                    self.dring, streqns, int (False), strders,
                    BMI_IX_undefined, notation, 0, 0)
        if bmi_c.bmi_sage_is_error (A):
            mesgerr = bmi_c.bmi_sage_mesgerr (A)
            bmi_c.bmi_balsa_clear_ALGEB (A)
            raise RuntimeError, mesgerr
        result = self._eval_sage(A.value)
        if not isinstance (equations, list):
            result = result [0]
        bmi_c.bmi_balsa_clear_ALGEB (A)
        return result

# differential_prem
    def differential_prem (
             self,
             polynomial,
             list redset,
             char* mode = BMI_IX_fully,
             char* notation = BMI_IX_undefined):
        r"""
        Ritt's reduction algorithm by a set of polynomials.

        INPUT:

        - ``polynomial`` -- a polynomial
        - ``redset``     -- a list of polynomials with integer coefficients
        - ``mode``       -- (default: 'full') a string
        - ``notation``   -- (default: 'undefined') a string

        OUTPUT:

        A pair `(h, r)` such that `h*p = r` modulo the differential
        ideal generated by ``redset`` (denoting `p` for ``polynomial``).

        The optional argument ``mode = 'full'``, ``'partial'`` or ``'algebraic'``.
        If ``'full'``, the remainder `r` is fully reduced w.r.t. ``redset``
        If ``'partial'``, the remainder `r` is partially reduced w.r.t.
        ``redset``. If ``'algebraic'``, the remainder is algebraically
        reduced w.r.t. ``redset``.

        In all cases, `h` is a product of powers of the initials and
        separants involved in the reduction process.

        NOTES:

        So far, there is a single notation available, because of
        limitations of function.

        EXAMPLES::

            sage: from sage.calculus.DifferentialAlgebra import DifferentialRing, RegularDifferentialChain, BaseFieldExtension
            sage: x,y = var ('x,y')
            sage: u,v = function ('u,v')
            sage: R = DifferentialRing (derivations = [x,y], blocks = [[v,u]])
            sage: eq1 = diff(u(x,y),x)^2 - 4*u(x,y)
            sage: eq2 = diff(u(x,y),x,y)*diff(v(x,y),y) - u(x,y) + 1
            sage: eq3 = diff(v(x,y),x,x) - diff(u(x,y),x)
            sage: redset = [eq1, eq2, eq3]

        Here is the reduction set::

            sage: R.equations (redset, solved = true)
            [diff(u(x, y), x)^2 == 4*u(x, y),
             diff(u(x, y), x, y) == (u(x, y) - 1)/diff(v(x, y), y),
             diff(v(x, y), x, x) == diff(u(x, y), x)]
            sage: ideal = R.RosenfeldGroebner (redset)
            sage: ideal
            [Regular differential chain (diff(u(x, y), y)^2 - 2*u(x, y),
                 diff(u(x, y), x)^2 - 4*u(x, y),
                -u(x, y)*diff(u(x, y), x)*diff(u(x, y), y) + diff(u(x, y), x)*diff(u(x, y), y)
                     + 4*u(x, y)*diff(v(x, y), y),
                -diff(u(x, y), x) + diff(v(x, y), x, x))]

        Here is a description of the radical differential ideal
        generated by the reduction set::

            sage: C = ideal [0]

        The first polynomial is fully reduced w.r.t. ``redset``::

            sage: poly = (1/7)*u(x,y)
            sage: h, r = R.differential_prem (poly, redset)
            sage: h, r
            (1, 1/7*u(x, y))
            sage: C.normal_form (h * poly - r)
            0

        The next polynomial is not fully reduced w.r.t ``redset``::

            sage: poly = (1/7)*diff (v(x,y),x,y) + diff (u(x,y),x,x)
            sage: h, r = R.differential_prem (poly, redset)
            sage: h, r
            (2*diff(u(x, y), x), 2/7*(diff(v(x, y), x, y) + 14)*diff(u(x, y), x))
            sage: C.normal_form (h * poly - r)
            0

        It is not partially reduced either::

            sage: h, r = R.differential_prem (poly, redset, mode = 'partial')
            sage: h, r
            (2*diff(u(x, y), x), 2/7*(diff(v(x, y), x, y) + 14)*diff(u(x, y), x))
            sage: C.normal_form (h * poly - r)
            0

        However, it is algebraically reduced w.r.t ``redset``::

            sage: h, r = R.differential_prem (poly, redset, mode = 'algebraic')
            sage: h, r
            (1, diff(u(x, y), x, x) + 1/7*diff(v(x, y), x, y))
            sage: C.normal_form (h * poly - r)
            0
        """
        cdef result
        cdef bytes streqns, strredset, mesgerr
        cdef bmi_c.ALGEB_string L
        streqns = bytes (polynomial)
        streqns = self._translate_str(streqns)
        strredset = bytes (redset)
        strredset = self._translate_str(strredset)
        L = bmi_c.bmi_sage_differential_prem2 (
            streqns, strredset, mode, self.dring,
            BMI_IX_undefined, notation, 0, 0)
        if bmi_c.bmi_sage_is_error (L):
            mesgerr = bmi_c.bmi_sage_mesgerr (L)
            bmi_c.bmi_balsa_clear_ALGEB (L)
            raise RuntimeError, mesgerr
        result = self._eval_sage(L.value)
        bmi_c.bmi_balsa_clear_ALGEB (L)
        return result


    def normal_form (
            self,
            equations,
            char* notation = BMI_IX_undefined):
        r"""
        The normal form of a rational differential fraction.

        INPUT:

        - ``equations`` -- a rational fraction or a list of
                           rational fractions
        - ``notation``  -- (default: 'undefined') a string

        OUTPUT:

        The list of normal forms of the equations.

        NOTES:

        Mathematically, this function does not do anything.
        Were there more than one notation, it could perform
        changes of notations. At least, it is useful to
        recognize rational fractions of differential polynomials
        equal to 0.

        So far, there is a single notation available, because of
        limitations of function.

        EXAMPLES::

            sage: from sage.calculus.DifferentialAlgebra import DifferentialRing
            sage: t = var ('t')
            sage: u,v = function ('u,v')
            sage: R = DifferentialRing (derivations = [t], blocks = [u,v])

        A complicated rational fraction of differential polynomials::

            sage: eq = diff (u(t) + (v(t)^2 + t)/(3*u(t)+diff(v(t),t)), t)
            sage: L = R.integrate (eq, t)

        This should be zero::

            sage: zero = eq - diff (L[1], t)
            sage: zero
            -(6*u(t)*diff(u(t), t) + 2*v(t)*diff(v(t), t) + diff(u(t), t)*diff(v(t), t) + u(t)*diff(v(t), t, t) + 1)/(3*u(t) + diff(v(t), t)) + (2*v(t)*diff(v(t), t) + 1)/(3*u(t) + diff(v(t), t)) + (3*u(t)^2 + v(t)^2 + u(t)*diff(v(t), t) + t)*(3*diff(u(t), t) + diff(v(t), t, t))/(3*u(t) + diff(v(t), t))^2 - (v(t)^2 + t)*(3*diff(u(t), t) + diff(v(t), t, t))/(3*u(t) + diff(v(t), t))^2 + diff(u(t), t)
            sage: R.normal_form (zero)
            0
        """
        cdef result
        cdef arg
        cdef bytes streqns, strders, mesgerr
        cdef bmi_c.ALGEB_string A

        if isinstance (equations, list):
            streqns = bytes (equations)
        else:
            streqns = bytes ([equations])
        streqns = self._translate_str(streqns)
        strders = bytes ([1])
        A = bmi_c.bmi_sage_differentiate (
                    self.dring, streqns, int (False), strders,
                    BMI_IX_undefined, notation, 0, 0)
        if bmi_c.bmi_sage_is_error (A):
            mesgerr = bmi_c.bmi_sage_mesgerr (A)
            bmi_c.bmi_balsa_clear_ALGEB (A)
            raise RuntimeError, mesgerr
        result = self._eval_sage(A.value)
        if not isinstance (equations, list):
            result = result [0]
        bmi_c.bmi_balsa_clear_ALGEB (A)
        return result

    def sort (
            self,
            list polynomials,
            char* mode = BMI_IX_ascending,
            char* notation = BMI_IX_undefined):
        r"""
        Sort a list of polynomials according to their leading monomials

        INPUT:

        - ``polynomials`` -- a list of polynomials
        - ``mode``        -- (default: 'increasing') a string
        - ``notation``    -- (default: 'undefined') a string

        OUTPUT:

        The list polynomials, sorted by increasing leading monomial,
        if ``mode`` is ``'increasing'``, by decreasing leading monomial,
        if ``mode`` is ``'decreasing'``.

        The leading monomial of a nonzero numeric polynomial is 1.
        The leading monomial of a nonzero, non numeric polynomial,
        is the product of the leading rank of the polynomial, by the
        leading monomial of its initial. Leading monomials are
        compared lexicographically, w.r.t. the ranking.
        This ordering is extended to 0, which is considered lower
        than any other polynomial.

        NOTES:

        In the future, this function will be extended to rational
        fractions.

        So far, there is a single notation available, because of
        limitations of function.

        EXAMPLES::

            sage: from sage.calculus.DifferentialAlgebra import DifferentialRing
            sage: a,x,y = var ('a,x,y')
            sage: u,v = function ('u,v')
            sage: R = DifferentialRing (derivations = [x,y], blocks = [u,v,a], parameters = [a,v(y)])
            sage: L = [diff (v(y),y,y), diff (v(y),y,y)^3, 0, a*u(x,y)^2, x]

        The list ``L`` is sorted by increasing leading monomial::

            sage: R.sort (L)
            [0, x, diff(v(y), y, y), diff(v(y), y, y)^3, a*u(x, y)^2]
            sage: R.sort (L, 'ascending')
            [0, x, diff(v(y), y, y), diff(v(y), y, y)^3, a*u(x, y)^2]

        By decreasing leading monomial::

            sage: R.sort (L, 'descending')
            [a*u(x, y)^2, diff(v(y), y, y)^3, diff(v(y), y, y), x, 0]

        Another example::

            sage: L = [ u(x,y), u(x,y)*v(y), u(x,y)*v(y)^2, u(x,y)*v(y)*a ]
            sage: R.sort (L, 'ascending')
            [u(x, y), u(x, y)*v(y), a*u(x, y)*v(y), u(x, y)*v(y)^2]
            sage: R.sort (L, 'descending')
            [u(x, y)*v(y)^2, a*u(x, y)*v(y), u(x, y)*v(y), u(x, y)]
        """
        cdef list result
        cdef bytes streqns, mesgerr
        cdef bmi_c.ALGEB_string L
        streqns = bytes (polynomials)
        streqns = self._translate_str(streqns)
        L = bmi_c.bmi_sage_sort (
                streqns, mode, self.dring, BMI_IX_undefined, notation, 0, 0)
        if bmi_c.bmi_sage_is_error (L):
            mesgerr = bmi_c.bmi_sage_mesgerr (L)
            bmi_c.bmi_balsa_clear_ALGEB (L)
            raise RuntimeError, mesgerr
        result = self._eval_sage (L.value)
        bmi_c.bmi_balsa_clear_ALGEB (L)
        return result


    def RosenfeldGroebner (
            self,
            equations,
            BaseFieldExtension basefield = None,
            list attributes = [
                    'differential', 'autoreduced', 'primitive',
                    'squarefree', 'normalized', 'coherent'],
            char* singsol = BMI_IX_all,
            char* dimlb = BMI_IX_safecase,
            char* notation = BMI_IX_undefined,
            int timeout = 0,
            int memout = 0):
        r"""
        The Rosenfeld-Groebner algorithm.

        INPUT:

        - ``equations``  -- a list of differential polynomials, differential
                            polynomial equations or inequations
        - ``basefield``  -- (default: None) a BaseFieldExtension
        - ``attributes`` -- (default: ['differential', 'autoreduced',
                            'primitive', 'squarefree', 'normalized',
                            'coherent']) a list of strings
        - ``singsol``    -- (default: 'all') a string
        - ``dimlb``      -- (default: 'safecase') a string
        - ``notation``   -- (default: 'undefined') a string
        - ``timeout``    -- (default: 0) a nonnegative integer
        - ``memout``     -- (default: 0) a nonnegative integer

        OUTPUT:

        A representation of the radical of the differential ideal generated
        by ``equations`` (assuming, for simplicity, that ``equations``
        only involves differential polynomials), as an intersection of
        radical differential ideals, presented by regular differential
        chains, with respect to the ranking of the 
        :class:`~sage.calculus.DifferentialAlgebra.DifferentialRing`.

        The list ``equations`` may involve differential polynomials,
        differential rational fractions, and, even relational expressions
        between differential rational fractions, built using the ``==`` and
        ``!=`` operators. In the presence of inequations (denominators of
        rational fractions, or relational expression built using ``!=``),
        the defined ideal is the radical of the differential ideal
        generated by the equations of ``equations``, saturated by
        the multiplicative family generated by the inequations of
        ``equations``.

        If ``basefield`` is not specified, then any case implying a
        relations between independent variables is discarded. If
        ``basefield`` is specified, then, moreover, any case implying
        the vanishing of some nonzero element of ``basefield`` is
        discarded also.

        The returned list may be empty, meaning the ideal is the unit one
        and that ``equations`` has no solution.

        If ``equations`` is empty then the output involves a single
        regular differential chain, which describes the zero ideal.

        The optional parameter ``attributes`` permit to specify the
        attributes of the regular differential chains to be computed.

        The optional parameter ``singsol`` = 'all', 'essential', 'none'.
        It controls the splittings performed by RosenfeldGroebner.

        - 'all' is the default value.

        - 'essential' ensures that the returned decomposition is
          irredundant. It only applies in the case of a single
          equation, without any parameter.

        - 'none' makes RosenfeldGroebner return at most one regular
          differential chain. This chain is the first one that would have
          been computed without this option. In many cases, this
          first chain could be considered as the general component,
          though this notion is not always mathematically well-defined.

        The optional parameter ``dimlb`` = 'nocase', 'safecase', 'odecase',
        'pdecase' controls the splittings performed by RosenfeldGroebner
        by discarding any chain whose dimension (differential dimension
        in the differential case) is lower than the number of input
        equations.

        - 'nocase' disables this option.

        - 'safecase' is the default value. The option is only activated
          in the non-differential case and in the case of a single input
          equation. These two cases are theoretically proven.

        - 'odecase' applies also the option to general ODE systems. The
          option implements a conjecture in this case (see Ritt's book,
          Questions for investigation, 10).

        - 'pdecase' applies also the option to general PDE systems. The
          option implements a conjecture in this case. It is still
          experimental.

        If not defined using ``notation``, the default notation of the
        returned regular differential chains is the one of the
        DifferentialRing.

        ``timeout`` and ``memout`` permit to limit in seconds and
        megabytes, the resources allowed for the computation.

        NOTES:

        The first version of the RosenfeldGroebner algorithm
        appeared in 1994. It benefited from many improvements,
        coming from many different people. The current version
        is close to the one described in [Boulier06]. It does
        not rely on any Gröbner basis computation.

        The RosenfeldGroebner method contains also an implementation
        of the Low Power Theorem, which is a major result of Ritt
        and Kolchin. See also [Hubert99]. The implemented version
        applies over differential base fields defined by generators
        and relations.

        REFERENCES:

        [Boulier06] François Boulier. Réécriture algébrique
        dans les systèmes d'équations différentielles polynomiales
        en vue d'application dans les Sciences du Vivant. Mémoire
        d'Habilitation à Diriger les Recherches. Université Lille I.
        2006. http://tel.archives-ouvertes.fr/tel-00137153

        [Hubert99] Évelyne Hubert. Essential Components of an
        Algebraic Differential Equation. Journal of Symbolic
        Computation 28(4-5), pp. 657-680. 1999.

        EXAMPLES:

        The first example features one nonlinear ODE. According to
        the Low Power Theorem, its solution `y(x) = 0` is a singular
        solution::

            sage: from sage.calculus.DifferentialAlgebra import DifferentialRing
            sage: x = var('x')
            sage: y = function ('y')
            sage: R = DifferentialRing (derivations = [x], blocks = [y])
            sage: L = R.RosenfeldGroebner ([diff(y(x),x)^2-4*y(x)], singsol = 'essential')
            sage: L
            [Regular differential chain (diff(y(x), x)^2 - 4*y(x)),
             Regular differential chain (y(x))]
            sage: L[0].equations ()
            [diff(y(x), x)^2 - 4*y(x)]
            sage: L[1].equations ()
            [y(x)]

        The second example features a system on nonlinear PDE. This
        academic example illustrates what a system of differential
        polynomials may look like, in general. Some integrability
        conditions are computed::

            sage: from sage.calculus.DifferentialAlgebra import DifferentialRing
            sage: x,y = var('x,y')
            sage: u,v = function('u,v')
            sage: R = DifferentialRing (derivations = [x,y], blocks = [[v,u]])
            sage: eq1 = diff(u(x,y),x)^2 - 4*u(x,y) == 0
            sage: eq2 = diff(diff(u(x,y),x),y) * diff(v(x,y),y) - u(x,y) + 1 == 0
            sage: eq3 = diff(diff(v(x,y),x),x) - diff(u(x,y),x) == 0
            sage: L = R.RosenfeldGroebner ([eq1, eq2, eq3])
            sage: L
            [Regular differential chain (diff(u(x, y), y)^2 - 2*u(x, y),
                 diff(u(x, y), x)^2 - 4*u(x, y),
                -u(x, y)*diff(u(x, y), x)*diff(u(x, y), y) + diff(u(x, y), x)*diff(u(x, y), y)
                     + 4*u(x, y)*diff(v(x, y), y),
                -diff(u(x, y), x) + diff(v(x, y), x, x))]
            sage: L[0].equations (solved=true)
            [diff(u(x, y), y)^2 == 2*u(x, y),
             diff(u(x, y), x)^2 == 4*u(x, y),
             diff(v(x, y), y) == 1/4*(u(x, y)*diff(u(x, y), x)*diff(u(x, y), y) - diff(u(x, y), x)*diff(u(x, y), y))/u(x, y),
             diff(v(x, y), x, x) == diff(u(x, y), x)]
        """
#
        cdef list result, eqns, ineq
        cdef BaseFieldExtension F
        cdef RegularDifferentialChain C
        cdef bmi_c.ALGEB_list L
        cdef bytes mesgerr, streqns, strineq, strgens, strrels, strattr
        if isinstance (equations, list):
            eqns, ineq = self.__process_equations (equations)
        else:
            eqns, ineq = self.__process_equations ([equations])
        streqns = bytes (eqns)
        strineq = bytes (ineq)
        strattr = string.replace (bytes (attributes), "'", "")
        streqns = self._translate_str(streqns)
        strineq = self._translate_str(strineq)
        if basefield == None:
            F = BaseFieldExtension ()
        else:
            F = basefield
        strgens = bytes (F.generators ())
        strrels = F.__relations ()
        strgens = self._translate_str(strgens)
        strrels = self._translate_str(strrels)
        L = bmi_c.bmi_sage_RosenfeldGroebner (
                    streqns, strineq, strgens, strrels,
                    strattr, self.dring, singsol, dimlb, 1,
                    BMI_IX_undefined, notation, timeout, memout)
        if bmi_c.bmi_sage_is_error (L):
            mesgerr = bmi_c.bmi_sage_mesgerr (L)
            bmi_c.bmi_balsa_clear_ALGEB (L)
            raise RuntimeError, mesgerr
        result = []
        for i in range (1, L.value.size) :
            C = RegularDifferentialChain.__new__ (RegularDifferentialChain)
            C.dring = self
            C.regchain = L.value.tab [i]
            bmi_c.bmi_balsa_increment_nbref (C.regchain)
            result.append (C)
        bmi_c.bmi_balsa_clear_ALGEB (L)
        return result

# __process_equations
    def __process_equations (self, list L):
        r"""
        Split a list into equations and inequations

        INPUT:

        - ``L`` -- a list

        OUTPUT:

        A pair of two lists: a first list of polynomials that should
        be considered as equations and a second one that should be
        considered as inequations.

        The elements of ``L`` can be rational differential fractions. Their
        denominators are then considered as inequations. The elements of
        ``L`` may also be relational expressions, built using ``==`` and ``!=``.

        NOTES:

        This function is used by :meth:`RosenfeldGroebner`, mostly.

        EXAMPLES::

            sage: from sage.calculus.DifferentialAlgebra import DifferentialRing
            sage: a,b = var ('a,b')
            sage: y = function ('y')
            sage: R = DifferentialRing (derivations = [a,b], blocks = [y])
            sage: L = [a/b, a+1/diff(y(a,b),a) == 0, b^2+b+1 != 0]
            sage: L
            [a/b, a + 1/diff(y(a, b), a) == 0, b^2 + b + 1 != 0]
            sage: R.__process_equations (L)
            ([a, a*diff(y(a, b), a) + 1], [b, diff(y(a, b), a), b^2 + b + 1])
        """
        cdef bmi_c.ALGEB_listof_string P
        cdef bytes L_as_string, mesgerr
        cdef list eqns, ineq
        L_as_string = bytes(L)
        L_as_string = self._translate_str(L_as_string)
        P = bmi_c.bmi_sage_process_equations (
                    L_as_string, self.dring,
                    BMI_IX_undefined, BMI_IX_undefined, 0, 0)
        if bmi_c.bmi_sage_is_error (P):
            mesgerr = bmi_c.bmi_sage_mesgerr (P)
            bmi_c.bmi_balsa_clear_ALGEB (P)
            raise RuntimeError, mesgerr

        eqns = self._eval_sage (P.value.tab [1].value)
        ineq = self._eval_sage (P.value.tab [2].value)

        bmi_c.bmi_balsa_clear_ALGEB (P)
        return (eqns, ineq)

    def __ranking (self):
        r"""
        The ranking of a DifferentialRing

        INPUT:

        Nothing

        OUTPUT:

        A string describing the ranking of the differential ring.

        NOTES:

        This function is used by the change_ranking method of
        the RegularDifferentialChain class.

        EXAMPLES::

            sage: from sage.calculus.DifferentialAlgebra import DifferentialRing
            sage: a,b = var ('a,b')
            sage: y = function ('y')
            sage: R = DifferentialRing (derivations = [a,b], blocks = [y])
            sage: R.__ranking ()
            'ranking (derivations = [a, b], blocks = [grlexA[y]])'
        """
        cdef bytes result
        cdef bytes mesgerr
        cdef bmi_c.ALGEB_string L
        L = bmi_c.bmi_sage_ranking (self.dring)
        if bmi_c.bmi_sage_is_error (L):
            mesgerr = bmi_c.bmi_sage_mesgerr (L)
            bmi_c.bmi_balsa_clear_ALGEB (L)
            raise RuntimeError, mesgerr
        result = bytes (L.value)
        bmi_c.bmi_balsa_clear_ALGEB (L)
        return result


cdef class RegularDifferentialChain:
    r"""
    The class RegularDifferentialChain implements differential regular chains.

        A regular differential chain is a set of polynomial differential
        equations in some simplified form. Regular differential chains
        belong to (polynomial) differential rings. They are usually
        produced, from raw systems of polynomial differential equations,
        by :class:`~meth.calculus.DifferentialAlgebra.DifferentialRing.RosenfeldGroebner`.

        Regular differential chains are slight generalizations of Ritt
        characteristic sets. A self-contained definition of regular
        differential chains is provided in [BL10].

    REFERENCES:

        [BL10] François Boulier and François Lemaire. A Normal
        Form Algorithm for Regular Differential Chains. Mathematics
        in Computer Science 4(2), pp. 185-201. 2010.
        http://dx.doi.org/10.1007/s11786-010-0060-3

    EXAMPLES:

    Create :class:`~sage.calculus.DifferentialAlgebra.DifferentialRing`,
    then a :class:`~sage.calculus.DifferentialAlgebra.RegularDifferentialChain`::

        sage: from sage.calculus.DifferentialAlgebra import DifferentialRing, RegularDifferentialChain
        sage: x,a,b = var ('x,a,b')
        sage: y = function ('y')
        sage: R = DifferentialRing (derivations = [x], blocks = [y,a,b], parameters = [a,b])
        sage: C = RegularDifferentialChain ([diff(y(x),x)^2 - a*y(x)], R)
        sage: C
        Regular differential chain (-a*y(x) + diff(y(x), x)^2)
    """

    cdef DifferentialRing dring
    """ The differential ring the chain belongs to """
    cdef bmi_c.ALGEB regchain
    """ Pointer to the BALSA internal data structure """
# __init__
    def __init__ (
            self,
            list equations,
            DifferentialRing DRing,
            list attributes = [
                'differential', 'autoreduced', 'primitive',
                'squarefree', 'normalized', 'coherent'],
            bool pretend = True,
            char* notation = BMI_IX_undefined,
            int timeout = 0,
            int memout = 0):
        r"""
        Constructor of RegularDifferentialChain

        INPUT:

        - ``equations``  -- a list of differential polynomial or
                            differential polynomial equations
        - ``DRing``      -- a DifferentialRing
        - ``attributes`` -- (default: ['differential', 'autoreduced',
                            'primitive', 'squarefree', 'normalized',
                            'coherent']) a list of strings
        - ``pretend``    -- (default: True) a boolean
        - ``notation``   -- (default: 'undefined') a string
        - ``timeout``    -- (default: 0) a nonnegative integer
        - ``memout``     -- (default: 0) a nonnegative integer

        OUTPUT:

        Create a regular differential chain of the differential
        polynomial ring ``DRing``, from ``equations``. The attributes
        of the chain depend on ``attributes``. By default, almost no test
        is performed, to check that the equations do constitute a chain,
        unless ``pretend`` is set to False.

        If not defined using ``notation``, the default notation is the
        one of ``DRing``.

        ``timeout`` and ``memout`` permit to limit in seconds and
        megabytes, the resources allowed for the computation.

        NOTES:

        So far, there is a single notation available, because of
        limitations of function.

        EXAMPLES:

        Create a regular differential chain with three PDE::

            sage: from sage.calculus.DifferentialAlgebra import DifferentialRing, RegularDifferentialChain
            sage: x,y = var ('x,y')
            sage: u,v = function ('u,v')
            sage: R = DifferentialRing (derivations = [x,y], blocks = [u,v])
            sage: eq1 = diff(u(x,y),x) == v(x,y)
            sage: eq2 = diff(u(x,y),y) == 0
            sage: eq3 = diff(v(x,y),y) == 0
            sage: C = RegularDifferentialChain ([eq1,eq2,eq3], R)
            sage: C.equations (solved=true)
            [diff(v(x, y), y) == 0, diff(u(x, y), y) == 0, diff(u(x, y), x) == v(x, y)]
        """
        cdef list eqns, ineq
        cdef bmi_c.ALGEB A
        cdef bytes streqns, strattr, mesgerr
        eqns, ineq = DRing.__process_equations (equations)
        streqns = bytes (eqns)
        streqns = DRing._translate_str (streqns)
        strattr = string.replace (bytes (attributes), "'", "")
        A = bmi_c.bmi_sage_pretend_regular_differential_chain (
                streqns, DRing.dring, strattr, int (pretend),
                BMI_IX_undefined, notation, timeout, memout)
        if bmi_c.bmi_sage_is_error (A):
            mesgerr = bmi_c.bmi_sage_mesgerr (A)
            bmi_c.bmi_balsa_clear_ALGEB (A)
            raise RuntimeError, mesgerr
        self.dring = DRing
        self.regchain = A
        bmi_c.bmi_balsa_increment_nbref (self.regchain)

    def __dealloc__ (self):
        """ The destructor """
        bmi_c.bmi_balsa_decrement_nbref (self.regchain)
        bmi_c.bmi_balsa_clear_ALGEB (self.regchain)

    def __repr__ (self):
        """ The external representation """
        eqs = self.equations()
        attrs = self.attributes()

        if 'differential' in attrs:
            attrstr = 'differential '
        else:
            attrstr = ''

        if len(eqs) > 0:
            return 'Regular ' + attrstr + 'chain (' + ', '.join([repr(eq) for eq in eqs]) + ')'
        else:
            return 'Regular ' + attrstr + 'chain (1)'

    def _latex_ (self):
        """ The LaTeX representation """
        eqs = self.equations()
        if len(eqs) > 0:
            return '\\left(' + ', '.join([latex(eq) for eq in eqs]) + '\\right)'
        else:
            return '(1)'

# differential_ring
    def differential_ring (self):
        r"""
        The :class:`~sage.calculus.DifferentialAlgebra.DifferentialRing` the chain belongs to

        INPUT:

        Nothing

        OUTPUT:

        A :class:`~sage.calculus.DifferentialAlgebra.DifferentialRing`

        EXAMPLES::

            sage: from sage.calculus.DifferentialAlgebra import DifferentialRing, RegularDifferentialChain
            sage: x,y = var ('x,y')
            sage: u,v = function ('u,v')
            sage: R = DifferentialRing (derivations = [x,y], blocks = [u,v])
            sage: eq1 = diff(u(x,y),x) == v(x,y)
            sage: eq2 = diff(u(x,y),y) == 0
            sage: eq3 = diff(v(x,y),y) == 0
            sage: C = RegularDifferentialChain ([eq1,eq2,eq3], R)

        Here, one recovers the ring from the chain::

            sage: S = C.differential_ring ()
            sage: S
            Differential Ring over Rational Field
                with dependents u(x, y), v(x, y) and independents x, y
            sage: S.indets ()
            [u(x, y), v(x, y)]
        """
        return self.dring

    def equations (
            self,
            bool solved = False,
            selection = None,
            char* notation = BMI_IX_undefined):
        r"""
        The list of equations of a regular differential chain.

        INPUT:

        - ``solved``    -- (default: False) a boolean
        - ``selection`` -- (default: None) a relational expression
        - ``notation``  -- (default: 'undefined') a string

        OUTPUT:

        A list of the differential polynomials which constitute the
        chain, sorted by increasing rank.

        If ``solved`` is ``True``, the differential polynomials are
        displayed as equations, solved with respect to their
        leading rank.

        The optional parameter ``selection`` is a relational
        expression combining one of the keywords ``'order'``, ``'leader'``,
        ``'rank'`` on the one side, with a value on the other side,
        using one of the relational operators. Only the equations
        which satisfy the condition are selected.

        - If the keyword is ``'order'``, then the value must be a
          nonnegative integer.

        - If the keyword is ``'rank'``, then the value must be a
          rank, including the special ranks 0 and 1.

        - If the keyword is ``'leader'``, and the relational
          operator is ``>=``, ``>``, ``<=`` or ``<``, then the value
          must be a derivative.

        - If the keyword is ``'leader'``, and the relational
          operator is ``==`` or ``!=``, then the value may also have
          the form: 'modifier' (derivative), where 'modifier' is
          one of the keywords ``'derivative'`` or ``'proper'``.

        If not defined using ``notation``, the notation is the
        chain default notation.

        NOTES:

        So far, there is a single notation available, because of
        limitations of function.

        EXAMPLES::

            sage: from sage.calculus.DifferentialAlgebra import DifferentialRing, RegularDifferentialChain
            sage: x,a,b = var('x,a,b')
            sage: y = function ('y')
            sage: R = DifferentialRing (derivations = [x], blocks = [y,[a,b]], parameters = [a,b])
            sage: C = RegularDifferentialChain ([diff(y(x),x)^2-a*y(x)], R)
            sage: C.equations ()
            [-a*y(x) + diff(y(x), x)^2]
            sage: C.equations (solved = true)
            [diff(y(x), x)^2 == a*y(x)]

        The next example illustrates the use of ``selection``::

            sage: from sage.calculus.DifferentialAlgebra import DifferentialRing, RegularDifferentialChain
            sage: x,y = var ('x,y')
            sage: u,v = function ('u,v')
            sage: R = DifferentialRing (derivations = [x,y], blocks = [u,v])
            sage: eq1 = diff(u(x,y),x) == v(x,y)
            sage: eq2 = diff(u(x,y),y) == 0
            sage: eq3 = diff(v(x,y),y) == 0
            sage: C = RegularDifferentialChain ([eq1,eq2,eq3], R)
            sage: C.equations (solved=True)
            [diff(v(x, y), y) == 0, diff(u(x, y), y) == 0, diff(u(x, y), x) == v(x, y)]
            sage: leader,rank,order = var ('leader,rank,order')
            sage: derivative,proper = function ('derivative,proper')

        The equations of order greater than `0`. Observe that the ranking
        is not necessarily orderly::

            sage: C.equations (selection = order > 0)
            [diff(v(x, y), y), diff(u(x, y), y), -v(x, y) + diff(u(x, y), x)]

        The equation whose leader is ``diff (u(x,y),x)``::

            sage: C.equations (selection = leader == diff (u(x,y),x))
            [-v(x, y) + diff(u(x, y), x)]

        The equations, whose leaders are derivatives of ``diff (u(x,y),x)``::

            sage: C.equations (selection = leader == derivative (diff (u(x,y),x)))
            [-v(x, y) + diff(u(x, y), x)]
        """
        cdef list result
        cdef bytes strsel, mesgerr
        cdef bmi_c.ALGEB_string A
        if selection == None:
            A = bmi_c.bmi_sage_equations (
                    self.regchain, int (False), int (solved),
                    BMI_IX_undefined, notation, 0, 0)
        else:
            strsel = bytes (selection)
            strsel = self.dring._translate_str(strsel)
            A = bmi_c.bmi_sage_equations_with_criterion_RDC (
                    self.regchain, int (False), int (solved), strsel,
                    BMI_IX_undefined, notation, 0, 0)
        if bmi_c.bmi_sage_is_error (A):
            mesgerr = bmi_c.bmi_sage_mesgerr (A)
            bmi_c.bmi_balsa_clear_ALGEB (A)
            raise RuntimeError, mesgerr

        result = self.dring._eval_sage (A.value)

        bmi_c.bmi_balsa_clear_ALGEB (A)
        return result

    def attributes (self):
        r"""
        The list of attributes of the :class:`Parent`.

        INPUT:

        Nothing

        OUTPUT:

        A list of strings, providing the attributes of
        a regular differential chain

        EXAMPLES::

            sage: from sage.calculus.DifferentialAlgebra import DifferentialRing, RegularDifferentialChain
            sage: x,y = var ('x,y')
            sage: u,v = function ('u,v')
            sage: R = DifferentialRing (derivations = [x,y], blocks = [u,v])
            sage: eq1 = diff(u(x,y),x) == v(x,y)
            sage: eq2 = diff(u(x,y),y) == 0
            sage: eq3 = diff(v(x,y),y) == 0
            sage: C = RegularDifferentialChain ([eq1,eq2,eq3], R)
            sage: C.attributes ()
            ['differential', 'prime', 'autoreduced', 'primitive', 'squarefree', 'coherent', 'normalized']
        """
        cdef list result
        cdef bytes mesgerr
        cdef dict dico
        cdef bmi_c.ALGEB_string A
        A = bmi_c.bmi_sage_attributes (self.regchain)
        if bmi_c.bmi_sage_is_error (A):
            mesgerr = bmi_c.bmi_sage_mesgerr (A)
            bmi_c.bmi_balsa_clear_ALGEB (A)
            raise RuntimeError, mesgerr
        dico = dict (
                [
                        ('differential', 'differential'),
                        ('prime','prime'),
                        ('autoreduced','autoreduced'),
                        ('squarefree','squarefree'),
                        ('coherent','coherent'),
                        ('primitive','primitive'),
                        ('normalized','normalized')
                ])
        result = eval (A.value, dico)
        bmi_c.bmi_balsa_clear_ALGEB (A)
        return result

# leading_derivative
    def leading_derivative (
            self,
            bool fullset = False,
            char* notation = BMI_IX_undefined):
        r"""
        The leading derivatives of a regular differential chain

        INPUT:

        - ``notation``  -- (default: 'undefined') a string

        OUTPUT:

        The list of the leading derivatives of the elements
        of the chain. See the documentation of :class:`~sage.calculus.DifferentialAlgebra.DifferentialRing`.

        NOTES:

        So far, there is a single notation available, because of
        limitations of function.

        EXAMPLES::

            sage: from sage.calculus.DifferentialAlgebra import DifferentialRing, RegularDifferentialChain
            sage: x,y = var ('x,y')
            sage: u,v = function ('u,v')
            sage: R = DifferentialRing (derivations = [x,y], blocks = [[v,u]])
            sage: eq1 = diff(u(x,y),x)^2 - 4*u(x,y)
            sage: eq2 = diff(u(x,y),x,y)*diff(v(x,y),y) - u(x,y) + 1
            sage: eq3 = diff(v(x,y),x,x) - diff(u(x,y),x)
            sage: ideal = R.RosenfeldGroebner ([eq1,eq2,eq3])
            sage: ideal
            [Regular differential chain (diff(u(x, y), y)^2 - 2*u(x, y),
                 diff(u(x, y), x)^2 - 4*u(x, y),
                -u(x, y)*diff(u(x, y), x)*diff(u(x, y), y) + diff(u(x, y), x)*diff(u(x, y), y)
                     + 4*u(x, y)*diff(v(x, y), y),
                 -diff(u(x, y), x) + diff(v(x, y), x, x))]
            sage: C = ideal [0]
            sage: C.leading_derivative()
            [diff(u(x, y), y), diff(u(x, y), x), diff(v(x, y), y), diff(v(x, y), x, x)]
        """
        cdef result
        cdef bytes streqns, mesgerr
        cdef bmi_c.ALGEB_string A
        streqns = bytes ("")
        A = bmi_c.bmi_sage_leading_derivative (
                    streqns, self.regchain, int (fullset),
                    BMI_IX_undefined, notation, 0, 0)
        if bmi_c.bmi_sage_is_error (A):
            mesgerr = bmi_c.bmi_sage_mesgerr (A)
            bmi_c.bmi_balsa_clear_ALGEB (A)
            raise RuntimeError, mesgerr
        result = self.dring._eval_sage(A.value)
        bmi_c.bmi_balsa_clear_ALGEB (A)
        return result

# leading_rank
    def leading_rank (
            self,
            bool fullset = False,
            bool listform = False,
            char* notation = BMI_IX_undefined):
        r"""
        The leading ranks of a regular differential chain

        INPUT:

        - ``listform``  -- (default: False) a bool
        - ``notation``  -- (default: 'undefined') a string

        OUTPUT:

        The list of the leading ranks of the elements
        of the chain. See the documentation of
        :class:`~sage.calculus.DifferentialAlgebra.DifferentialRing`.

        If the optional parameter ``listform`` is ``True``, then
        the ranks are returned as lists rather than monomials.

        NOTES:

        So far, there is a single notation available, because of
        limitations of function.

        EXAMPLES::

            sage: from sage.calculus.DifferentialAlgebra import DifferentialRing, RegularDifferentialChain
            sage: x,y = var ('x,y')
            sage: u,v = function ('u,v')
            sage: R = DifferentialRing (derivations = [x,y], blocks = [[v,u]])
            sage: eq1 = diff(u(x,y),x)^2 - 4*u(x,y)
            sage: eq2 = diff(u(x,y),x,y)*diff(v(x,y),y) - u(x,y) + 1
            sage: eq3 = diff(v(x,y),x,x) - diff(u(x,y),x)
            sage: ideal = R.RosenfeldGroebner ([eq1,eq2,eq3])
            sage: ideal
            [Regular differential chain (diff(u(x, y), y)^2 - 2*u(x, y),
                 diff(u(x, y), x)^2 - 4*u(x, y),
                -u(x, y)*diff(u(x, y), x)*diff(u(x, y), y) + diff(u(x, y), x)*diff(u(x, y), y)
                     + 4*u(x, y)*diff(v(x, y), y),
                -diff(u(x, y), x) + diff(v(x, y), x, x))]
            sage: C = ideal [0]
            sage: C.leading_rank()
            [diff(u(x, y), y)^2, diff(u(x, y), x)^2, diff(v(x, y), y), diff(v(x, y), x, x)]
        """
        cdef result
        cdef bytes streqns, mesgerr
        cdef bmi_c.ALGEB_string A
        streqns = bytes ("")
        A = bmi_c.bmi_sage_leading_rank (
                    streqns, self.regchain, int (fullset), int (listform),
                    BMI_IX_undefined, notation, 0, 0)
        if bmi_c.bmi_sage_is_error (A):
            mesgerr = bmi_c.bmi_sage_mesgerr (A)
            bmi_c.bmi_balsa_clear_ALGEB (A)
            raise RuntimeError, mesgerr
        result = self.dring._eval_sage(A.value)
        bmi_c.bmi_balsa_clear_ALGEB (A)
        return result

# initial
    def initial (
            self,
            bool fullset = False,
            char* notation = BMI_IX_undefined):
        r"""
        The initials of a regular differential chain

        INPUT:

        - ``notation``  -- (default: 'undefined') a string

        OUTPUT:

        The list of initials of the regular differential chain.
        See the documentation of :class:`~sage.calculus.DifferentialAlgebra.DifferentialRing`.

        NOTES:

        So far, there is a single notation available, because of
        limitations of function.

        EXAMPLES::

            sage: from sage.calculus.DifferentialAlgebra import DifferentialRing, RegularDifferentialChain
            sage: x,y = var ('x,y')
            sage: u,v = function ('u,v')
            sage: R = DifferentialRing (derivations = [x,y], blocks = [[v,u]])
            sage: eq1 = diff(u(x,y),x)^2 - 4*u(x,y)
            sage: eq2 = diff(u(x,y),x,y)*diff(v(x,y),y) - u(x,y) + 1
            sage: eq3 = diff(v(x,y),x,x) - diff(u(x,y),x)
            sage: ideal = R.RosenfeldGroebner ([eq1,eq2,eq3])
            sage: ideal
            [Regular differential chain (diff(u(x, y), y)^2 - 2*u(x, y),
                 diff(u(x, y), x)^2 - 4*u(x, y),
                -u(x, y)*diff(u(x, y), x)*diff(u(x, y), y) + diff(u(x, y), x)*diff(u(x, y), y)
                     + 4*u(x, y)*diff(v(x, y), y),
                -diff(u(x, y), x) + diff(v(x, y), x, x))]
            sage: C = ideal [0]
            sage: C.initial ()
            [1, 1, 4*u(x, y), 1]
        """
        cdef result
        cdef bytes streqns, strvar, mesgerr
        cdef bmi_c.ALGEB_string A
        streqns = bytes ("")
        strvar = bytes (0)
        A = bmi_c.bmi_sage_leading_coefficient (
                    streqns, self.regchain, int (fullset), strvar,
                    BMI_IX_undefined, notation, 0, 0)
        if bmi_c.bmi_sage_is_error (A):
            mesgerr = bmi_c.bmi_sage_mesgerr (A)
            bmi_c.bmi_balsa_clear_ALGEB (A)
            raise RuntimeError, mesgerr
        result = self.dring._eval_sage(A.value)
        bmi_c.bmi_balsa_clear_ALGEB (A)
        return result

# leading_coefficient
    def leading_coefficient (
            self,
            variable = None,
            bool fullset = False,
            char* notation = BMI_IX_undefined):
        r"""
        The leading coefficients of the regular differential chain.

        INPUT:

        - ``variable``  -- (default: None) a variable
        - ``notation``  -- (default: 'undefined') a string

        OUTPUT:

        The list of leading coefficients w.r.t. ``variable``
        of the regular differential chain. If ``variable``
        is omitted, the leading coefficients are taken w.r.t.
        the leading derivatives. See the documentation of
        :class:`~sage.calculus.DifferentialAlgebra.DifferentialRing`.

        NOTES:

        So far, there is a single notation available, because of
        limitations of function.

        EXAMPLES::

            sage: from sage.calculus.DifferentialAlgebra import DifferentialRing, RegularDifferentialChain
            sage: x,y = var ('x,y')
            sage: u,v = function ('u,v')
            sage: R = DifferentialRing (derivations = [x,y], blocks = [[v,u]])
            sage: eq1 = diff(u(x,y),x)^2 - 4*u(x,y)
            sage: eq2 = diff(u(x,y),x,y)*diff(v(x,y),y) - u(x,y) + 1
            sage: eq3 = diff(v(x,y),x,x) - diff(u(x,y),x)
            sage: ideal = R.RosenfeldGroebner ([eq1,eq2,eq3])
            sage: ideal
            [Regular differential chain (diff(u(x, y), y)^2 - 2*u(x, y),
                 diff(u(x, y), x)^2 - 4*u(x, y),
                -u(x, y)*diff(u(x, y), x)*diff(u(x, y), y) + diff(u(x, y), x)*diff(u(x, y), y)
                     + 4*u(x, y)*diff(v(x, y), y),
                -diff(u(x, y), x) + diff(v(x, y), x, x))]
            sage: C = ideal [0]

        If ``variable`` is omitted, the initials are returned::

            sage: C.leading_coefficient ()
            [1, 1, 4*u(x, y), 1]

        The leading coefficients w.r.t. ``u(x,y)``::

            sage: C.leading_coefficient (u(x,y))
            [-2,
             -4,
             -diff(u(x, y), x)*diff(u(x, y), y) + 4*diff(v(x, y), y),
             -diff(u(x, y), x) + diff(v(x, y), x, x)]
        """
        cdef result
        cdef bytes streqns, strvar, mesgerr
        cdef bmi_c.ALGEB_string A
        streqns = bytes ("")
        if variable == None:
            strvar = bytes (0)
        else:
            strvar = bytes (variable)
        A = bmi_c.bmi_sage_leading_coefficient (
                    streqns, self.regchain, int (fullset), strvar,
                    BMI_IX_undefined, notation, 0, 0)
        if bmi_c.bmi_sage_is_error (A):
            mesgerr = bmi_c.bmi_sage_mesgerr (A)
            bmi_c.bmi_balsa_clear_ALGEB (A)
            raise RuntimeError, mesgerr
        result = self.dring._eval_sage(A.value)
        bmi_c.bmi_balsa_clear_ALGEB (A)
        return result

# tail
    def tail (
            self,
            variable = None,
            bool fullset = False,
            char* notation = BMI_IX_undefined):
        r"""
        The tails of the regular differential chain.

        INPUT:

        - ``variable``  -- (default: None) a variable
        - ``notation``  -- (default: 'undefined') a string

        OUTPUT:

        The list of the tails, w.r.t. ``variable``, of the
        regular differential chain. If ``variable`` is omitted,
        the tails are taken w.r.t. the leading derivatives.
        See the documentation of :class:`~sage.calculus.DifferentialAlgebra.DifferentialRing`.

        NOTES:

        So far, there is a single notation available, because of
        limitations of function.

        EXAMPLES::

            sage: from sage.calculus.DifferentialAlgebra import DifferentialRing, RegularDifferentialChain
            sage: x,y = var ('x,y')
            sage: u,v = function ('u,v')
            sage: R = DifferentialRing (derivations = [x,y], blocks = [[v,u]])
            sage: eq1 = diff(u(x,y),x)^2 - 4*u(x,y)
            sage: eq2 = diff(u(x,y),x,y)*diff(v(x,y),y) - u(x,y) + 1
            sage: eq3 = diff(v(x,y),x,x) - diff(u(x,y),x)
            sage: ideal = R.RosenfeldGroebner ([eq1,eq2,eq3])
            sage: ideal
            [Regular differential chain (diff(u(x, y), y)^2 - 2*u(x, y),
                 diff(u(x, y), x)^2 - 4*u(x, y),
                -u(x, y)*diff(u(x, y), x)*diff(u(x, y), y) + diff(u(x, y), x)*diff(u(x, y), y)
                     + 4*u(x, y)*diff(v(x, y), y),
                 -diff(u(x, y), x) + diff(v(x, y), x, x))]
            sage: C = ideal [0]

        If ``variable`` is omitted, the tails are taken w.r.t. the
        leading derivatives::

            sage: C.tail ()
            [-2*u(x, y),
             -4*u(x, y),
             -u(x, y)*diff(u(x, y), x)*diff(u(x, y), y) + diff(u(x, y), x)*diff(u(x, y), y),
             -diff(u(x, y), x)]

        The tails w.r.t. ``u(x,y)``::

            sage: C.tail (u(x,y))
            [diff(u(x, y), y)^2, diff(u(x, y), x)^2, diff(u(x, y), x)*diff(u(x, y), y), 0]
        """
        cdef result
        cdef bytes streqns, strvar, mesgerr
        cdef bmi_c.ALGEB_string A
        streqns = bytes ("")
        if variable == None:
            strvar = bytes (0)
        else:
            strvar = bytes (variable)
        A = bmi_c.bmi_sage_tail (
                    streqns, self.regchain, int (fullset), strvar,
                    BMI_IX_undefined, notation, 0, 0)
        if bmi_c.bmi_sage_is_error (A):
            mesgerr = bmi_c.bmi_sage_mesgerr (A)
            bmi_c.bmi_balsa_clear_ALGEB (A)
            raise RuntimeError, mesgerr
        result = self.dring._eval_sage(A.value)
        bmi_c.bmi_balsa_clear_ALGEB (A)
        return result


    def separant (
            self,
            variable = None,
            bool fullset = False,
            char* notation = BMI_IX_undefined):
        r"""
        The separants of the regular differential chain.

        INPUT:

        - ``variable``  -- (default: ``None``) a variable
        - ``notation``  -- (default: 'undefined') a string

        OUTPUT:

        The list of the separants, w.r.t. ``variable``, of
        the regular differential chain. If ``variable`` is
        omitted, the separants are taken w.r.t. the leading
        derivatives. See the documentation of :class:`~sage.calculus.DifferentialAlgebra.DifferentialRing`.

        NOTES:

        So far, there is a single notation available, because of
        limitations of function.

        EXAMPLES::

            sage: from sage.calculus.DifferentialAlgebra import DifferentialRing, RegularDifferentialChain
            sage: x,y = var ('x,y')
            sage: u,v = function ('u,v')
            sage: R = DifferentialRing (derivations = [x,y], blocks = [[v,u]])
            sage: eq1 = diff(u(x,y),x)^2 - 4*u(x,y)
            sage: eq2 = diff(u(x,y),x,y)*diff(v(x,y),y) - u(x,y) + 1
            sage: eq3 = diff(v(x,y),x,x) - diff(u(x,y),x)
            sage: ideal = R.RosenfeldGroebner ([eq1,eq2,eq3])
            sage: ideal
            [Regular differential chain (diff(u(x, y), y)^2 - 2*u(x, y),
                 diff(u(x, y), x)^2 - 4*u(x, y),
                -u(x, y)*diff(u(x, y), x)*diff(u(x, y), y)
                     + diff(u(x, y), x)*diff(u(x, y), y)
                     + 4*u(x, y)*diff(v(x, y), y),
                -diff(u(x, y), x) + diff(v(x, y), x, x))]
            sage: C = ideal [0]

        If ``variable`` is omitted, the separants are taken w.r.t. the
        leading derivatives::

            sage: C.separant ()
            [2*diff(u(x, y), y), 2*diff(u(x, y), x), 4*u(x, y), 1]

        The separants w.r.t. ``u(x,y)``::

            sage: C.separant (u(x,y))
            [-2, -4, -diff(u(x, y), x)*diff(u(x, y), y) + 4*diff(v(x, y), y), 0]
        """
        cdef result
        cdef bytes streqns, strvar, mesgerr
        cdef bmi_c.ALGEB_string A
        streqns = bytes ("")
        if variable == None:
            strvar = bytes (0)
        else:
            strvar = bytes (variable)
        A = bmi_c.bmi_sage_separant (
                    streqns, self.regchain, int (fullset), strvar,
                    BMI_IX_undefined, notation, 0, 0)
        if bmi_c.bmi_sage_is_error (A):
            mesgerr = bmi_c.bmi_sage_mesgerr (A)
            bmi_c.bmi_balsa_clear_ALGEB (A)
            raise RuntimeError, mesgerr
        result = self.dring._eval_sage(A.value)
        bmi_c.bmi_balsa_clear_ALGEB (A)
        return result


    def indets (
                self,
                char* selection = BMI_IX_derivs,
                bool fullset = False,
                char* notation = BMI_IX_undefined):
        r"""
        The indeterminates of the regular differential chain.

        INPUT:

        - ``selection``  -- (default: 'derivatives') a string
        - ``notation``   -- (default: 'undefined') a string

        OUTPUT:

        A list of variable and function, selected according to
        ``selection``. See the examples below.

        NOTES:

        So far, the keyword 'constants' is not allowed to build
        selections.

        So far, there is a single notation available, because of
        limitations of function.

        EXAMPLES::

            sage: from sage.calculus.DifferentialAlgebra import DifferentialRing, RegularDifferentialChain
            sage: x,y = var ('x,y')
            sage: u,v = function ('u,v')
            sage: R = DifferentialRing (derivations = [x,y], blocks = [[v,u]])
            sage: eq1 = diff(u(x,y),x)^2 - 4*u(x,y)
            sage: eq2 = diff(u(x,y),x,y)*diff(v(x,y),y) - u(x,y) + 1
            sage: eq3 = diff(v(x,y),x,x) - diff(u(x,y),x)
            sage: ideal = R.RosenfeldGroebner ([eq1,eq2,eq3])
            sage: ideal
            [Regular differential chain (diff(u(x, y), y)^2 - 2*u(x, y),
                 diff(u(x, y), x)^2 - 4*u(x, y),
                -u(x, y)*diff(u(x, y), x)*diff(u(x, y), y)
                     + diff(u(x, y), x)*diff(u(x, y), y)
                     + 4*u(x, y)*diff(v(x, y), y),
                -diff(u(x, y), x) + diff(v(x, y), x, x))]
            sage: C = ideal [0]

        With no selection, or with ``selection='derivatives'``, return the list of
        the derivatives occuring in the elements of ``C``::

            sage: C.indets ()
            [u(x, y),
             diff(u(x, y), x),
             diff(v(x, y), y),
             diff(v(x, y), x, x),
             diff(u(x, y), y)]
            sage: C.indets (selection = 'derivatives')
            [u(x, y),
             diff(u(x, y), x),
             diff(v(x, y), y),
             diff(v(x, y), x, x),
             diff(u(x, y), y)]

        With ``selection='dependent'``, return the list of the differential
        indeterminates, or dependent variables, occuring in the elements of ``C``::

            sage: C.indets (selection = 'dependent')
            [u(x, y), v(x, y)]

        With ``selection = 'parameters'``, return the list of the parameters
        occuring in the elements of $C$::

            sage: C.indets( selection = 'parameters' )
            []

        With ``selection = 'all'``, return the list of the independent and the dependent
        variables occuring in the elements of ``C``::

            sage: C.indets (selection = 'all')
            [u(x, y),
             diff(u(x, y), x),
             diff(v(x, y), y),
             diff(v(x, y), x, x),
             diff(u(x, y), y)]
        """
        cdef list result
        cdef bytes streqns, strfull, mesgerr
        cdef bmi_c.ALGEB_string L
        streqns = bytes ("")
        if fullset:
            strfull = bytes ("true")
        else:
            strfull = bytes ("false")
        L = bmi_c.bmi_sage_indets (
                        streqns, self.regchain, selection, strfull,
                        BMI_IX_undefined, notation, 0, 0)
        if bmi_c.bmi_sage_is_error (L):
            mesgerr = bmi_c.bmi_sage_mesgerr (L)
            bmi_c.bmi_balsa_clear_ALGEB (L)
            raise RuntimeError, mesgerr
        result = self.dring._eval_sage(L.value)
        bmi_c.bmi_balsa_clear_ALGEB (L)
        return result


    def differential_prem (
             self,
             polynomial,
             char* mode = BMI_IX_fully,
             char* notation = BMI_IX_undefined):
        r"""
        Ritt's reduction algorithm by a regular differential chain

        INPUT:

        - ``polynomial`` -- a polynomial
        - ``mode``       -- (default: 'full') a string
        - ``notation``   -- (default: 'undefined') a string

        OUTPUT:

        A pair `(h, r)` such that `h*p = r` modulo the differential
        ideal generated by the equations of the regular differential
        chain (denoting `p` for ``polynomial``).

        The optional argument ``mode = 'full'``, ``'partial'`` or ``'algebraic'``.
        If ``'full'``, the remainder ``r`` is fully reduced w.r.t. the
        chain equations. If ``'partial'``, the remainder ``r`` is partially
        reduced w.r.t. the chain equations. If ``'algebraic'``, the
        remainder is algebraically reduced w.r.t. the chain equations.

        In all cases, ``h`` is a product of powers of the initials and
        separants involved in the reduction process.

        NOTES:

        So far, there is a single notation available, because of
        limitations of function.

        EXAMPLES::

            sage: from sage.calculus.DifferentialAlgebra import DifferentialRing, RegularDifferentialChain
            sage: x,y = var ('x,y')
            sage: u,v = function ('u,v')
            sage: R = DifferentialRing (derivations = [x,y], blocks = [[v,u]])
            sage: eq1 = diff(u(x,y),x)^2 - 4*u(x,y)
            sage: eq2 = diff(u(x,y),x,y)*diff(v(x,y),y) - u(x,y) + 1
            sage: eq3 = diff(v(x,y),x,x) - diff(u(x,y),x)
            sage: ideal = R.RosenfeldGroebner ([eq1,eq2,eq3])
            sage: ideal
            [Regular differential chain (diff(u(x, y), y)^2 - 2*u(x, y),
                 diff(u(x, y), x)^2 - 4*u(x, y),
                -u(x, y)*diff(u(x, y), x)*diff(u(x, y), y)
                     + diff(u(x, y), x)*diff(u(x, y), y)
                     + 4*u(x, y)*diff(v(x, y), y),
                -diff(u(x, y), x) + diff(v(x, y), x, x))]
            sage: C = ideal [0]

        Here are the chain equations::

            sage: C.equations (solved = true)
            [diff(u(x, y), y)^2 == 2*u(x, y),
             diff(u(x, y), x)^2 == 4*u(x, y),
             diff(v(x, y), y) == 1/4*(u(x, y)*diff(u(x, y), x)*diff(u(x, y), y) - diff(u(x, y), x)*diff(u(x, y), y))/u(x, y),
             diff(v(x, y), x, x) == diff(u(x, y), x)]

        The first polynomial is already fully reduced w.r.t. the chain
        equations::

            sage: poly = (1/7)*u(x,y)
            sage: h, r = C.differential_prem (poly)
            sage: h, r
            (1, 1/7*u(x, y))
            sage: C.normal_form (h * poly - r)
            0

        The next polynomial is not reduced w.r.t. the chain equations::

            sage: poly = (1/7)*diff (v(x,y),x,y) + diff (u(x,y),x,x)
            sage: h, r = C.differential_prem (poly)
            sage: h, r
            (16*u(x, y)^2*diff(u(x, y), x)*diff(u(x, y), y),
             32/7*(u(x, y) + 7*diff(u(x, y), y))*u(x, y)^2*diff(u(x, y), x))
            sage: C.normal_form (h * poly - r)
            0

        Let us try a partial reduction, instead of a full one::

            sage: h, r = C.differential_prem (poly, mode = 'partial')
            sage: h, r
            (16*u(x, y)*diff(u(x, y), x)*diff(u(x, y), y),
             4/7*(diff(u(x, y), x)^2*diff(u(x, y), y)^2 + u(x, y)*diff(u(x, y), x)^2 + 2*u(x, y)*diff(u(x, y), y)^2 - 4*diff(u(x, y), x)*diff(u(x, y), y)*diff(v(x, y), y) - diff(u(x, y), x)^2 + 56*u(x, y)*diff(u(x, y), y) - 2*diff(u(x, y), y)^2)*diff(u(x, y), x))
            sage: C.normal_form (h * poly - r)
            0

        The polynomial actually is algebraically reduced w.r.t the chain::

            sage: h, r = C.differential_prem (poly, mode = 'algebraic')
            sage: h, r
            (1, diff(u(x, y), x, x) + 1/7*diff(v(x, y), x, y))
            sage: C.normal_form (h * poly - r)
            0
        """
        cdef result
        cdef bytes streqns, mesgerr
        cdef bmi_c.ALGEB_string L
        streqns = bytes (polynomial)
        streqns = self.dring._translate_str(streqns)
        L = bmi_c.bmi_sage_differential_prem (
            streqns, mode, self.regchain, BMI_IX_undefined, notation, 0, 0)
        if bmi_c.bmi_sage_is_error (L):
            mesgerr = bmi_c.bmi_sage_mesgerr (L)
            bmi_c.bmi_balsa_clear_ALGEB (L)
            raise RuntimeError, mesgerr
        result = self.dring._eval_sage(L.value)
        bmi_c.bmi_balsa_clear_ALGEB (L)
        return result

    def change_ranking (
             self,
             DifferentialRing R = None,
             bool prime = False,
             char* notation = BMI_IX_undefined,
             int timeout = 0,
             int memout = 0):
        r"""
        Change of ranking for regular differential chains defining prime ideals

        INPUT:

        - ``R``          -- (default: None) a differential ring
        - ``prime``      -- (default: False) a bool
        - ``notation``   -- (default: 'undefined') a string
        - ``timeout``    -- (default: 0) a nonnegative integer
        - ``memout``     -- (default: 0) a nonnegative integer

        OUTPUT:

        A regular differential chain defining the same ideal
        as ``self``, w.r.t. the ranking of ``R``.

        The differential ring ``R`` must define the same
        mathematical differential ring as the underlying
        ring of ``self`` (only the ranking should be different).

        The optional parameter ``prime`` permit to force the
        application of the algorithm over ``self``, even if
        ``self`` does not define a prime ideal.

        ``timeout`` and ``memout`` permit to limit in seconds and
        megabytes, the resources allowed for the computation.

        NOTES:

        The implemented algorithm is the PARDI algorithm. See
        references below.

        So far, there is a single notation available, because of
        limitations of function.

        REFERENCES:

        [BLM10] François Boulier, François Lemaire and Marc Moreno
        Maza. Computing differential characteristic sets by change
        of ordering. Journal of Symbolic Comp. 45(1), pp. 124-149. 2010

        EXAMPLES::

            sage: from sage.calculus.DifferentialAlgebra import DifferentialRing, RegularDifferentialChain
            sage: x,y = var ('x,y')
            sage: u,v = function ('u,v')
            sage: R = DifferentialRing (derivations = [x,y], blocks = [[v,u]])
            sage: eq1 = diff(u(x,y),x)^2 - 4*u(x,y)
            sage: eq2 = diff(u(x,y),x,y)*diff(v(x,y),y) - u(x,y) + 1
            sage: eq3 = diff(v(x,y),x,x) - diff(u(x,y),x)
            sage: ideal = R.RosenfeldGroebner ([eq1,eq2,eq3])
            sage: ideal
            [Regular differential chain (diff(u(x, y), y)^2 - 2*u(x, y),
                 diff(u(x, y), x)^2 - 4*u(x, y),
                -u(x, y)*diff(u(x, y), x)*diff(u(x, y), y)
                     + diff(u(x, y), x)*diff(u(x, y), y)
                     + 4*u(x, y)*diff(v(x, y), y),
                -diff(u(x, y), x) + diff(v(x, y), x, x))]

        The regular differential chain ``C`` defines a prime ideal.
        Its ranking is orederly::

            sage: C = ideal [0]
            sage: C.equations (solved = True)
            [diff(u(x, y), y)^2 == 2*u(x, y),
             diff(u(x, y), x)^2 == 4*u(x, y),
             diff(v(x, y), y) == 1/4*(u(x, y)*diff(u(x, y), x)*diff(u(x, y), y) - diff(u(x, y), x)*diff(u(x, y), y))/u(x, y),
             diff(v(x, y), x, x) == diff(u(x, y), x)]

        The ring ``Rbar`` is mathematically equivalent to ``R``, but w.r.t.
        an elimination ranking::

            sage: Rbar = DifferentialRing (derivations = [x,y], blocks = [u,v])

        Perform the change of ranking, claiming that the ideal is prime::

            sage: Cbar = C.change_ranking (Rbar, prime=True)
            sage: Cbar.equations (solved = True)
            [diff(v(x, y), y, y)^4 == 2*diff(v(x, y), y)^2 + 2*diff(v(x, y), y, y)^2 - 1,
             diff(v(x, y), x, y) == (diff(v(x, y), y, y)^3 - diff(v(x, y), y, y))/diff(v(x, y), y),
             diff(v(x, y), x, x) == 2*diff(v(x, y), y, y),
             u(x, y) == diff(v(x, y), y, y)^2]
        """
        cdef bmi_c.ALGEB L
        cdef bytes mesgerr, ranking
        cdef RegularDifferentialChain C
        ranking = string.replace (R.__ranking (), "'", "")
        L = bmi_c.bmi_sage_pardi (
                self.regchain, ranking, int (prime),
                BMI_IX_undefined, notation, timeout, memout)
        if bmi_c.bmi_sage_is_error (L):
            mesgerr = bmi_c.bmi_sage_mesgerr (L)
            bmi_c.bmi_balsa_clear_ALGEB (L)
            raise RuntimeError, mesgerr
        C = RegularDifferentialChain.__new__ (RegularDifferentialChain)
        C.dring = R
        C.regchain = L
        bmi_c.bmi_balsa_increment_nbref (C.regchain)
        bmi_c.bmi_balsa_clear_ALGEB (L)
        return C


    def normal_form (
            self,
            equations,
            bool casesplit = False,
            char* notation = BMI_IX_undefined):
        r"""
        The normal form of a rational differential fraction

        INPUT:

        - ``equations`` -- a rational fraction or a list of rational
                           fractions
        - ``casesplit`` -- (default: False) a bool
        - ``notation``  -- (default: 'undefined') a string

        OUTPUT:

        The normal form of the equation w.r.t. the regular
        differential chain. If ``equations`` is a list, a
        list of normal form is returned.

        Given a regular differential chain, not all rational
        fractions have normal forms: the denominator should not
        be zero, or even a zerod-ivisor, modulo the ideal defined
        by the regular differential chain.

        The optional argument ``casesplit`` permit to handle
        the case of rational fractions which do not have normal
        forms, by splitting cases.

        NOTES:

        The argument ``casesplit = True`` is not yet implemented.

        The current implementation of normal forms may fail, even
        if the normal form does exist (very rare but possible situation).
        A complete algorithm is available but not yet implemented.
        See references below.

        So far, there is a single notation available, because of
        limitations of function.

        REFERENCES:

        [BLS11] François Boulier, François Lemaire and Alexandre
        Sedoglavic. On the Regularity Property of Differential
        Polynomials Modulo Regular Differential Chains. In proc.
        of CASC 2011. LNCS 6885, pp. 61-72. 2011.
        http://hal.archives-ouvertes.fr/hal-00599440

        EXAMPLES::

            sage: from sage.calculus.DifferentialAlgebra import DifferentialRing, RegularDifferentialChain
            sage: x,y = var ('x,y')
            sage: u,v = function ('u,v')
            sage: R = DifferentialRing (derivations = [x,y], blocks = [[v,u]])
            sage: eq1 = diff(u(x,y),x)^2 - 4*u(x,y)
            sage: eq2 = diff(u(x,y),x,y)*diff(v(x,y),y) - u(x,y) + 1
            sage: eq3 = diff(v(x,y),x,x) - diff(u(x,y),x)
            sage: ideal = R.RosenfeldGroebner ([eq1,eq2,eq3])
            sage: ideal
            [Regular differential chain (diff(u(x, y), y)^2 - 2*u(x, y),
                 diff(u(x, y), x)^2 - 4*u(x, y),
                -u(x, y)*diff(u(x, y), x)*diff(u(x, y), y)
                     + diff(u(x, y), x)*diff(u(x, y), y)
                     + 4*u(x, y)*diff(v(x, y), y),
                -diff(u(x, y), x) + diff(v(x, y), x, x))]
            sage: C = ideal [0]

        A few basic examples::

            sage: p = diff(u(x,y),y,y) / diff(u(x,y),x)
            sage: C.normal_form (p)
            1/4*diff(u(x, y), x)/u(x, y)
            sage: C.normal_form (1/p)
            diff(u(x, y), x)

        An expected behaviour of normal forms::

            sage: C.normal_form (C.normal_form (p) * C.normal_form (1/p))
            1
        """
        cdef result, elt
        cdef bytes streqns, mesgerr
        cdef bmi_c.ALGEB_string A
        if casesplit:
            raise NotImplementedError, bytes ("casesplit = true")
        if isinstance (equations, list):
            streqns = bytes (equations)
        else:
            streqns = bytes ([equations])
        streqns = self.dring._translate_str(streqns)
        A = bmi_c.bmi_sage_normal_form (
                    self.regchain, streqns,
                    BMI_IX_undefined, notation, 0, 0)
        if bmi_c.bmi_sage_is_error (A):
            mesgerr = bmi_c.bmi_sage_mesgerr (A)
            bmi_c.bmi_balsa_clear_ALGEB (A)
            raise RuntimeError, mesgerr
        result = self.dring._eval_sage(A.value)
        result = [elt[0] for elt in result]
        if not isinstance (equations, list):
            result = result [0]
        bmi_c.bmi_balsa_clear_ALGEB (A)
        return result

    def preparation_equation (
            self,
            polynomial,
            BaseFieldExtension basefield = None,
            bool congruence = False,
            char* zstring = 'z%d',
            char* notation = BMI_IX_undefined):
        r"""
        The preparation equation of a polynomial.

        INPUT:

        - ``polynomial`` -- a polynomial
        - ``basefield``  -- (default: None) a BaseFieldExtension
        - ``congruence`` -- (default: False) a bool
        - ``zstring``    -- (default: 'z%d') a string
        - ``notation``   -- (default: 'undefined') a string

        OUTPUT:

        Return a preparation equation [Kolchin73, chapter IV, section 13]
        of ``polynomial`` w.r.t. the regular differential chain.

        Preparation equations are an important tool of the Low Power
        Theorem [Ritt50, chapter III, Kolchin73, chapter IV, section 15].
        The Low Power Theorem is applied by the
        :meth:`~sage.calculus.DifferentialAlgebra.DifferentialRing.RosenfeldGroebner`
        method, when there is a single equation to process and
        the optional parameter ``singsol = True``.

        Let `I` the differential ideal defined by the chain, denote
        ``A[1], ..., A[r]`` the differential polynomials which constitute
        the chain. Introduce `r` new dependent variables ``z[i]``. Each
        variable ``z[i]`` representing the differential polynomial ``A[i]``.

        The returned preparation equation is an expression having the
        form ``h*f = c[1]*t[1] + c[2]*t[2] + ... + c[n]*t[n]``, where ``f``
        stands for ``polynomial``. The differential polynomial ``h`` is
        a power product of initials and separants of the ``A[i]``. The
        coefficients ``c[i]`` are reduced w.r.t. the chain and regular
        with respect to `I`. The monomials ``t[i]`` are power products of
        the ``z[k]`` variables and their derivatives. They satisfy some
        further properties, described in [Kolchin73, chapter IV,
        section 13]. If each ``z[k]`` is replaced by the corresponding
        polynomial ``A[k]``, then the preparation equation becomes an
        equality.

        The optional parameter ``congruence``. In the right hand-side
        of the preparation equation, denote ``q`` the minimum total
        degree (the low power, actually) of the monomials ``t[i]``.
        If ``congruence`` is true, then all the terms ``c[i]*t[i]``
        such that the total degree of ``t[i]`` is strictly greater
        than ``q`` are removed. The remaining sum is call a preparation
        congruence of ``f``.

        The optional parameter ``basefield`` is useful in conjunction
        with ``congruence = true``. Reductions by the ``A[i]`` which
        belong to the base field are not taken into account for
        computing the preparation congruence of ``f``. The terms ``t[k]``
        which depend on such ``z[i]`` and their derivatives are not
        considered for computing the degree ``q`` and do not appear
        in the preparation congruence.

        The optional parameter ``zstring`` permits to denote the
        polynomial ``A[i]`` by symbols different from ``z[i]``. This
        may be useful since the method creates one function (in
        the SAGE sense) for each ``z[i]``. The ``zstring`` should
        correspond to a valid function identifier and involve
        the substring ``'%d'``.

        NOTES:

        So far, there is a single notation available, because of
        limitations of function.

        EXAMPLES::

            sage: from sage.calculus.DifferentialAlgebra import DifferentialRing, RegularDifferentialChain, BaseFieldExtension
            sage: t = var('t')
            sage: u = function('u')
            sage: R = DifferentialRing (derivations = [t], blocks = [u])

        A famous example of Ritt. Two cases are considered::

            sage: poly = diff(u(t),t)^2 - 4*u(t)
            sage: ideal = R.RosenfeldGroebner ([poly])
            sage: [ C.equations (solved = true) for C in ideal ]
            [[diff(u(t), t)^2 == 4*u(t)], [u(t) == 0]]

        Is the solution ``u(t) = 0`` an essential component::

            sage: C = ideal [1]

        Here is a preparation equation for poly, w.r.t. the chain ``C``.
        Replacing ``z1(t)`` by ``u(t)``, it is easy to see check the result::

            sage: C.preparation_equation (poly)
            diff(u(t), t)^2 - 4*u(t) == diff(z1(t), t)^2 - 4*z1(t)

        Here is the preparation congruence (``q = 1``)::

            sage: C.preparation_equation (poly, congruence = true)
            diff(u(t), t)^2 - 4*u(t) == -4*z1(t)

        By [Kolchin73, chapter IV, section 15, theorem 6], the
        solution is, indeed, essential. This result can be interpreted
        as follows: the solution ``u(t) = 0`` is not a particular case
        of the general solution `u(t) = (t + c)^2`, where `c` is an arbitrary
        constant::

            sage: ideal = R.RosenfeldGroebner ([poly], singsol = 'essential')
            sage: [ C.equations (solved = true) for C in ideal ]
            [[diff(u(t), t)^2 == 4*u(t)], [u(t) == 0]]

        Another famous example. Two similar cases are generated::

            sage: poly = diff(u(t),t)^2 - 4*u(t)^3
            sage: ideal = R.RosenfeldGroebner ([poly])
            sage: ideal
            [Regular differential chain (-4*u(t)^3 + diff(u(t), t)^2),
             Regular differential chain (u(t))]
            sage: [ C.equations (solved = true) for C in ideal ]
            [[diff(u(t), t)^2 == 4*u(t)^3], [u(t) == 0]]

        Is the solution `u(t) = 0` an essential component ?::

            sage: C = ideal [1]

        Here is a preparation equation for poly, w.r.t. the chain ``C``.
        Replacing ``z1(t)`` by ``u(t)``, it is easy to see check the result::

            sage: C.preparation_equation (poly)
            -4*u(t)^3 + diff(u(t), t)^2 == -4*z1(t)^3 + diff(z1(t), t)^2

        Here is the preparation congruence (`q = 2`)::

            sage: C.preparation_equation (poly, congruence = true)
            -4*u(t)^3 + diff(u(t), t)^2 == diff(z1(t), t)^2

        By Ritt's theorem, the solution is particular. This result
        can be interpreted as follows: the solution ``u(t) = 0`` is
        the limit, when the arbitrary constant `c` tends towards infinity,
        of the general solution `u(t) = 1/(t + c)^2`::

            sage: ideal = R.RosenfeldGroebner ([poly], singsol = 'essential')
            sage: [ C.equations (solved = true) for C in ideal ]
            [[diff(u(t), t)^2 == 4*u(t)^3]]

        The following example just illustrates the use of base fields.
        The F is a differential field involving two functions `s(t), c(t)`,
        which can be viewed as `\sin(t)` and `\cos(t)`::

            sage: t = var('t')
            sage: u,s,c = function ('u,s,c')
            sage: R = DifferentialRing (derivations = [t], blocks = [u,[s,c]])
            sage: C = RegularDifferentialChain ([diff (c(t),t) + s(t), s(t)^2 + c(t)^2 - 1], R)
            sage: F = BaseFieldExtension (relations = C)

        Observe that, over ``F``, ``poly`` and ``poly2`` are the same polynomial.
        They thus should have the same preparation congruence::

            sage: poly = diff (u(t),t)^2 - 4*u(t)^3
            sage: poly2 = diff (u(t),t)^2 - 4*u(t)^3 + s(t)^2 + c(t)^2 - 1

        Two cases are generated::

            sage: ideal = R.RosenfeldGroebner ([poly], basefield = F)
            sage: [ C.equations (solved = true) for C in ideal ]
            [[s(t)^2 == -c(t)^2 + 1, diff(c(t), t) == -s(t), diff(u(t), t)^2 == 4*u(t)^3], [s(t)^2 == -c(t)^2 + 1, diff(c(t), t) == -s(t), u(t) == 0]]

        Is the solution `u(t) = 0` essential?::

            sage: C = ideal[1]
            sage: C.equations (solved = true)
            [s(t)^2 == -c(t)^2 + 1, diff(c(t), t) == -s(t), u(t) == 0]

        The preparation congruence of poly (playing with ``zstring``)::

            sage: C.preparation_equation (poly, zstring = 'C%d')
            -4*u(t)^3 + diff(u(t), t)^2 == -4*C3(t)^3 + diff(C3(t), t)^2

        The preparation congruence of poly shows that the solution
        is not essential (see the previous example)::

            sage: C.preparation_equation (poly, congruence=true)
            -4*u(t)^3 + diff(u(t), t)^2 == diff(z3(t), t)^2

        Let us now consider ``poly2``, which is equal to poly, over ``F``::

            sage: C.preparation_equation (poly2)
            -4*u(t)^3 + c(t)^2 + s(t)^2 + diff(u(t), t)^2 - 1 == -4*z3(t)^3 + diff(z3(t), t)^2 + z1(t)

        The term ``t[i]`` of minimal degree correspond to a reduction
        by a field equation::

            sage: C.preparation_equation (poly2, congruence=true)
            -4*u(t)^3 + c(t)^2 + s(t)^2 + diff(u(t), t)^2 - 1 == z1(t)

        Specifying that we should consider ``poly2`` over ``F``, the preparation
        congruence changes and we conclude similarly as for poly::

            sage: C.preparation_equation (poly2, congruence=true, basefield=F)
            -4*u(t)^3 + c(t)^2 + s(t)^2 + diff(u(t), t)^2 - 1 == diff(z3(t), t)^2

        We check below that :meth:`~sage.calculus.DifferentialAlgebra.DifferentialRing.RosenfeldGroebner`
        recognizes the solution as not essential::

            sage: ideal = R.RosenfeldGroebner ([poly], singsol = 'essential', basefield = F)
            sage: [ C.equations (solved = true) for C in ideal ]
            [[s(t)^2 == -c(t)^2 + 1, diff(c(t), t) == -s(t), diff(u(t), t)^2 == 4*u(t)^3]]
        """
        cdef result, nbz, i
        cdef bytes mesgerr, strpoly, strgens, strrels, newval
        cdef bmi_c.ALGEB_string A

        strpoly = bytes (polynomial)
        strpoly = self.dring._translate_str(strpoly)
        if basefield == None:
            F = BaseFieldExtension ()
        else:
            F = basefield
        strgens = bytes (F.generators ())
        strrels = F.__relations ()
        strgens = self.dring._translate_str(strgens)
        strrels = self.dring._translate_str(strrels)
        A = bmi_c.bmi_sage_preparation_equation (
                strpoly, self.regchain, strgens, strrels, int (congruence),
                zstring, BMI_IX_undefined, notation, 0, 0)
        if bmi_c.bmi_sage_is_error (A):
            mesgerr = bmi_c.bmi_sage_mesgerr (A)
            bmi_c.bmi_balsa_clear_ALGEB (A)
            raise RuntimeError, mesgerr
        nbz = self.__number_of_equations ()
        local_names = dict()
        for i in range (nbz):
            local_names[zstring % (i+1)] = symbolic_function (zstring % (i+1))
        newval = bytes (A.value)
        newval = string.replace (newval, "=", "==")
        result = self.dring._eval_sage(newval, local_names)
        bmi_c.bmi_balsa_clear_ALGEB (A)
        return result

# __number_of_equations
    def __number_of_equations (self):
        r"""
        The number of equations of a regular differential chain

        INPUT:

        Nothing

        OUTPUT:

        The number of equations of a regular differential chain

        EXAMPLES::

            sage: from sage.calculus.DifferentialAlgebra import DifferentialRing, RegularDifferentialChain, BaseFieldExtension
            sage: t = var('t')
            sage: u = function('u')
            sage: R = DifferentialRing (derivations = [t], blocks = [u])
            sage: poly = diff(u(t),t)^2 - 4*u(t)
            sage: ideal = R.RosenfeldGroebner ([poly])
            sage: ideal
            [Regular differential chain (diff(u(t), t)^2 - 4*u(t)),
             Regular differential chain (u(t))]
            sage: C = ideal [0]
            sage: C.equations ()
            [diff(u(t), t)^2 - 4*u(t)]
            sage: C.__number_of_equations ()
            1
        """
        cdef bmi_c.ALGEB_string A
        cdef bytes mesgerr
        cdef result
        A = bmi_c.bmi_sage_number_of_equations (self.regchain)
        if bmi_c.bmi_sage_is_error (A):
            mesgerr = bmi_c.bmi_sage_mesgerr (A)
            bmi_c.bmi_balsa_clear_ALGEB (A)
            raise RuntimeError, mesgerr
        result = self.dring._eval_sage(A.value)
        bmi_c.bmi_balsa_clear_ALGEB (A)
        return result


# BaseFieldExtension
cdef class BaseFieldExtension:
    r"""
    The class BaseFieldExtension implements base field extensions of
    differential rings.

    These field extensions are defined by generators and relations.
    They permit to control the splittings performed by
    :meth:`~sage.calculus.DifferentialAlgebra.DifferentialRing.RosenfeldGroebner`.

    Observe that a base field extension `F` can be defined without
    specifying the differential ring `R = K\{U\}` whose base field `K`
    is going to be extended. Indeed, the actual extension of `K`
    is performed by the :meth:`~sage.calculus.DifferentialAlgebra.DifferentialRing.RosenfeldGroebner`
    algorithm.

    The generators of the base field extension of a differential ring
    must belong to the lowest blocks of the ranking (i.e. the block list).
    If the differential ring is not specified at the base field extension
    definition, the test is performed by the
    :meth:`~sage.calculus.DifferentialAlgebra.DifferentialRing.RosenfeldGroebner` algorithm.

    EXAMPLES:

    First create a differential ring with two parameters ``a``, ``b``.
    The base field of ``R`` is `K = \QQ(x)`::

        sage: from sage.calculus.DifferentialAlgebra import DifferentialRing, RegularDifferentialChain, BaseFieldExtension
        sage: x,a,b = var('x,a,b')
        sage: y = function ('y')
        sage: R = DifferentialRing (derivations = [x], blocks = [y,[a,b]], parameters = [a,b])

    Create now `F = \QQ(a,b)` a purely transcendental field extension
    of the field `\QQ`::

        sage: from sage.calculus.DifferentialAlgebra import DifferentialRing, RegularDifferentialChain, BaseFieldExtension
        sage: F = BaseFieldExtension (generators = [a,b])
        sage: F.generators ()
        [a, b]
        sage: F.relations ()
        Regular chain (1)
        sage: F.relations ().equations ()
        []

    Alternatively, since we know in advance that ``F`` is going to be used
    as a base field extension of the sole ring ``R``, it is safe to specify
    the ring as an optional parameter, in order to check that ``a`` and ``b``
    lie at the bottom of the ranking of ``R``::

        sage: F = BaseFieldExtension (generators = [a,b], ring = R)
        sage: F
        differential_field

    In the next example (``F`` is not yet used), the
    :meth:`~sage.calculus.DifferentialAlgebra.DifferentialRing.RosenfeldGroebner`
    algorithm considers three cases.::

        sage: L = R.RosenfeldGroebner ([diff(y(x),x)^2-a*y(x)])
        sage: [ C.equations () for C in L ]
        [[-a*y(x) + diff(y(x), x)^2], [a, diff(y(x), x)], [y(x)]]

    In this variant, the :meth:`~sage.calculus.DifferentialAlgebra.DifferentialRing.RosenfeldGroebner`
    algorithm is applied over the base field ``F(x) = K(a,b) = Q(a,b,x)``, obtained by
    extending the base field ``K`` of ``R`` with the two generators ``a`` and ``b``.
    The case which considers the vanishing of ``a`` is omitted: as a
    nonzero element of the base field, ``a`` cannot vanish.::

        sage: L = R.RosenfeldGroebner ([diff(y(x),x)^2-a*y(x)], basefield = F)
        sage: [ C.equations () for C in L ]
        [[-a*y(x) + diff(y(x), x)^2], [y(x)]]

    The next example illustrates field extensions defined by generators
    and relations. Denote ``f(x,y)`` = `\tan(x+y^2)`. Then ``f(x,y)`` satisfies
    the equations of the regular differential chain ``C``.
    The differential field ``F`` is generated (in the differential sense)
    by ``x``, ``y`` and ``f(x,y)``.  The generator ``f(x,y)`` satisfies the
    relations of ``C``::

        sage: from sage.calculus.DifferentialAlgebra import DifferentialRing, RegularDifferentialChain, BaseFieldExtension
        sage: x,y = var ('x,y')
        sage: f,u = function ('f,u')
        sage: R = DifferentialRing (derivations = [x,y], blocks = [u,f], parameters = [u(x)])
        sage: eqns = [diff(f(x,y),x) == f(x,y)^2 + 1, diff(f(x,y),y) == 2*y*(f(x,y)^2 + 1)]
        sage: C = RegularDifferentialChain (eqns, R, pretend=false)
        sage: F = BaseFieldExtension (relations = C)
        sage: F.generators ()
        [f(x, y)]
        sage: F.relations ().equations ()
        [-2*y*f(x, y)^2 - 2*y + diff(f(x, y), y), -f(x, y)^2 + diff(f(x, y), x) - 1]

    Again, the base field ``F`` permits to avoid some cases in the output of
    the :meth:`~sage.calculus.DifferentialAlgebra.DifferentialRing.RosenfeldGroebner` 
    algorithm. In the next example (``F`` is not used), two cases are computed::

        sage: eqn = diff(u(x),x)^2-(f(x,y)^2+1)*u(x) == 0
        sage: L = R.RosenfeldGroebner ([eqns[0], eqns[1], eqn])
        sage: [ C.equations () for C in L ]
        [[f(x, y)^2 + 1, diff(u(x), x)],
         [-2*y*f(x, y)^2 - 2*y + diff(f(x, y), y),
          -f(x, y)^2 + diff(f(x, y), x) - 1,
          u(x)]]

    In the next variant, a single case is computed. 
    The :meth:`~sage.calculus.DifferentialAlgebra.DifferentialRing.RosenfeldGroebner`
    algorithm is run over ``F``. As a nonzero element of ``F``, the polynomial
    `f(x, y)^2 + 1` cannot vanish. The only solution of the
    equations is thus `u(x) = 0`::

        sage: L = R.RosenfeldGroebner ([eqn], basefield = F)
        sage: [ C.equations () for C in L ]
        [[-2*y*f(x, y)^2 - 2*y + diff(f(x, y), y),
          -f(x, y)^2 + diff(f(x, y), x) - 1,
          u(x)]]
    """

    cdef bool rels_are_provided
    cdef list gens
    cdef RegularDifferentialChain rels
    """ rels_are_provided is set to true if relations were provided """
# __init__
    def __init__ (
            self,
            list generators = [],
            RegularDifferentialChain relations = None,
            DifferentialRing ring = None):
        r"""
        The constructor

        INPUT:

        - ``generators`` -- (default: []) a list of dependent variables
        - ``relations``  -- (optional) a regular differential chain
        - ``ring``       -- (optional) a differential ring

        OUTPUT:

        An extension F of the base field K of a differential
        ring, defined by ``generators`` and ``relations``.
        If ``ring`` is not specified, F actually is the
        description of a base field extension, which will be
        performed by a future call to the
        :meth:`~sage.calculus.DifferentialAlgebra.DifferentialRing.RosenfeldGroebner`
        algorithm.

        If ``relations`` is not specified, then F is the
        differential field, obtained by extending the field of
        the rational numbers, by the elements of ``generators``.

        If ``relations`` is specified, then F is the differential
        field of fractions of Q{X}/I where X denotes the variables
        (dependent and independent) present in ``generators``,
        ``relations``, and I denotes the differential ideal (which
        needs to be prime) defined by ``relations``.

        NOTES:

        The variables present in ``generators`` and ``relations``
        should belong to the lowest blocks of the ranking of the
        differential ring, whose base field is extended.

        EXAMPLES:

        The field $F$ of the rational numbers::

            sage: from sage.calculus.DifferentialAlgebra import DifferentialRing, RegularDifferentialChain, BaseFieldExtension
            sage: F = BaseFieldExtension ()

        The field `F = Q(a,b)`::

            sage: a,b = var ('a,b')
            sage: F = BaseFieldExtension (generators = [a,b])

        The field `\QQ(a,b)` again. The differential ring `R` is provided.
        This permits to check that `[a,b]` is the lowest block of the
        ranking of `R`::

            sage: y = function ('y')
            sage: R = DifferentialRing (derivations = [], blocks = [y,[a,b]])
            sage: F = BaseFieldExtension (generators = [a,b], ring = R)

        The field of fractions of `\QQ[a,b]/(a^2 + 1)`::

            sage: C = RegularDifferentialChain ([a^2 + 1], R)
            sage: F = BaseFieldExtension (generators = [a,b], relations = C)
        """
        cdef DifferentialRing R
        cdef bytes strgens, strrels, mesgerr
        cdef bmi_c.ALGEB_string A
        if (relations == None and ring == None):
            R = DifferentialRing (derivations = [], blocks = generators)
            self.rels_are_provided = False
            self.rels = RegularDifferentialChain ([], R)
            self.gens = generators
        elif relations != None:
            if not isinstance (relations, RegularDifferentialChain):
                raise TypeError, 'RegularDifferentialChain expected for relations keyword parameter'
            self.rels_are_provided = True
            self.rels = relations
            self.gens = generators
            strrels = self.__relations ()
            strgens = bytes (generators)
            strrels = relations.dring._translate_str(strrels)
            strgens = relations.dring._translate_str(strgens)
            A = bmi_c.bmi_sage_base_field_generators (
                        strgens, strrels, self.rels.regchain,
                        BMI_IX_undefined, BMI_IX_undefined, 0, 0)
            if bmi_c.bmi_sage_is_error (A):
                mesgerr = bmi_c.bmi_sage_mesgerr (A)
                bmi_c.bmi_balsa_clear_ALGEB (A)
                raise RuntimeError, mesgerr
            self.gens = relations.dring._eval_sage(A.value)
        elif ring != None:
            if not isinstance (ring, DifferentialRing):
                raise TypeError, 'DifferentialRing expected for ring keyword parameter'
            R = ring
            self.rels_are_provided = False
            self.rels = RegularDifferentialChain ([], R)
            self.gens = generators
            strrels = self.__relations ()
            strgens = bytes (generators)
            strrels = R._translate_str(strrels)
            strgens = R._translate_str(strgens)
            A = bmi_c.bmi_sage_base_field_generators (
                        strgens, strrels, R.dring,
                        BMI_IX_undefined, BMI_IX_undefined, 0, 0)
            if bmi_c.bmi_sage_is_error (A):
                mesgerr = bmi_c.bmi_sage_mesgerr (A)
                bmi_c.bmi_balsa_clear_ALGEB (A)
                raise RuntimeError, mesgerr
            self.gens = R._eval_sage(A.value)

    def __repr__ (self):
        """ The external representation """
        return 'differential_field'

    def _latex_ (self):
        """ The LaTeX representation """
        return 'differential\\_field'

    def generators (self):
        r"""
        The list of generators of a base field extension

        INPUT:

        Nothing

        OUTPUT:

        The list of generators

        EXAMPLES::

            sage: from sage.calculus.DifferentialAlgebra import DifferentialRing, RegularDifferentialChain, BaseFieldExtension
            sage: a,b = var ('a,b')
            sage: R = DifferentialRing (derivations = [], blocks = [a,b])
            sage: C = RegularDifferentialChain ([a^2 + 1], R)
            sage: F = BaseFieldExtension (generators = [a,b], relations = C)
            sage: F.generators ()
            [a, b]
        """
        return self.gens

    def relations (self):
        """
        The regular differential chain defining the relations of
        a base field extension

        INPUT:

        Nothing

        OUTPUT:

        A regular differential chain. If no relations were provided
        at the field construction, the regular differential chain
        is the zero chain.

        EXAMPLES::

            sage: from sage.calculus.DifferentialAlgebra import DifferentialRing, RegularDifferentialChain, BaseFieldExtension
            sage: a,b = var ('a,b')
            sage: R = DifferentialRing (derivations = [], blocks = [a,b])
            sage: C = RegularDifferentialChain ([a^2 + 1], R)
            sage: F = BaseFieldExtension (generators = [a,b], relations = C)
            sage: F.relations ()
            Regular chain (a^2 + 1)
            sage: F.relations ().equations ()
            [a^2 + 1]
        """
        return self.rels


    def __relations (self):
        r"""

        The relations, as a string, in a format that can be parsed by
        BMI after processing by DifferentialRing._translate_str()

        INPUT:

        Nothing

        OUTPUT:

        A string

        EXAMPLES::

            sage: from sage.calculus.DifferentialAlgebra import DifferentialRing, RegularDifferentialChain, BaseFieldExtension
            sage: a,b = var ('a,b')
            sage: R = DifferentialRing (derivations = [], blocks = [a,b])
            sage: C = RegularDifferentialChain ([a^2 + 1], R)
            sage: F = BaseFieldExtension (generators = [a,b], relations = C)
            sage: F.__relations ()
            'regchain ([a^2 + 1], [autoreduced, primitive, squarefree, normalized])'
        """
        cdef bytes strrels, strattr
        if self.rels_are_provided:
            strrels = bytes (self.rels.equations ())
            strattr = bytes (self.rels.attributes ())
            strattr = string.replace (bytes (strattr), "'", "")
            return 'regchain (%s, %s)' % (strrels, strattr)
        else:
            return 'regchain ([], [])'

