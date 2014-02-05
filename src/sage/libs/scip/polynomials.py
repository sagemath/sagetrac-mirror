"""
(Boolean) polynomial system solving with SCIP.

AUTHOR:

- Martin Albrecht (2010-13, initial version)
"""

################################################################################
#    Copyright (C) 2010-2013 Martin Albrecht <martinralbrecht@googlemail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
################################################################################

from sage.libs.scip.scip import SCIP
from sage.misc.updown import UpDown
from sage.numerical.mip import MIPSolverException

def default_objective_weight_callback(var):
    """
    The default objective weight for any variable is 1.0

    EXAMPLE::

        sage: from sage.libs.scip.polynomials import default_objective_weight_callback
        sage: default_objective_weight_callback(None)
        1.0
    """
    return 1.0

def solve(F, maximization=True, objective_weight_callback=None, parameters=None):
    """
    Solve (Boolean) polynomial system ``F`` and solve using SCIP
    maximizing/minimizing the cost function given by
    ``objective_weight_callback``.

    INPUT:

    - ``F`` - an instance of :class:`PolynomialSequence
      <sage.rings.polynomial.multi_polynomial_sequence>`

    - ``maximization`` - if `True` the objective sense is to maximize (default:
      ``True``)

    - ``objective_weight_callback`` - called on each variable and must return a
      floating point number which is the coefficient of this variable in the
      objective function. If ``None`` is given, the coefficient is 1.0 for all
      variables (default: ``None``)

    - ``parameters`` - passed as parameters to the SCIP solver.

    OUTPUT: A dictionary with a solution for ``F`` or ``False``.

    EXAMPLE::

        sage: B.<a,b,c> = BooleanPolynomialRing()
        sage: F = Sequence([a*b + c + 1, c+1]); F
        [a*b + c + 1, c + 1]

        sage: from sage.libs.scip.polynomials import solve as scip_solve # optional - SCIP
        sage: scip_solve(F)                                              # optional - SCIP
        {b: 0, c: 1, a: 1}
        sage: scip_solve(F,maximization=False)                           # optional - SCIP
        {b: 0, c: 1, a: 0}

    ::

        sage: from sage.libs.scip.polynomials import solve as scip_solve # optional - SCIP
        sage: P = PolynomialRing(GF(7),5,'x')
        sage: I = sage.rings.ideal.Katsura(P)
        sage: I.variety()
        [{x2: 0, x1: 0, x0: 1, x4: 0, x3: 0},
         {x2: 0, x1: 0, x0: 5, x4: 5, x3: 0},
         {x2: 3, x1: 0, x0: 6, x4: 5, x3: 0}]
        sage: F = I.gens()

        sage: scip_solve(F)                    # optional - SCIP
        {x2: 3, x1: 0, x0: 6, x4: 5, x3: 0}
        sage: scip_solve(F,maximization=False) # optional - SCIP
        {x2: 0, x1: 0, x0: 1, x4: 0, x3: 0}
    """
    scip = SCIP(maximization=maximization)

    if F.ring().base_ring().order() == 2:
        Q = boolean_polynomials(scip, F, use_xor=True, objective_weight_callback=objective_weight_callback)
    else:
        Q = polynomials(scip, F, objective_weight_callback=objective_weight_callback)

    if parameters:
        for k,v in parameters.iteritems():
            scip[k] = v

    K = F.ring().base_ring()

    try:
        scip.solve()
        s = dict()
        for v in F.ring().gens():
            s[v] = K(int(scip.get_variable_value(Q.down(v))))
        return s
    except MIPSolverException:
        return False

def boolean_polynomials(scip, F, use_xor=True, objective_weight_callback=None):
    """
    Read Boolean polynomial system ``F`` into SCIP instance ``scip``.

    INPUT:

    - ``F`` - an instance of :class:`PolynomialSequence
      <sage.rings.polynomial.multi_polynomial_sequence>`

    - ``use_xor`` - if ``False`` linear constraints are used to model
      logical-XOR chains, othwerweise XOR constraints are used (default:
      ``True``)

    - ``objective_weight_callback`` - called on each variable and must return a
      floating point number which is the coefficient of this variable in the
      objective function. If ``None`` is given, the coefficient is 1.0 for all
      variables (default: ``None``)

    OUTPUT: A dictionary mapping polynomial ring variables to SCIP variables and back.
    """
    V = sorted(F.variables())
    M = sorted(F.monomials())

    if objective_weight_callback is None:
        objective_weight_callback = default_objective_weight_callback

    Q = UpDown()
    Q1 = UpDown()

    for v in V:
        obj = objective_weight_callback(v)
        v1 = scip.add_variable(name=str(v), binary=True, obj=obj)
        Q[v] = v1
        Q1[v] = v1

    for m in M:
        if m.degree() <= 1:
            continue

        m1 = scip.add_variable(name=str(m), binary=True, obj=0.0)
        scip.add_and_constraint(m1, [Q.down(v) for v in m.variables()])
        Q[m] = m1

    if use_xor:
        for i,f in enumerate(F):
            cc = f.constant_coefficient()
            f1 = [Q.down(m) for m in (f-cc).monomials()]
            scip.add_xor_constraint(float(cc), f1)
    else:
        for i,f in enumerate(F):
            cc = f.constant_coefficient()
            ub = len(f.monomials())
            ub = (ub - ub%2)
            m1 = scip.add_variable(name="xor%04d"%i, integer=True, obj=0.0, lower_bound=0.0, upper_bound=ub)
            coeffs  = [1.0 for m in (f-cc).monomials()] + [-2.0]
            indices = [Q.down(m) for m in (f-cc).monomials()] + [m1]
            scip.add_linear_constraint(zip(indices, coeffs), cc, cc)

    return Q1

def polynomials(scip, F, objective_weight_callback=None):
    """
    Read polynomial system ``F`` into SCIP instance ``scip``.

    INPUT:

    - ``F`` - an instance of :class:`PolynomialSequence
      <sage.rings.polynomial.multi_polynomial_sequence>`

    - ``objective_weight_callback`` - called on each variable and must return a
      floating point number which is the coefficient of this variable in the
      objective function. If ``None`` is given, the coefficient is 1.0 for all
      variables (default: ``None``)

    OUTPUT: A dictionary mapping polynomial ring variables to SCIP variables and back.
    """
    R = F.ring()
    p = R.base_ring().order()
    assert(p.is_prime())

    if objective_weight_callback is None:
        objective_weight_callback = default_objective_weight_callback

    V = sorted(F.variables())

    Q = UpDown()

    assert(max([f.degree() for f in F]) <= 2)

    for v in V:
        obj = objective_weight_callback(v)
        v1 = scip.add_variable(name=str(v), integer=True, obj=obj, lower_bound=0, upper_bound=p-1)
        Q[v] = v1

    for i,f in enumerate(F):
        m1 = scip.add_variable(name="multp%d"%i, integer=True, obj=0.0, lower_bound=0)
        l,q = [(m1,-p)],[]
        cc = 0
        for c, m in f:
            if m.degree() == 2:
                v = m.variables()
                if len(v) == 1:
                    v = v+v
                q.append( tuple(map(Q.down,v)) + (c,) )
            elif m.degree() == 1:
                l.append( (Q.down(m),c) )
            elif m.degree() == 0:
                cc = c
        if q:
            scip.add_quadratic_constraint(l, q, -cc, -cc)
        elif l:
            scip.add_linear_constraint(l, -cc, -cc)
        else:
            pass
    return Q
