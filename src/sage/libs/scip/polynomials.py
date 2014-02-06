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
    Solve Boolean polynomial system ``F`` and solve using SCIP
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
    """
    scip = SCIP(maximization=maximization)

    if F.ring().base_ring().order() == 2:
        Q = boolean_polynomials(scip, F, use_xor=True, objective_weight_callback=objective_weight_callback)
    else:
        raise NotImplementedError

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

    - ``scip`` - a SCIP instance
    
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


    EXAMPLE:

    We find the minimal solution for an underdetermined system of boolean polynomials::

        sage: from sage.libs.scip.polynomials import boolean_polynomials # optional - SCIP
        sage: from sage.libs.scip import SCIP                            # optional - SCIP
        sage: scip = SCIP(maximization=False)                            # optional - SCIP
        sage: B.<a,b,c,d,e,f> = BooleanPolynomialRing()                  # optional - SCIP
        sage: F = Sequence([B.random_element() for _ in range(4)])       # optional - SCIP
        sage: updown = boolean_polynomials(scip, F)                      # optional - SCIP
        sage: scip.solve()                                               # optional - SCIP
        sage: dict([(x,scip.get_variable_value(updown.down(x))) for x in B.gens()]) # optional - SCIP
        {f: 1.0, e: 0.0, d: 0.0, b: 0.0, c: 1.0, a: 0.0}
        sage: F.ideal().variety()                                        # optional - SCIP
        [{f: 0, e: 0, d: 0, b: 1, c: 1, a: 1},
         {f: 0, e: 0, d: 1, b: 1, c: 1, a: 1},
         {f: 1, e: 0, d: 0, b: 0, c: 1, a: 0},
         {f: 1, e: 0, d: 1, b: 1, c: 1, a: 1}]
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
