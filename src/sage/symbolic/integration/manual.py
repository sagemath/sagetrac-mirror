"""
Manual Symbolic Integration for educational use
"""

# ****************************************************************************
#       Copyright (C) 2019 Miguel Marco <mmarco@unizar.es>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************`

from sage.structure.sage_object import SageObject
from sage.symbolic.integration.integral import integral
from sympy.integrals.manualintegrate import (integral_steps, AlternativeRule,
                                    _manualintegrate,
                                    ConstantTimesRule,
                                    RewriteRule,
                                    AddRule,
                                    ReciprocalRule,
                                    PowerRule,
                                    URule,
                                    ExpRule,
                                    PartsRule,
                                    TrigRule,
                                    CyclicPartsRule)
from sage.symbolic.ring import SR
from sage.symbolic.expression import Expression
from operator import eq
from sage.functions.log import exp, log
from sage.misc.flatten import flatten

class ManualIntegral(SageObject):
    r"""
    A class to model the process of manual integrate an expression.
    It is created by passing an expression and a variable to integrate.

    We can then ask for hints to proceed with integration steps, either one by
    one or the whole step by step process from the beginning to the end.

    Hints can either help continue the process initiated by a previous hint, or
    provide an alternative approach.

    For example, let's try to compute the integral `\int x\sin(x)dx`. First
    we have to construct the corresponding manual integral object::

        sage: from sage.symbolic.integration.manual import ManualIntegral
        sage: M = ManualIntegral(x*sin(x), x)

    Now we can ask it for hints about how to proceed::

        sage: M.hint()
        (u == x,
         dv == sin(x),
         du == 1,
         v == integrate(sin(x), x),
         integrate(x*sin(x), x) == v*x - integrate(v, x))

    The first hint we are given is to use integration by parts. This reduces the
    problem of computing our original integral, to two subproblems: first we
    need to compute `v=\int\sin(x)dx` and once done that, we need to compute
    `\int vdx`. We can keep asking for hints. First it will give us hints to
    resolve the first subproblem, and once it has been solved, proceed with the
    second one::

        sage: M.hint()
        (integrate(sin(x), x) == -cos(x),)

    The computation of `v` was immediate, so it is solved with just one hint.
    Now  we ask for hints for our second subproblem::

        sage: M.hint()
        (integrate(-cos(x), x) == -integrate(cos(x), x),)
        sage: M.hint()
        (integrate(cos(x), x) == sin(x),)


    So we have to take the sign out of the integral, and then we just have an
    immediate integral. Since all the possible hints are given, if we ask for
    a next one, we get an error::

        sage: M.hint()
        Traceback (most recent call last):
        ...
        StopIteration

    But after that, we can start again with the first one::

        sage: M.hint()
        (u == x,
         dv == sin(x),
         du == 1,
         v == integrate(sin(x), x),
         integrate(x*sin(x), x) == v*x - integrate(v, x))

    We can also ask for all the hints at once::

        sage: M.all_hints()
        [(u == x,
          dv == sin(x),
          du == 1,
          v == integrate(sin(x), x),
          integrate(x*sin(x), x) == v*x - integrate(v, x)),
         (integrate(sin(x), x) == -cos(x),),
         (integrate(-cos(x), x) == -integrate(cos(x), x),),
         (integrate(cos(x), x) == sin(x),)]

    We have already seen examples of hints consisting on taking a sign out of
    the integrand, immediate integrals, and suggestions for integration by parts.
    There are also hints suggesting a change of variables::

        sage: N = ManualIntegral(cos(exp(x))*exp(x),x)
        sage: N.hint()
        (u == e^x, integrate(cos(e^x)*e^x, x) == integrate(cos(u), u))

    to split the integral as a sum::

        sage: L = ManualIntegral(2*sin(x) + 3*exp(x), x)
        sage: L.hint()
        (integrate(3*e^x + 2*sin(x), x) == integrate(3*e^x, x) + integrate(2*sin(x), x),)


    to take constants out of the integral::

        sage: L.hint()
        (integrate(3*e^x, x) == 3*integrate(e^x, x),)

    or just to rewrite the integrand::

        sage: O = ManualIntegral((x+2)/(x^2-2*x+1), x)
        sage: O.hint()
        (integrate((x + 2)/(x^2 - 2*x + 1), x) == integrate(1/(x - 1) + 3/(x - 1)^2, x),)

    When, after several application of the integration by parts, we get back
    again the same integral, we also get a hint telling us to isolate it::

        sage: P = ManualIntegral(sin(x)*exp(x), x)
        sage: P.hint()
        (u == sin(x),
         dv == e^x,
         du == cos(x),
         v == integrate(e^x, x),
         integrate(e^x*sin(x), x) == v*sin(x) - integrate(v*cos(x), x))
        sage: P.hint()
        (integrate(e^x, x) == e^x,)
        sage: P.hint()
        (u == cos(x),
         dv == e^x,
         du == -sin(x),
         v == integrate(e^x, x),
         integrate(cos(x)*e^x, x) == v*cos(x) - integrate(-v*sin(x), x))
        sage: P.hint()
        (integrate(e^x, x) == e^x,)
        sage: P.hint()
        (integrate(e^x*sin(x), x) == -cos(x)*e^x + e^x*sin(x) - integrate(e^x*sin(x), x),
         2*integrate(e^x*sin(x), x) == -cos(x)*e^x + e^x*sin(x))

    For some integrals, there might be different methods that can be suggested.
    In that case, we first get the full resolution using each one, and then start
    with the next one::

        sage: O.all_hints()
        [(integrate((x + 2)/(x^2 - 2*x + 1), x) == integrate(1/(x - 1) + 3/(x - 1)^2, x),),
         (integrate(1/(x - 1) + 3/(x - 1)^2, x) == integrate(1/(x - 1), x) + integrate(3/(x - 1)^2, x),),
         (u == x - 1, integrate(1/(x - 1), x) == integrate(1/u, u)),
         (integrate(1/u, u) == log(abs(u)),),
         (integrate(3/(x - 1)^2, x) == 3*integrate((x - 1)^(-2), x),),
         (u == x - 1, integrate((x - 1)^(-2), x) == integrate(u^(-2), u)),
         (integrate(u^(-2), u) == -1/u,),
         (integrate((x + 2)/(x^2 - 2*x + 1), x) == integrate(x/(x^2 - 2*x + 1) + 2/(x^2 - 2*x + 1), x),),
         (integrate(x/(x^2 - 2*x + 1) + 2/(x^2 - 2*x + 1), x) == integrate(x/(x^2 - 2*x + 1), x) + integrate(2/(x^2 - 2*x + 1), x),),
         (integrate(x/(x^2 - 2*x + 1), x) == integrate(1/(x - 1) + 1/(x - 1)^2, x),),
         (integrate(1/(x - 1) + 1/(x - 1)^2, x) == integrate(1/(x - 1), x) + integrate((x - 1)^(-2), x),),
         (u == x - 1, integrate(1/(x - 1), x) == integrate(1/u, u)),
         (integrate(1/u, u) == log(abs(u)),),
         (u == x - 1, integrate((x - 1)^(-2), x) == integrate(u^(-2), u)),
         (integrate(u^(-2), u) == -1/u,),
         (integrate(2/(x^2 - 2*x + 1), x) == 2*integrate(1/(x^2 - 2*x + 1), x),),
         (integrate(1/(x^2 - 2*x + 1), x) == integrate((x - 1)^(-2), x),),
         (u == x - 1, integrate((x - 1)^(-2), x) == integrate(u^(-2), u)),
         (integrate(u^(-2), u) == -1/u,)]

    """
    def __init__(self, f, x):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: from sage.symbolic.integration.manual import ManualIntegral
            sage: M = ManualIntegral(x*sin(x), x)
            sage: M
            integrate(x*sin(x), x)
            sage: TestSuite(M).run()

        """
        self._f = f
        self._x = x
        self._rule = integral_steps(f._sympy_(), x._sympy_())
        self._rules_list = create_rules_list(self._rule)
        self._current_index = -1

    def __iter__(self):
        r"""
        The iterator of an iterable object should be itself

        EXAMPLES::

            sage: from sage.symbolic.integration.manual import ManualIntegral
            sage: M = ManualIntegral(x*sin(x), x)
            sage: M
            integrate(x*sin(x), x)
            sage: M.__iter__()
            integrate(x*sin(x), x)

        """
        return self

    def __eq__(self, other):
        r"""
        Two objects of this class are equal if they compute the integral of equal
        expressions with respect to equal variables.

        EXAMPLES::

            sage: from sage.symbolic.integration.manual import ManualIntegral
            sage: M = ManualIntegral(x*sin(x), x)
            sage: N = ManualIntegral(cos(exp(x))*exp(x), x)
            sage: O = ManualIntegral(x*sin(x), x)
            sage: M == N
            False
            sage: M == O
            True
            sage: M.__eq__(N)
            False
            sage: M.__eq__(O)
            True

        """
        if not type(self) == type(other):
            return False
        return (self._f, self._x) == (other._f, other._x)

    def _repr_(self):
        r"""
        Objects of this class are represented as the integral they compute.

        EXAMPLES::

            sage: from sage.symbolic.integration.manual import ManualIntegral
            sage: M = ManualIntegral(x*sin(x), x)
            sage: M._repr_()
            'integrate(x*sin(x), x)'

        """
        return integral(self._f, self._x, hold=True)._repr_()

    def _latex_(self):
        r"""
        The latex expression of an object of this class is the same as the
        integral it computes.

        EXAMPLES::

            sage: from sage.symbolic.integration.manual import ManualIntegral
            sage: M = ManualIntegral(x*sin(x), x)
            sage: M._latex_()
            '\\int x \\sin\\left(x\\right)\\,{d x}'

        """
        return integral(self._f, self._x, hold=True)._latex_()

    def _ascii_art_(self):
        r"""
        The ascii_art expression of an object of this class is the same as the
        integral it computes.

        EXAMPLES::

            sage: from sage.symbolic.integration.manual import ManualIntegral
            sage: M = ManualIntegral(x*sin(x), x)
            sage: M._ascii_art_()
              /
             |
             | x*sin(x) dx
             |
            /
        """
        return integral(self._f, self._x, hold=True)._ascii_art_()

    def _unicode_art_(self):
        r"""
        The unicode_art expression of an object of this class is the same as the
        integral it computes.

        EXAMPLES::

            sage: from sage.symbolic.integration.manual import ManualIntegral
            sage: M = ManualIntegral(x*sin(x), x)
            sage: M._unicode_art_()
            ⌠
            ⎮ x⋅sin(x) dx
            ⌡
        """
        return integral(self._f, self._x, hold=True)._unicode_art_()

    def next(self):
        r"""
        A ManualIntegral can be iterated over its hints. Each time the method
        ``hint`` or ``next`` is called, a new hint is returned.

        EXAMPLES::

            sage: from sage.symbolic.integration.manual import ManualIntegral
            sage: M = ManualIntegral(x**2+3*x, x)
            sage: M.hint()
            (integrate(x^2 + 3*x, x) == integrate(x^2, x) + integrate(3*x, x),)
            sage: M.hint()
            (integrate(x^2, x) == 1/3*x^3,)

        The methods ``next`` or ``hint``can be used indistinctly::

            sage: M.next()
            (integrate(3*x, x) == 3*integrate(x, x),)

        When all the hints have been returned, and a new hint is requested,
        an exception is raised::

            sage: M.next()
            (integrate(x, x) == 1/2*x^2,)
            sage: M.next()
            Traceback (most recent call last):
            ...
            StopIteration

        After that, the process can start again::

            sage: M.hint()
            (integrate(x^2 + 3*x, x) == integrate(x^2, x) + integrate(3*x, x),)

        """
        self._current_index += 1
        if self._current_index >= len(self._rules_list):
            self._current_index = -1
            raise StopIteration
        else:
            return integration_hint(self._rules_list[self._current_index])

    hint = next

    def all_hints(self):
        r"""
        Return a list with all the hints to compute this integral.

        EXAMPLES::

            sage: from sage.symbolic.integration.manual import ManualIntegral
            sage: M = ManualIntegral(x*sin(x), x)
            sage: M.all_hints()
            [(u == x,
              dv == sin(x),
              du == 1,
              v == integrate(sin(x), x),
              integrate(x*sin(x), x) == v*x - integrate(v, x)),
             (integrate(sin(x), x) == -cos(x),),
             (integrate(-cos(x), x) == -integrate(cos(x), x),),
             (integrate(cos(x), x) == sin(x),)]

        """
        return [integration_hint(r) for r in self._rules_list]


def create_rules_list(rule):
    r"""
    Return a list with rules for all the substeps of a rule for integration.

    INPUT:

    - `rule` -- a sympy `manualintegrate` Rule

    OUTPUT: a list of sympy rules, in the appropriate order.

    EXAMPLES:

    In the case of `CyclicPartsRule`, the different integration by parts steps
    are put before the actual rule for collecting the cyclic parts::

        sage: from sage.symbolic.integration.manual import ManualIntegral, create_rules_list
        sage: M = ManualIntegral(exp(x)*sin(x), x)
        sage: M._rule
        CyclicPartsRule(parts_rules=[PartsRule(u=sin(x), dv=exp(x), v_step=ExpRule(base=E, exp=x, context=exp(x), symbol=x), second_step=None, context=None, symbol=None), PartsRule(u=cos(x), dv=exp(x), v_step=ExpRule(base=E, exp=x, context=exp(x), symbol=x), second_step=None, context=None, symbol=None)], coefficient=-1, context=exp(x)*sin(x), symbol=x)
        sage: create_rules_list(M._rule)
        [PartsRule(u=sin(x), dv=exp(x), v_step=ExpRule(base=E, exp=x, context=exp(x), symbol=x), second_step=None, context=exp(x)*sin(x), symbol=x),
         ExpRule(base=E, exp=x, context=exp(x), symbol=x),
         PartsRule(u=cos(x), dv=exp(x), v_step=ExpRule(base=E, exp=x, context=exp(x), symbol=x), second_step=None, context=exp(x)*cos(x), symbol=x),
         ExpRule(base=E, exp=x, context=exp(x), symbol=x),
         CyclicPartsRule(parts_rules=[PartsRule(u=sin(x), dv=exp(x), v_step=ExpRule(base=E, exp=x, context=exp(x), symbol=x), second_step=None, context=None, symbol=None), PartsRule(u=cos(x), dv=exp(x), v_step=ExpRule(base=E, exp=x, context=exp(x), symbol=x), second_step=None, context=None, symbol=None)], coefficient=-1, context=exp(x)*sin(x), symbol=x)]

    For `AddRule` the different subrules are put after it::

        sage: N = ManualIntegral(2*x + 1/x, x)
        sage: N._rule
        AddRule(substeps=[ConstantTimesRule(constant=2, other=x, substep=PowerRule(base=x, exp=1, context=x, symbol=x), context=2*x, symbol=x), ReciprocalRule(func=x, context=1/x, symbol=x)], context=2*x + 1/x, symbol=x)
        sage: create_rules_list(N._rule)
        [AddRule(substeps=[ConstantTimesRule(constant=2, other=x, substep=PowerRule(base=x, exp=1, context=x, symbol=x), context=2*x, symbol=x), ReciprocalRule(func=x, context=1/x, symbol=x)], context=2*x + 1/x, symbol=x),
         ConstantTimesRule(constant=2, other=x, substep=PowerRule(base=x, exp=1, context=x, symbol=x), context=2*x, symbol=x),
         PowerRule(base=x, exp=1, context=x, symbol=x),
         ReciprocalRule(func=x, context=1/x, symbol=x)]

    For `AlternativeRule`, the different alternatives are put one after another::

        sage: O = ManualIntegral((x+1)/(x**2+2*x+2),x)
        sage: O._rule
        AlternativeRule(alternatives=[URule(u_var=_u, u_func=x**2 + 2*x + 2, constant=1/2, substep=ConstantTimesRule(constant=1/2, other=1/_u, substep=ReciprocalRule(func=_u, context=1/_u, symbol=_u), context=1/_u, symbol=_u), context=(x + 1)/(x**2 + 2*x + 2), symbol=x)], context=(x + 1)/(x**2 + 2*x + 2), symbol=x)
        sage: create_rules_list(O._rule)
        [URule(u_var=_u, u_func=x**2 + 2*x + 2, constant=1/2, substep=ConstantTimesRule(constant=1/2, other=1/_u, substep=ReciprocalRule(func=_u, context=1/_u, symbol=_u), context=1/_u, symbol=_u), context=(x + 1)/(x**2 + 2*x + 2), symbol=x),
         ConstantTimesRule(constant=1/2, other=1/_u, substep=ReciprocalRule(func=_u, context=1/_u, symbol=_u), context=1/_u, symbol=_u),
         ReciprocalRule(func=_u, context=1/_u, symbol=_u)]

    If the rule has a substep, the corresponding rules are put after the rule itself::

        sage: P = ManualIntegral(cos(x)^2,x)
        sage: P._rule
        RewriteRule(rewritten=cos(2*x)/2 + 1/2, substep=AddRule(substeps=[ConstantTimesRule(constant=1/2, other=cos(2*x), substep=URule(u_var=_u, u_func=2*x, constant=1/2, substep=ConstantTimesRule(constant=1/2, other=cos(_u), substep=TrigRule(func='cos', arg=_u, context=cos(_u), symbol=_u), context=cos(_u), symbol=_u), context=cos(2*x), symbol=x), context=cos(2*x)/2, symbol=x), ConstantRule(constant=1/2, context=1/2, symbol=x)], context=cos(2*x)/2 + 1/2, symbol=x), context=cos(x)**2, symbol=x)
        sage: create_rules_list(P._rule)
        [RewriteRule(rewritten=cos(2*x)/2 + 1/2, substep=AddRule(substeps=[ConstantTimesRule(constant=1/2, other=cos(2*x), substep=URule(u_var=_u, u_func=2*x, constant=1/2, substep=ConstantTimesRule(constant=1/2, other=cos(_u), substep=TrigRule(func='cos', arg=_u, context=cos(_u), symbol=_u), context=cos(_u), symbol=_u), context=cos(2*x), symbol=x), context=cos(2*x)/2, symbol=x), ConstantRule(constant=1/2, context=1/2, symbol=x)], context=cos(2*x)/2 + 1/2, symbol=x), context=cos(x)**2, symbol=x),
         AddRule(substeps=[ConstantTimesRule(constant=1/2, other=cos(2*x), substep=URule(u_var=_u, u_func=2*x, constant=1/2, substep=ConstantTimesRule(constant=1/2, other=cos(_u), substep=TrigRule(func='cos', arg=_u, context=cos(_u), symbol=_u), context=cos(_u), symbol=_u), context=cos(2*x), symbol=x), context=cos(2*x)/2, symbol=x), ConstantRule(constant=1/2, context=1/2, symbol=x)], context=cos(2*x)/2 + 1/2, symbol=x),
         ConstantTimesRule(constant=1/2, other=cos(2*x), substep=URule(u_var=_u, u_func=2*x, constant=1/2, substep=ConstantTimesRule(constant=1/2, other=cos(_u), substep=TrigRule(func='cos', arg=_u, context=cos(_u), symbol=_u), context=cos(_u), symbol=_u), context=cos(2*x), symbol=x), context=cos(2*x)/2, symbol=x),
         URule(u_var=_u, u_func=2*x, constant=1/2, substep=ConstantTimesRule(constant=1/2, other=cos(_u), substep=TrigRule(func='cos', arg=_u, context=cos(_u), symbol=_u), context=cos(_u), symbol=_u), context=cos(2*x), symbol=x),
         ConstantTimesRule(constant=1/2, other=cos(_u), substep=TrigRule(func='cos', arg=_u, context=cos(_u), symbol=_u), context=cos(_u), symbol=_u),
         TrigRule(func='cos', arg=_u, context=cos(_u), symbol=_u),
         ConstantRule(constant=1/2, context=1/2, symbol=x)]

    """
    if not rule:
        return []
    elif isinstance(rule, AlternativeRule):
        return flatten([create_rules_list(r) for r in rule.alternatives], max_level=1)
    elif isinstance(rule, PartsRule):
        res = [rule]
        res += create_rules_list(rule.v_step)
        res += create_rules_list(rule.second_step)
        return res
    elif isinstance(rule, CyclicPartsRule):
        auxrules = []
        context = rule.context
        symbol = rule.symbol
        for r in rule.parts_rules:
            u = r.u
            dv = r.dv
            vs = r.v_step
            ss = r.second_step
            nr = PartsRule(u=u, dv=dv, context=u*dv, symbol=symbol, v_step=vs, second_step=ss)
            auxrules += create_rules_list(nr)
        return auxrules + [rule]
    elif isinstance(rule, AddRule):
        return [rule] + flatten([create_rules_list(r) for r in rule.substeps], max_level=1)
    elif hasattr(rule, 'substep'):
        return [rule] + create_rules_list(rule.substep)
    else:
        return [rule]


def integration_hint(rule):
    r"""
    Construct a tuple of equations giving the integration hint for a certain rule.

    INPUT:

    - ``rule`` -- A sympy manualintegrate rule

    OUTPUT:

    - A tuple of equations, with the corresponding hint.

    EXAMPLES::

        sage: from sympy import S
        sage: from sympy.integrals.manualintegrate import integral_steps
        sage: from sage.symbolic.integration.manual import integration_hint
        sage: f = exp(x)*x
        sage: fs = S(f)
        sage: xs = S(x)
        sage: rule = integral_steps(fs, xs)
        sage: rule
        PartsRule(u=x, dv=exp(x), v_step=ExpRule(base=E, exp=x, context=exp(x), symbol=x), second_step=ExpRule(base=E, exp=x, context=exp(x), symbol=x), context=x*exp(x), symbol=x)
        sage: integration_hint(rule)
        (u == x,
         dv == e^x,
         du == 1,
         v == integrate(e^x, x),
         integrate(x*e^x, x) == v*x - integrate(v, x))


    TESTS::

        sage: from sympy import S
        sage: from sage.symbolic.integration.manual import integration_hint
        sage: from sympy.integrals.manualintegrate import integral_steps
        sage: integration_hint(integral_steps(S(5*x), S(x)))
        (integrate(5*x, x) == 5*integrate(x, x),)
        sage: integration_hint(integral_steps(S(5 + x), S(x)))
        (integrate(x + 5, x) == integrate(5, x) + integrate(x, x),)
        sage: integration_hint(integral_steps(S(log(x)/x), S(x)).alternatives[0])
        (u == (1/x), integrate(log(x)/x, x) == integrate(-log(1/u)/u, u))
        sage: integration_hint(integral_steps(S(cos(x)), S(x)))
        (integrate(cos(x), x) == sin(x),)
        sage: integration_hint(integral_steps(S(exp(x)), S(x)))
        (integrate(e^x, x) == e^x,)
        sage: integration_hint(integral_steps(S(2^(x)), S(x)))
        (integrate(2^x, x) == e^x/log(2),)
        sage: integration_hint(integral_steps(S(exp(x)*x), S(x)))
        (u == x,
         dv == e^x,
         du == 1,
         v == integrate(e^x, x),
         integrate(x*e^x, x) == v*x - integrate(v, x))
        sage: integration_hint(integral_steps(S(exp(x)*cos(x)), S(x)))
        (integrate(cos(x)*e^x, x) == cos(x)*e^x + e^x*sin(x) - integrate(cos(x)*e^x, x),
         2*integrate(cos(x)*e^x, x) == cos(x)*e^x + e^x*sin(x))
    """
    f = SR(rule.context)
    x = SR(rule.symbol)
    lhs = integral(f, x, hold=True)
    if isinstance(rule, ConstantTimesRule):
        constant = SR(rule.constant)
        other = SR(rule.other)
        rhs = constant * integral(other, x, hold = True)
        lhs = integral(constant * other, x, hold=True)
        return (lhs == rhs,)
    elif isinstance(rule, RewriteRule):
        rewritten = SR(rule.rewritten)
        rhs = integral(rewritten, x, hold=True)
        return (lhs == rhs,)
    elif isinstance(rule, AddRule):
        steps = rule.substeps
        rhs = sum(integral(SR(s.context), x, hold=True) for s in steps)
        return (lhs == rhs,)
    elif isinstance(rule, ReciprocalRule):
        rhs = log(abs(SR(rule.func)))
        return (lhs == rhs,)
    elif isinstance(rule, PowerRule):
        rhs = SR(rule.base**(rule.exp+1)/(rule.exp+1))
        return (lhs == rhs,)
    elif isinstance(rule, TrigRule):
        result = SR(_manualintegrate(rule))
        return (lhs == result,)
    elif isinstance(rule, ExpRule):
        base = SR(rule.base)
        exponent = SR(rule.exp)
        return (lhs == exp(exponent)/log(base),)
    elif isinstance(rule, URule):
        u = SR(rule.u_var)
        uf = SR(rule.u_func)
        if rule.constant is None:
            const = SR(1)
        else:
            const = SR(rule.constant)
        subcontext = SR(rule.substep.context)*const
        if u == x:
            U = SR.var('U')
            subcontext = subcontext.subs({u:U})
            u = U
        return (u == uf, lhs == integral(subcontext, u, hold=True))
    elif type(rule) == PartsRule:
        u = SR(rule.u)
        dv = SR(rule.dv)
        du = SR(u.diff(x))
        v = integral(dv, x, hold=True)
        uv = SR.var('u')
        dvv = SR.var('dv')
        vv = SR.var('v')
        duv = SR.var('du')
        if str(x) == 'u':
            uv = SR.var('U')
            vv = SR.var('V')
            dvv = SR.var('dV')
            duv = SR.var('dU')
        return (uv == u, dvv == dv, duv == du, vv == v, lhs == u*vv - integral(vv * du, x, hold=True))
    elif type(rule) == CyclicPartsRule:
        rhs = 0
        coefficient = 1
        for r in rule.parts_rules:
            v = _manualintegrate(r.v_step)
            rhs += coefficient * r.u * v
            coefficient *= -1
        res0 = SR(rhs) + SR(rule.coefficient) * integral(f, x, hold=True)
        return ( lhs == res0, SR(1-rule.coefficient)*lhs ==  rhs)
