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
    """
    def __init__(self, f, x):
        self._f = f
        self._x = x
        self._rule = integral_steps(f._sympy_(), x._sympy_())
        self._rules_list = create_rules_list(self._rule)
        self._current_index = -1

    def __iter__(self):
        return self

    def _repr_(self):
        return integral(self._f, self._x, hold=True)._repr_()

    def _latex_(self):
        return integral(self._f, self._x, hold=True)._latex_()

    def next(self):
        self._current_index += 1
        if self._current_index >= len(self._rules_list):
            self._current_index = -1
            raise StopIteration
        else:
            return integration_hint(self._rules_list[self._current_index])

    hint = next

    def all_hints(self):
        return [integration_hint(r) for r in self._rules_list]


def create_rules_list(rule):
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




def hint_to_integrate(*args):
    r"""
    Give possible hints to integrate a function

    The input must be either the symbolic expression to integrate, and the
    variable to integrate with respect to; or a symbolic expression that contains
    an integral in it. If the expression is an equality, the right hand side
    is used.

    The output is a list of hints. Each hint is a tuple of equations. If the
    hint consist on a substitution, the first equation gives the substitution
    to perform. If the hint is to integrate by parts, the first two equations
    give the ``u`` and ``dv``. The last equation contains the integral after the
    hint is applied.

    EXAMPLES::

        sage: from sage.symbolic.integration.manual import hint_to_integrate
        sage: hint_to_integrate((x^2+3*x+1)/(x^3-x),x)
        [(integrate((x^2 + 3*x + 1)/(x^3 - x), x) == integrate(-1/2/(x + 1) + 5/2/(x - 1) - 1/x, x),),
        (integrate((x^2 + 3*x + 1)/(x^3 - x), x) == integrate(x^2/(x^3 - x) + 3*x/(x^3 - x) + 1/(x^3 - x), x),)]
        sage: hint_to_integrate(log(x)/x, x)
        [(u == (1/x), integrate(log(x)/x, x) == integrate(-log(1/u)/u, u)),
        (u == log(x), integrate(log(x)/x, x) == integrate(u, u))]


    We can see how the hint to perform integration by parts is given::

        sage: h = hint_to_integrate(sin(x)*exp(x),x)
        sage: h
        [(u == sin(x),
        dv == e^x,
        integrate(e^x*sin(x), x) == e^x*sin(x) - integrate(cos(x)*e^x, x))]

    The output can be used to get a new hint for the next step to follow::

        sage: new_expression = h[0][2]
        sage: new_expression
        integrate(e^x*sin(x), x) == e^x*sin(x) - integrate(cos(x)*e^x, x)
        sage: hint_to_integrate(new_expression)
        [(u == cos(x),
        dv == e^x,
        integrate(cos(x)*e^x, x) == cos(x)*e^x - integrate(-e^x*sin(x), x))]
        sage: hint_to_integrate(new_expression.rhs())
        [(u == cos(x),
        dv == e^x,
        integrate(cos(x)*e^x, x) == cos(x)*e^x - integrate(-e^x*sin(x), x))]

    Sometimes the hint is as simple as separating summands, or taking constants
    out of the expression::

        sage: hint_to_integrate(cos(x) + x, x)
        [(integrate(x + cos(x), x) == integrate(x, x) + integrate(cos(x), x),)]
        sage: hint_to_integrate(5*exp(x), x)
        [(integrate(5*e^x, x) == 5*integrate(e^x, x),)]

    """
    if len(args) == 2:
        f = args[0]
        x = args[1]
    elif len(args) == 1:
        if not isinstance(args[0], Expression):
            raise ValueError("input must be a symbolic expression, maybe with a variable")
        w0 = SR.wild(0)
        w1 = SR.wild(1)
        if args[0].operator() == eq:
            expr = args[0].rhs()
        else:
            expr = args[0]
        if not expr.has(integral(w0, w1)):
            raise ValueError("no symbolic integral in input")
        integ = expr.find(integral(w0, w1))[0]
        f, x = integ.operands()
    else:
        raise ValueError("input must be a symbolic expression, maybe with a variable")
    rule = integral_steps(f._sympy_(), x._sympy_())
    if isinstance(rule, AlternativeRule):
        return [integration_hint(r) for r in rule.alternatives]
    return [integration_hint(rule)]

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


class IntegrationHint(SageObject):
    r"""
    Container class to represent the hints.
    """
    def __init__(self, rule):

        self._f = SR(rule.context)
        self._x = SR(rule.symbol)
        self._rule = rule
        self._lhs = integral(_f, _x, hold=True)

    def _repr_(self):
        rule = self._rule
        if isinstance(rule, ConstantTimesRule):
            constant = SR(rule.constant)
            other = SR(rule.other)
            rhs = constant * integral(other, self._x, hold = True)
            return (self._lhs == rhs)._repr_()
        if isinstance(rule, RewriteRule):
            rewritten = SR(rule.rewritten)
            rhs = integral(rewritten, self._x, hold=True)
            return (self._lhs == rhs)._repr_()
        if isinstance(rule, AddRule):
            steps = rule.substeps
            rhs = sum(integral(SR(s.context), self._x, hold=True) for s in steps)
            return (self._lhs == rhs)._repr_()
        if isinstance(rule, TrigRule):
            result = SR(_manualintegrate(rule))
            return (lhs == result)._repr_()
        if isintance(rule, ExpRule):
            base = SR(rule.base)
            exponent = SR(rule.exp)
            string = 'The integral of an exponential function is itself divided by the natural logarithm of the base, so: '
            string += (rhs == exp(exponent)/log(base))._repr_()
            return string

def rule_string(rule, depth):

    if rule is None:
        return ''
    f = SR(rule.context)
    x = SR(rule.symbol)
    if type(rule) == AlternativeRule:
        string = 'There are several ways to compute $\displaystyle\int {} d{}$<br>'.format(latex(f),latex(x))
        alternatives = []
        for r in rule.alternatives:
            if type(r) == AlternativeRule:
                alternatives += r.alternatives
            else:
                alternatives.append(r)
        for (i,srule) in enumerate(alternatives):
            string += '{})<br>'.format(i+1)
            string += rule_string(srule, depth)+'<br>'
        return string
    if type(rule) == ConstantTimesRule:
        return rule_string(rule.substep, depth-1)
    string = 'Hint: to compute $$\int {} d{}$$<br>'.format(latex(f), latex(x))
    if type(rule) == RewriteRule:
        rewritten = SR(rule.rewritten)
        string += '$$ \int {} d{} = \int {} d{}$$<br>'.format(latex(f),latex(x),latex(rewritten),latex(x))
        string += rule_string(rule.substep, depth-1)
        return string
    if type (rule) == AddRule:
        substeps = rule.substeps
        string += '$$\int {} d{} = '.format(latex(f),latex(x))
        for substep in substeps:
            string += '\int {} d{} +'.format(latex(SR(substep.context)),latex(x))
        string = string[:-1] + '$$<br>'
        finalstring = ''
        for substep in substeps:
            string += rule_string(substep, depth-1)
            finalstring += latex(SR(_manualintegrate(substep))) + '+'
        string += 'So $$\int {} d{} = '.format(latex(f),latex(x)) + finalstring[:-1]+'+C$$<br>'
        return string
    if type(rule) == URule:
        u = SR(rule.u_var)
        uf = SR(rule.u_func)
        if rule.constant is None:
            const = SR(1)
        else:
            const = SR(rule.constant)
        subcontext = SR(rule.substep.context)*const
        string += 'Substitute ${}={}$, $d{}={}d{}$<br>'.format(latex(u),latex(uf),latex(u),latex(uf.diff(x)),latex(x))
        string += '$$\int {}d{}=\int {}d{}$$<br><br>'.format(latex(f),latex(x),latex(subcontext),latex(u))
        string += rule_string(rule.substep, depth-1)
        # string += 'So, $$\int {} d {} = {} + C$$<br>'.format(latex(f),latex(x),latex(SR(manualintegrate(f._sympy_(), x._sympy_()))))
        return string
    if type(rule) == ExpRule:
        r
        string += 'The integral of an exponential function is itself divided by the natural logarithm of the base, so<br>'
        string += '$$ \int {} d{} = {} + C$$<br>'.format(latex(f),latex(x),latex(exp(exponent)/log(base)))
        return string
    if type(rule) == PartsRule:
        u = SR(rule.u)
        dv = SR(rule.dv)
        du = SR(u.diff(x))
        v = SR(_manualintegrate(rule.v_step))
        string += 'Integration by parts, taking $u = {}$, $dv = {}d{}$<br>'.format(latex(u),latex(dv),latex(x))
        string += 'So $du={}$, and $\displaystyle v=\int {} d{}$<br>'.format(latex(du),latex(dv),latex(x))
        string += rule_string(rule.v_step, depth-1)
        string += 'So $$\int {} d{} =  {} - \int {} d{}$$<br>'.format(latex(f),latex(x),latex(u*v),latex(v*du),latex(x))
        string += rule_string(rule.second_step,depth-1)
        return string
    if type(rule) == TrigRule:
        result = SR(_manualintegrate(rule))
        string += '$$\int {} d{} = {} + C$$'.format(latex(f), latex(x), latex(result))
        return string
    if type(rule) == CyclicPartsRule:
        partialstring = '\int {} d{} = '.format(latex(f),latex(x))
        sign = -1
        resul = 0
        count = 0
        for srule in rule.parts_rules:
            u = SR(srule.u)
            dv = SR(srule.dv)
            du = SR(u.diff(x))
            v = SR(_manualintegrate(srule.v_step))
            string += 'Integrate by parts, taking $u = {}$, $dv= {}d{}$ <br>'.format(latex(u),latex(dv),latex(x))
            string += 'So $du={}$, and $\displaystyle v=\int {} d{}$ <br>'.format(latex(du),latex(dv),latex(x))
            string += rule_string(srule.v_step,depth-1)
            if sign == -1:
                resul += u*v
                if count :
                    partialstring += '+{}'.format(latex(u*v))
                else:
                    partialstring += '{}'.format(latex(u*v))
                string += 'So $${} - \int {} d{}$$'.format(partialstring, latex(v*du), latex(x))
            else:
                resul -= u*v
                partialstring += '-{}'.format(latex(u*v))
                string += 'So $${} + \int {} d{}$$<br>'.format(partialstring, latex(v*du), latex(x))
            sign *= -1
            count += 1
        resul = resul / (1-SR(rule.coefficient))
        string += 'Notice the integrand has repeated, so moving to one side we get:'
        string += '$$\int {} d{} = {} +C$$'.format(latex(f), latex(x), latex(resul))
        return string
    else:
        return ''
