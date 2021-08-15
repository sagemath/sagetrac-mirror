r"""
Symbolic Boolean Operators
"""

from sage.misc.latex import latex
from sage.symbolic.function import BuiltinFunction
from sage.symbolic.expression import is_Expression


def _trivial_bool(x):
    if is_Expression(x):
        if x.is_trivial_zero():
            return False
        if (x - 1).is_trivial_zero():
            return True
        return x
    else:
        return bool(x)


class AndSymbolic(BuiltinFunction):

    def __init__(self):
        r"""

        EXAMPLES::

            sage: and_symbolic(True, True)
            1
            sage: and_symbolic(False, True)
            0
            sage: and_symbolic(x>0, True)
            x > 0
            sage: and_symbolic(x>0, True, x<1)
            and_symbolic(x > 0, x < 1)
            sage: and_symbolic(x>0, False, x<1)
            0

            sage: var('y')
            y
            sage: and_symbolic(and_symbolic(x>0, x<1), and_symbolic(y>0, y<1))
            and_symbolic(x > 0, x < 1, y > 0, y < 1)

        TESTS:

        Conversions to other systems::

            sage: and_symbolic(x>0, x>1)._sympy_()
            (x > 0) & (x > 1)
        """
        BuiltinFunction.__init__(self, 'and_symbolic', nargs=0,
                                 conversions=dict(sympy='And'))

    def _eval_(self, *args):

        def clauses(args):
            for arg in args:
                if _trivial_bool(arg) is not True:
                    if is_Expression(arg) and isinstance(arg.operator(), AndSymbolic):
                        yield from clauses(arg.operands())
                    else:
                        yield arg

        flattened_args = list(clauses(args))

        if any(_trivial_bool(arg) is False for arg in flattened_args):
            return False

        if not flattened_args:
            # trivially true
            return True

        if len(args) == len(flattened_args) and all(x is y for x, y in zip(args, flattened_args)):
            # leave unevaluated
            return

        if len(flattened_args) == 1:
            return flattened_args[0]

        return and_symbolic(*flattened_args)

    def _print_latex_(self, *args):
        r"""
        EXAMPLES::

            sage: latex(and_symbolic(x>0, x<1))
            x > 0 \wedge x < 1

        """
        return r" \wedge ".join(latex(arg) for arg in args)

    def _giac_init_evaled_(self, *args):
        """
        TESTS::

            sage: and_symbolic(x>0, x<1)._giac_()
            ((sageVARx>0) and (1>sageVARx))
        """
        # modeled after Function_gamma_inc_lower._mathematica_init_evaled_
        args_giac = []
        for a in args:
            if isinstance(a, str):
                args_giac.append(a)
            elif hasattr(a, '_giac_init_'):
                args_giac.append(a._giac_init_())
            else:
                args_giac.append(str(a))
        return " and ".join('(' + arg + ')'
                            for arg in args_giac)

    def _maxima_init_evaled_(self, *args):
        """
        TESTS::

            sage: and_symbolic(x<0, x>1)._maxima_()
            _SAGE_VAR_x < 0 and _SAGE_VAR_x > 1
        """
        # modeled after Function_gamma_inc_lower._mathematica_init_evaled_
        args_maxima = []
        for a in args:
            if isinstance(a, str):
                args_maxima.append(a)
            elif hasattr(a, '_maxima_init_'):
                args_maxima.append(a._maxima_init_())
            else:
                args_maxima.append(str(a))
        return " and ".join('(' + arg + ')'
                           for arg in args_maxima)


and_symbolic = AndSymbolic()


class OrSymbolic(BuiltinFunction):

    def __init__(self):
        r"""

        EXAMPLES::

            sage: or_symbolic(True, True)
            1
            sage: or_symbolic(False, True)
            1
            sage: or_symbolic(x>0, False)
            x > 0
            sage: or_symbolic(x<0, False, x>1)
            or_symbolic(x < 0, x > 1)
            sage: or_symbolic(x<0, True, x>1)
            1

            sage: var('y')
            y
            sage: or_symbolic(or_symbolic(x<0, x>1), or_symbolic(y<0, y>1))
            or_symbolic(x < 0, x > 1, y < 0, y > 1)

        TESTS:

        Conversions to other systems::

            sage: or_symbolic(x<0, x>1)._sympy_()
            (x > 1) | (x < 0)
        """
        BuiltinFunction.__init__(self, 'or_symbolic', nargs=0,
                                 conversions=dict(sympy='Or'))

    def _eval_(self, *args):

        def clauses(args):
            for arg in args:
                if _trivial_bool(arg) is not False:
                    if is_Expression(arg) and isinstance(arg.operator(), OrSymbolic):
                        yield from clauses(arg.operands())
                    else:
                        yield arg

        flattened_args = list(clauses(args))

        if any(_trivial_bool(arg) is True for arg in flattened_args):
            return True

        if not flattened_args:
            # trivially false
            return False

        if len(args) == len(flattened_args) and all(x is y for x, y in zip(args, flattened_args)):
            # leave unevaluated
            return

        if len(flattened_args) == 1:
            return flattened_args[0]

        return or_symbolic(*flattened_args)

    def _print_latex_(self, *args):
        r"""
        EXAMPLES::

            sage: latex(or_symbolic(x>0, x<1))
            x > 0 \vee x < 1

        """
        return r" \vee ".join(latex(arg) for arg in args)

    def _giac_init_evaled_(self, *args):
        """
        TESTS::

            sage: or_symbolic(x<0, x>1)._giac_()
            ((0>sageVARx) or (sageVARx>1))
        """
        # modeled after Function_gamma_inc_lower._mathematica_init_evaled_
        args_giac = []
        for a in args:
            if isinstance(a, str):
                args_giac.append(a)
            elif hasattr(a, '_giac_init_'):
                args_giac.append(a._giac_init_())
            else:
                args_giac.append(str(a))
        return " or ".join('(' + arg + ')'
                           for arg in args_giac)

    def _maxima_init_evaled_(self, *args):
        """
        TESTS::

            sage: or_symbolic(x<0, x>1)._maxima_()
            _SAGE_VAR_x < 0 or _SAGE_VAR_x > 1
        """
        # modeled after Function_gamma_inc_lower._mathematica_init_evaled_
        args_maxima = []
        for a in args:
            if isinstance(a, str):
                args_maxima.append(a)
            elif hasattr(a, '_maxima_init_'):
                args_maxima.append(a._maxima_init_())
            else:
                args_maxima.append(str(a))
        return " or ".join('(' + arg + ')'
                           for arg in args_maxima)


or_symbolic = OrSymbolic()


class NotSymbolic(BuiltinFunction):

    def __init__(self):
        r"""

        EXAMPLES::

            sage: not_symbolic(True)
            0
            sage: not_symbolic(False)
            1
            sage: not_symbolic(x>0)
            not_symbolic(x > 0)
            sage: not_symbolic(not_symbolic(x>0))
            x > 0

        TESTS:

        Conversions to other systems::

            sage: var('P')
            P
            sage: not_symbolic(P)._sympy_()
            ~P
            sage: not_symbolic(P)._giac_()
            not(sageVARP)
            sage: not_symbolic(P)._maxima_()
            not _SAGE_VAR_P
        """
        BuiltinFunction.__init__(self, 'not_symbolic', nargs=1,
                                 conversions=dict(sympy='Not',
                                                  maxima='not',
                                                  giac='not'))

    def _eval_(self, arg):

        bool_arg = _trivial_bool(arg)

        if bool_arg is True:
            return False
        if bool_arg is False:
            return True
        if is_Expression(arg) and isinstance(arg.operator(), NotSymbolic):
            return arg.op[0]
        # leave unevaluated
        return

    def _print_latex_(self, arg):
        r"""
        EXAMPLES::

            sage: from sage.functions.boolean import not_symbolic
            sage: latex(not_symbolic(x>0))
            \neg( x > 0 )

        """
        return r"\neg(" + latex(arg) + r")"


not_symbolic = NotSymbolic()
