r"""
Symbolic Boolean Operators
"""

from sage.misc.latex import latex
from sage.symbolic.function import BuiltinFunction
from sage.symbolic.expression import Expression
from sage.symbolic.ring import SR

class AndSymbolic(BuiltinFunction):

    def __init__(self):
        r"""

        EXAMPLES::

            sage: from sage.functions.boolean import and_symbolic
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

            sage: and_symbolic(x>0, x>1)._sympy_()
            (x > 0) & (x > 1)

        """
        BuiltinFunction.__init__(self, 'and_symbolic', nargs=0,
                                 conversions=dict(sympy='And'))

    def _eval_(self, *args):

        if any(arg in (0, False) for arg in args):
            return False

        if not args:
            # trivially true
            return True

        if not any(arg in (1, True) for arg in args):
            # leave unevaluated
            return

        args = [arg for arg in args if arg not in (1, True)]

        if not args:
            # trivially true
            return True

        if len(args) == 1:
            return args[0]

        return and_symbolic(*args)

    def _print_latex_(self, *args):
        r"""
        EXAMPLES::

            sage: from sage.functions.boolean import and_symbolic
            sage: latex(and_symbolic(x>0, x<1))
            x > 0 \wedge x < 1

        """
        return r" \wedge ".join(latex(arg) for arg in args)

and_symbolic = AndSymbolic()
