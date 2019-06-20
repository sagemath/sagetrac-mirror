"""
Various oracle implementations for Boltzmann sampling.

Oracles are used to get (often approximate) values of generating functions.
Thanks to the symbolic method, functionnal equations can be derived from
grammar specifications. This module implements some mechanics to approximate
genrating functions based on these equations.

Currently two oracles are implemented:

- :class:`SimpleOracle` implements approximation by simple iteration of the
  equations.

- :class:`OracleFromFunctions` wraps an generating function given in the form
  of a python ore sage function as an oracle.

AUTHORS:
- Matthieu Dien (2019): initial version
- Martin Pépin (2019): initial version
"""

from sage.all import SR, latex


def _maximum_norm(y):
    return max(map(abs, y))


def _dict_diff(y, yp):
    return {k: y[k] - yp[k] for k in y.keys()}


class SimpleOracle:
    """Simple oracle for critical Boltzmann sampling based on iteration.

    EXAMPLES::

        sage: leaf = Atom("leaf", size=0)
        sage: z = Atom("z")
        sage: g = Grammar(rules={"B": Union(leaf, Product(z, "B", "B"))})
        sage: oracle = SimpleOracle(g)
        sage: oracle("z")  # abs tol 0.0001
        0.25

        sage: oracle("B")  # abs tol 0.0001
        2
    """

    def __init__(
        self, grammar, zstart=0, zmin=0.0, zmax=1.0, e1=0.0001, e2=0.0001
    ):
        """Create an oracle and annotate a grammar with the computed weights.

        INPUT:

        - ``grammar`` -- a Grammar

        - ``zstart`` -- number (default: 0); where to start the iteration

        - ``zmax`` -- number (default: 0); an upper bound on the singularity

        - ``zmin`` -- number (default: 0); a lower bound on the singularity

        - ``e1`` -- number (default: 0.001); TODO: explain

        - ``e2`` -- number (default: 0.002); TODO: explain
        """
        self.zstart = zstart
        self.zmin = zmin
        self.zmax = zmax
        self.e1 = e1
        self.e2 = e2

        self.grammar = grammar
        self.z = SR.symbol()
        self.combsys = self._normalize_combsys(grammar._to_combsys())
        self.weights = None
        self._compute_weights()
        self._register_in_grammar()

    def _normalize_combsys(self, combsys):
        cs = dict(combsys)
        variables = [v for x in cs.values() for v in x.variables()]
        for v in variables:
            if str(v)[:5] == "_var_":
                cs[v] = self.z
            if str(v)[:5] == "_eps_":
                cs[v] = SR(1)

        cs[self.z] = self.z
        return cs

    def _eval_combsys(self, z):
        y = {k: 0 for k in self.combsys.keys()}
        y[self.z] = z
        yp = {k: v.subs(y) for k, v in self.combsys.items()}
        yp = {k: v.subs(yp) for k, v in self.combsys.items()}

        while _maximum_norm(_dict_diff(y, yp)) > self.e2:
            y = yp
            yp = {k: v.subs(y) for k, v in self.combsys.items()}
        return yp

    def _diverge(self, y):
        return any(x < 0 or x > 1 / self.e1 for x in y)

    def _find_singularity(self):
        zstart = self.zstart
        zmax = self.zmax
        zmin = self.zmin

        while zmax - zmin > self.e1:
            y = self._eval_combsys(zstart)
            if self._diverge(y):
                zmax = zstart
                zstart = (zmin + zstart) / 2
            else:
                zmin = zstart
                zstart = (zmax + zstart) / 2

        return self._eval_combsys(zstart)

    def _compute_weights(self):
        self.weights = self._find_singularity()

    def _repr_(self):
        return "SimpleOracle for {}".format(latex(self.grammar))

    def __call__(self, rule):
        """Evaluate a rule at the computed main singularity."""
        if SR(rule) in self.weights:
            return self.weights[SR(rule)]
        elif SR("_var_" + rule) in self.weights:
            return self.weights[SR("_var_" + rule)]
        elif SR("_eps_" + rule) in self.weights:
            return self.weights[SR("_eps_" + rule)]
        else:
            raise KeyError(rule)

    def _register_in_grammar(self):
        self.grammar.annotate(self)


class OracleFromFunctions:
    """Wrapper for generating functions when they are known.

    In the case where the generating functions of all symbols in the grammar
    are knwon, this class wraps them as an oracle.
    """

    def __init__(self, variables, gen_funs):
        """Wrap generating functions as an oracle.

        INPUT:

        - ``variables`` -- dictionary mapping strings (atom names) to numbers
          (their values)

        - ``gen_funs`` -- dictionary mapping strings (non-terminal names) to
          functions (their generating series). These functions should accept
          the variables from the first argument as named argument.

        EXAMPLES::

            sage: B(z) = (1 - sqrt(1 - 4 * z)) / (2 * z)
            sage: oracle = OracleFromFunctions({"z": 1/4}, {"B": B})
            sage: oracle("z")  # abs tol 0.0000001
            0.25

            sage: oracle("B")  # abs tol 0.0000001
            2
        """
        self.variables = variables
        self.gfs = gen_funs

    def __call__(self, name):
        """Evaluate a rule at the computed main singularity."""
        if name in self.variables:
            return self.variables[name].n()
        elif name in self.gfs:
            return self.gfs[name](**self.variables).n()
        else:
            raise ValueError("Unknown name {}".format(name))
