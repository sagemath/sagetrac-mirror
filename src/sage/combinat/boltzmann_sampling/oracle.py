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

from sage.structure.sage_object import SageObject
from sage.rings.infinity import Infinity as oo
from sage.all import SR, latex, var, RR, vector

from .grammar import Grammar


def oracle(sys, **kargs):
    """TODO: document.

    EXAMPLES::

        sage: z = Atom("z")
        sage: eps = Atom("eps", size=0)
        sage: grammar = Grammar(rules={"B": Union(eps, Product(z, "B", "B"))})

        sage: # Build a sampler from a grammar
        sage: oracle(grammar)

        sage: # Build a sampler from a dictionary of function
        sage: oracle({"B": lambda z: (1 - sqrt(1 - 4 * z)) / (2 * z)})
    """
    if isinstance(sys, Grammar):
        return SimpleOracle(sys, **kargs)
    elif isinstance(sys, dict):
        return OracleFromFunctions(sys, **kargs)


class SimpleOracle(SageObject):
    """Simple oracle for critical Boltzmann sampling based on iteration.

    EXAMPLES::

        sage: leaf = Atom("leaf", size=0)
        sage: z = Atom("z")
        sage: g = Grammar(rules={"B": Union(leaf, Product(z, "B", "B"))})
        sage: oracle = SimpleOracle(g)
        sage: oracle.eval_rule("z", {"z":1/4}) # abs tol 0.01
        0.25

        sage: oracle.eval_rule("B", {"z":1/4}) # abs tol 0.01
        2
    """

    def __init__(self, grammar, precision=1e-6):
        """Create an oracle and annotate a grammar with the computed weights.

        INPUT:

        - ``grammar`` -- a Grammar

        - ``precision`` -- number (default: 1e-6); TODO: explain
        """
        self.precision = precision

        self.combsys = grammar.combsys()

        # non terminal names of the grammar i.e. combinatorial classes
        self.non_terminals = set(self.combsys.keys())
        # terminal names of the grammar i.e. atoms
        self.terminals = {str(var) for expr in self.combsys.values()
                          for var in expr.variables()
                          if str(var) not in self.non_terminals}
        # all atoms are represented by the same variable
        self.combsys.update({v : var(v) for v in self.terminals})

    def eval_combsys(self, z):
        values = {k: RR(0) for k in self.non_terminals}
        values.update(z)
        new_values = {k: RR(self.combsys[k].subs(**values)) for k in values.keys()}

        while vector((values[k] - new_values[k] for k in values.keys())).norm(oo) > self.precision:
            values = new_values
            new_values = {k: RR(self.combsys[k].subs(**values)) for k in values.keys()}
        return new_values

    def eval_rule(self, name, z):
        values = self.eval_combsys(z)
        return values[name]

    def _repr_(self):
        return "SimpleOracle for {}".format(latex(self.combsys))


def find_singularity(oracle, precision=1e-6, zstart=0., zmin=0., zmax=1., divergence=1e3):
    """Given an oracle for a combinatorial system try to find the singularity.

    The algorithm proceed by dichotomic search. The divergence parameter allows
    to decide of the divergence of system.

    EXAMPLE::

        sage: leaf = Atom("leaf", size=0)
        sage: z = Atom("z")
        sage: g = Grammar(rules={"B": Union(leaf, Product(z, "B", "B"))})
        sage: oracle = SimpleOracle(g)
        sage: find_singularity(oracle)["z"] # abs tol 1e-6
        0.25
    """

    y = None
    while zmax - zmin > precision:
        y = oracle.eval_combsys({v: zstart for v in oracle.terminals})
        if any((x < 0 or x > divergence for x in y.values())):
            zmax = zstart
            zstart = (zmin + zstart) / 2
        else:
            zmin = zstart
            zstart = (zmax + zstart) / 2

    return oracle.eval_combsys({v: zmin for v in oracle.terminals})


class OracleFromFunctions(SageObject):
    """Wrapper for generating functions when they are known.

    In the case where the generating functions of all symbols in the grammar
    are knwon, this class wraps them as an oracle.
    """

    def __init__(self, sys, precision=oo):
        """Wrap generating functions as an oracle.

        INPUT:

        - ``variables`` -- dictionary mapping strings (atom names) to numbers
          (their values)

        - ``gen_funs`` -- dictionary mapping strings (non-terminal names) to
          functions (their generating series). These functions should accept
          the variables from the first argument as named argument.

        EXAMPLES::

            sage: B(z) = (1 - sqrt(1 - 4 * z)) / (2 * z)
            sage: oracle = OracleFromFunctions({"z": z, "B": B})
            sage: oracle.eval_rule("z", {"z":1/4})
            1/4

            sage: oracle.eval_rule("B", {"z":1/4})
            2
        """
        self.sys = sys
        self.precision = precision

        # Scalar field
        self.SF = None
        if precision == oo:
            self.SF = SR
        else:
            self.SF = RR

    def eval_combsys(self, z):
        return {k: self.eval_rule(k, z) for k in self.sys.keys()}

    def eval_rule(self, name, z):
        return self.SF(self.sys[name].subs(**z))
