# coding: utf-8
"""
Various oracle implementations for Boltzmann sampling.

Oracles are used to get (often approximate) values of generating functions.
Thanks to the symbolic method, functionnal equations can be derived from
grammar specifications. This module implements some mechanics to approximate
generating functions based on these equations.

Currently two oracles are implemented:

- :class:`SimpleOracle` implements approximation by simple iteration of the
  equations.

- :class:`OracleFromFunctions` wraps an generating function given in the form
  of a python ore sage function as an oracle.

The entry point of these algorithms is the function oracle which determines
the oracle to use given its inputs.

AUTHORS:
- Matthieu Dien (2019): initial version
- Martin Pépin (2019): initial version
"""

from sage.structure.sage_object import SageObject
from sage.rings.infinity import Infinity as oo
from sage.all import SR, var, RR, vector
from .grammar import Grammar

def _diverge(dict, threshold=1e2):
    return any(x not in RR or x < 0 or x > threshold for x in dict.values())

def oracle(sys, **kargs):
    """Build different oracles given different inputs

    EXAMPLES::

        sage: from sage.combinat.boltzmann_sampling.oracle import SimpleOracle, OracleFromFunctions

        sage: leaf = Atom("leaf", size=0)
        sage: z = Atom("z")
        sage: g = Grammar(rules={"B": Union(leaf, Product(z, "B", "B"))})
        sage: oracle(g)
        SimpleOracle({'B': B^2*z + 1, 'z': z})

        sage: z = var('z')
        sage: oracle({'B': (1-sqrt(1-4*z))/(2*z), 'z': z})
        OracleFromFunctions({'B': -1/2*(sqrt(-4*z + 1) - 1)/z, 'z': z})
    """

    if isinstance(sys, Grammar):
        return SimpleOracle(sys, **kargs)
    elif isinstance(sys, dict):
        return OracleFromFunctions(sys, **kargs)


class SimpleOracle(SageObject):
    """Simple oracle for critical Boltzmann sampling based on iteration.

    EXAMPLES::

        sage: from sage.combinat.boltzmann_sampling.oracle import SimpleOracle

        sage: leaf = Atom("leaf", size=0)
        sage: z = Atom("z")
        sage: g = Grammar(rules={"B": Union(leaf, Product(z, "B", "B"))})
        sage: oracle = SimpleOracle(g)
        sage: oracle.eval_rule("z", {"z":1/4}) # abs tol 0.01
        0.25

        sage: oracle.eval_rule("B", {"z":1/4}) # abs tol 0.01
        2

        sage: oracle.eval_rule("B", {"z":1/17}) # abs tol 0.01
        1.06696562634075
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
        self.combsys.update({v: var(v) for v in self.terminals})

    def eval_combsys(self, z):
        """Compute a numerical evaluation of the combinatorial system
        at a given point ``z`` with the oracle's precision

        INPUT:

        - ``z`` -- a dictionary mapping strings (names of the variables) to
          numerical values

        OUTPUT: a dictionary associating symbols of the grammar
        to value of their generating functions at the input point.

        EXAMPLES::

            sage: leaf = Atom("leaf", size=0)
            sage: z = Atom("z")
            sage: g = Grammar(rules={"B": Union(leaf, Product(z, "B", "B"))})
            sage: o = oracle(g)
            sage: o.eval_combsys({"z": 1/17}) # abs tol 1e-3
            {'B': 1.06696559462842, 'z': 0.0588235294117647}
        """

        values = {k: RR(0) for k in self.non_terminals}
        values.update(z)
        new_values = {k: RR(self.combsys[k].subs(**values)) for k in values.keys()}
        
        while vector((values[k] - new_values[k] for k in values.keys())).norm(oo) > self.precision \
        and not _diverge(new_values):
            values = new_values
            new_values = {k: RR(self.combsys[k].subs(**values)) for k in values.keys()}
        return new_values

    def eval_rule(self, name, z):
        """Compute a numerical evaluation of the grammar rule ``name`` at a
        given point ``z``.

        INPUT:

        - ``name`` -- a string corresponding to a grammar non-terminal symbol

        - ``z`` -- a dictionary from string (name of the variables) to
          numerical value
        """

        values = self.eval_combsys(z)
        return values[name]

    def _repr_(self):
        return "SimpleOracle({})".format(self.combsys)


    def find_singularity(self, precision=1e-6, zmin=1e-9, zmax=1.-1e-9, divergence=1e2):
        """Given an oracle for a combinatorial system try to find the singularity.
        
        The algorithm proceed by dichotomic search. The divergence parameter allows
        to decide of the divergence of system.
    
        EXAMPLE::

            sage: from sage.combinat.boltzmann_sampling.oracle import SimpleOracle

            sage: leaf = Atom("leaf", size=0)
            sage: z = Atom("z")
            sage: g = Grammar(rules={"B": Union(leaf, Product(z, "B", "B"))})
            sage: o = SimpleOracle(g)
            sage: o.find_singularity()["z"] # abs tol 1e-6
            0.25
        """

        y = None
        zstart = zmin
        while zmax - zmin > precision:
            try:
                y = self.eval_combsys({v: zstart for v in self.terminals})
                if _diverge(y, threshold=divergence):
                    zmax = zstart
                    zstart = (zmin + zstart) / 2
                else:
                    zmin = zstart
                    zstart = (zmax + zstart) / 2
            except ValueError:
                zmax = zstart
                zstart = (zmin + zstart) / 2

        return self.eval_combsys({v: zmin for v in self.terminals})


class OracleFromFunctions(SageObject):
    """Wrapper for generating functions when they are known.

    In the case where the generating functions of all symbols in the grammar
    are knwon, this class wraps them as an oracle.
    """

    def __init__(self, sys, precision_ring=SR):
        """Wrap generating functions as an oracle.

        INPUT:

        - ``sys`` -- dictionary mapping strings (non-terminal names) to
          functions (their generating series).

        - ``precision_ring`` -- do the computation in the given ring (default:
          Symbolic Ring)

        EXAMPLES::

            sage: from sage.combinat.boltzmann_sampling.oracle import OracleFromFunctions

            sage: B(z) = (1 - sqrt(1 - 4 * z)) / (2 * z)
            sage: oracle = OracleFromFunctions({"z": z, "B": B})
            sage: oracle.eval_rule("z", {"z":1/4})
            1/4

            sage: oracle.eval_rule("B", {"z":1/4})
            2
        """
        self.sys = sys
        self.precision_ring = precision_ring
        self.terminals = {str(v) for f in sys.values() for v in f.variables()}
        self.singularity = None
        
    def eval_combsys(self, z):
        """Compute an evaluation of the combinatorial system
        at a given point ``z`` with the oracle's precision.

        INPUT:

        - ``z`` -- a dictionary mapping strings (name of the variables) to
          numerical values

        OUTPUT: a dictionary associating symbols of the grammar
        to value of their generating functions at the input point.
        """
        values = {k: self.eval_rule(k, z) for k in self.sys.keys()}
        values.update(z)
        return values

    def eval_rule(self, name, z):
        """Compute a evaluation of the grammar rule ``name``
        at a given point ``z`` with the oracle's precision.

        INPUT:

        - ``name`` -- a string corresponding to a grammar non-terminal symbol

        - ``z`` -- a dictionary from string (name of the variables) to
          numerical value
        """

        # return self.precision_ring(self.sys[name].subs(**z))
        return self.sys[name](**z)

    def _repr_(self):
        return "OracleFromFunctions({})".format(self.sys)

    def set_singularity(self, singularity):
        self.singularity = singularity
    
    def find_singularity(self, *args, **kargs):
        """Given an oracle for a combinatorial system try to find the singularity.
        
        The algorithm proceed by dichotomic search. The divergence parameter allows
        to decide of the divergence of system.
    
        EXAMPLE::

            sage: from sage.combinat.boltzmann_sampling.oracle import OracleFromFunctions

            sage: B(z) = (1 - sqrt(1 - 4 * z)) / (2 * z)
            sage: oracle = OracleFromFunctions({"z": z, "B": B})
            sage: oracle.set_singularity({"z":1/4})
            sage: oracle.find_singularity()
            {'B': 2, 'z': 1/4}
        """

        if self.singularity is None:
            raise NotImplementedError("You should provide the singularity using set_singularity.")
        else:
            return self.eval_combsys(self.singularity)
