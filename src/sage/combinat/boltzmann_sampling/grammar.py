"""
Context-free grammars for Boltzmann generation.

Grammars use the basic operators of the symbolic method of analytic
combinatorics to specify labelled or unlabelled combinatorial classes. For
instance, binary tree can be specified by ``B = leaf + Z * B * B`` which, using
the syntax implemented in this module, looks like:

EXAMPLE::

    sage: z = Atom("z")
    sage: leaf = Atom("leaf", size=0)
    sage: Grammar(rules={"B": Union(leaf, Product(z, "B", "B"))})
    B -> Union(leaf, Product(z, B, B))

Grammars are not limited to a single rule:

EXAMPLE::

    sage: z = Atom("z")
    sage: leaf = Atom("leaf", size=0)
    sage: Grammar(rules={
    ....:     "T": Product(z, "S"),
    ....:     "S": Union(leaf, Product("T", "S")),
    ....: })
    S -> Union(leaf, Product(T, S))
    T -> Product(z, S)


AUTHORS:

- Matthieu Dien (2019): initial version

- Martin Pépin (2019): initial version
"""

from sage.all import latex, var, SR
from sage.structure.sage_object import SageObject


class Rule(SageObject):
    """The super class of all grammar rules.

    Should not be instantiated directly.
    """


class Atom(Rule):
    """Terminal symbol of a grammar.

    EXAMPLES::

        sage: z = Atom("z")

        sage: z = Atom("z", size=4)

        sage: x = Atom("x", size=0)
    """

    def __init__(self, name, size=1):
        """Create an Atom.

        INPUT:

        - ``name`` -- string; names should be unique across a grammar

        - ``size`` -- integer (default: 1)
        """
        self.name = name
        self.size = size

    def _latex_(self):
        r"""Return the LaTeX representation of an atom.

        EXAMPLES::

            sage: z = Atom("z")
            sage: latex(z)
            z

            sage: z = Atom("longname")
            sage: latex(z)
            \textrm{longname}

            sage: z = Atom("z", size=4)
            sage: latex(z)
            z^4

            sage: x = Atom("x", size=0)
            sage: latex(x)
            1
        """
        nice_name = self.name
        if len(self.name) != 1:
            nice_name = r"\textrm{{{}}}".format(self.name)

        if self.size == 1:
            return nice_name
        elif self.size == 0:
            return "1"
        else:
            return "{}^{}".format(nice_name, self.size)

    def _repr_(self):
        return self.name

    def combsys(self):
        if self.size > 0:
            return var(self.name) ** self.size
        else:
            return SR(1)


class Ref(Rule):
    """Non terminal symbols of a grammar.

    Instances of this class represent recursive references to non-terminal
    symbols inside grammar rules. In general you should not use this class
    directly.
    """

    def __init__(self, name):
        """Create a reference to a non-terminal.

        INPUT:

        - ``name`` -- string; the name of the referenced terminal symbol
        """
        self.name = name

    def _latex_(self):
        r"""Return the LaTeX representation of a non-terminal symbol.

        EXAMPLES::

            sage: latex(Ref("X"))
            X

            sage: latex(Ref("LongName"))
            \textrm{LongName}
        """
        nice_name = self.name
        if len(self.name) != 1:
            nice_name = r"\textrm{{{}}}".format(self.name)
        return nice_name

    def _repr_(self):
        # shall it be copy-pastable?
        return self.name

    def combsys(self):
        return var(self.name)


def _to_rule(r):
    if isinstance(r, Rule):
        return r
    else:
        return Ref(r)


class Union(Rule):
    """Union of two or more rules.

    D = Union(A, B, C) corresponds to the following grammar in BNF syntax:
    ``D ::= A | B | C``

    EXAMPLES::

        sage: Union("A", "B", "C")
        Union(A, B, C)

        sage: z = Atom("z")
        sage: Union(z, "A")
        Union(z, A)
    """

    def __init__(self, *args):
        """Build a union of two or more rules.

        INPUT:

        - ``args`` -- list of strings or Rules; strings are interpreted as Refs
        """
        if len(args) < 2:
            if len(args) == 0:
                msg = "Empty unions are forbidden"
            if len(args) == 1:
                arg, = args
                msg = (
                    "Unions of a single element are forbidden, "
                    "use {} directly instead".format(args)
                )
            raise ValueError(msg)
        self.args = [_to_rule(arg) for arg in args]

    def _latex_(self):
        """Return the LaTeX representation of a union.

        EXAMPLES::

            sage: A, B = Ref("A"), Ref("B")
            sage: latex(Union(A, B))
            A + B

            sage: z = Atom("z", size=2)
            sage: latex(Union(z, "A"))
            z^2 + A
        """
        return " + ".join(map(latex, self.args))

    def _repr_(self):
        # shall it be copy-pastable?
        return "Union({})".format(", ".join(map(repr, self.args)))

    def combsys(self):
        res = SR(0)
        for rule in self.args:
            res += rule.combsys()
        return res


class Product(Rule):
    """Product of two or more rules.

    Product(A, B, C) corresponds to the following grammar in BNF syntax:
    x ::= A × B × C

    EXAMPLES::

        sage: Product("A", "B", "C")
        Product(A, B, C)

        sage: z = Atom("z")
        sage: Product(z, "A")
        Product(z, A)
    """

    def __init__(self, *args):
        """Build the product of two or more rules.

        INPUT:

        - ``args`` -- list of strings or Rules; strings are interpreted as Refs
        """
        if len(args) < 2:
            if len(args) == 0:
                msg = "Empty products are forbidden"
            if len(args) == 1:
                arg, = args
                msg = (
                    "Products of a single element are forbidden, "
                    "use {} directly instead".format(args)
                )
            raise ValueError(msg)
        self.args = [_to_rule(arg) for arg in args]

    def _latex_(self):
        r"""Return the LaTeX representation of a product.

        EXAMPLES::

            sage: A, B = Ref("A"), Ref("B")
            sage: latex(Product(A, B))
            A \times B

            sage: z = Atom("z", size=2)
            sage: latex(Product(z, "A"))
            z^2 \times A
        """
        def wrap_latex(rule):
            if isinstance(rule, Union):
                return "({})".format(latex(rule))
            else:
                return latex(rule)
        return r" \times ".join(map(wrap_latex, self.args))

    def _repr_(self):
        # shall it be copy-pastable?
        return "Product({})".format(", ".join(map(repr, self.args)))

    def combsys(self):
        res = SR(1)
        for rule in self.args:
            res *= rule.combsys()
        return res


class Seq(Rule):
    """Sequence (with order) of rule.
    """

    def __init__(self, arg):
        """

        INPUT:

        - ``arg`` -- a rule or the name of a grammar symbol (string)
        """
        #TODO : add
        self.arg = _to_rule(arg)

    def _latex_(self):
        r"""Return the LaTeX representation of sequence

        EXAMPLES::
            sage: A, B = Ref("A"), Ref("B")
            sage: latex(Seq(A))
            {\sc Seq}\left(A\right)

            sage: latex(Seq(Union("A", "B")))
            {\sc Seq}\left(A + B\right)
        """
        return r"{{\sc Seq}}\left({}\right)".format(latex(self.arg))
    
    def _repr_(self):
        # shall it be copy-pastable?
        return "Seq({})".format(self.arg)

    def combsys(self):
        return SR(1/(1-self.arg.combsys()))
    

class Grammar(SageObject):
    """Context free grammars."""

    def __init__(self, rules=None, labelled=False):
        r"""Create a grammar.

        INPUT:

        - ``rules`` (default: None) -- dictionary mapping strings (non-terminal
          names) to Rules

        - ``labelled`` (default: False) -- whether the atoms of the grammar are
          labelled or not. In the labelled case, the generator module will draw
          objects according to the labelled distribution (i.e.
          `P_x[a] = \frac{x^{|a|}{|a|!A(x)}`).

        EXAMPLES::

            sage: z = Atom("z")
            sage: eps = Atom("eps", size=0)

            sage: Grammar(rules={"S": Union(eps, Product(z, "S"))})
            S -> Union(eps, Product(z, S))

            sage: Grammar(rules={"B": Union(eps, Product(z, "B", "B"))})
            B -> Union(eps, Product(z, B, B))

            sage: g = Grammar()
            sage: g.set_rule("D", Union(z, Product(z, "S", z)))
            sage: g.set_rule("S", Union(eps, Product("D", "S")))
            sage: g
            D -> Union(z, Product(z, S, z))
            S -> Union(eps, Product(D, S))
        """
        self.labelled = labelled
        self.rules = {}
        rules = rules or {}
        for name, rule in rules.items():
            self.set_rule(name, rule)

    def set_rule(self, name, rule):
        """Add a rule to the grammar.

        INPUT:

        - ``name`` -- string; the name of the non-terminal being added. If
          a non-terminal with the same name was already present in the grammar,
          it is replaced.

        - ``rule`` -- a Rule

        EXAMPLES:

            sage: g = Grammar()
            sage: g.set_rule("A", Union("B", "C"))
            sage: g.rules["A"]
            Union(B, C)
        """
        rule = _to_rule(rule)
        self.rules[name] = rule

    def _latex_(self):
        r"""Return a LaTeX representation of the grammar.

        EXAMPLES::

            sage: y, z = Atom("y"), Atom("z")
            sage: leaf = Atom("leaf", size=0)
            sage: g = Grammar(rules={
            ....:     "T": Union(leaf, Product(y, "T"), Product(z, "T", "T"))
            ....:  })
            sage: latex(g)
            T = 1 + y \times T + z \times T \times T

            sage: g = Grammar(rules={
            ....:     "A": Ref("B"),
            ....:     "B": Product("C", "D")
            ....: })
            sage: latex(g)
            \begin{cases}
            A &= B \\
            B &= C \times D
            \end{cases}

        """
        if len(self.rules) <= 1:
            (name, rule), = self.rules.items()
            return "{} = {}".format(name, latex(rule))
        else:
            inside = " \\\\ \n".join(
                "{} &= {}".format(name, latex(rule))
                for name, rule in sorted(self.rules.items())
            )
            return "\\begin{{cases}}\n{}\n\\end{{cases}}".format(inside)

    def _repr_(self):
        return "\n".join(
            "{} -> {}".format(non_terminal, expr)
            for non_terminal, expr in sorted(self.rules.items())
        )

    def combsys(self):
        return {name: rule.combsys() for name, rule in self.rules.items()}
