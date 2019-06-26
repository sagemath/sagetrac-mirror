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

from functools import reduce

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

        sage: z = Atom("z", size=4, labelled=True)

        sage: x = Atom("x", size=0)
    """

    def __init__(self, name, size=1, labelled=False):
        r"""Create an Atom.

        INPUT:

        - ``name`` -- string; names should be unique across a grammar

        - ``size`` -- integer (default: 1)

        - ``labelled`` (default: False) -- whether the atom is labelled or not.
          In the labelled case, the generator module will draw objects
          according to the labelled distribution (i.e. `P_x[a] =
          \frac{x^{|a|}{|a|!A(x)}`).
        """
        self.name = name
        self.size = size
        self._labelled = labelled

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

    def labelled(self):
        return self._labelled

    def atoms(self):
        return {self}

    def __hash__(self):
        return hash((self.name, self.size, self._labelled))

    def __eq__(self, other):
        return (
            isinstance(other, Atom)
            and self.name == other.name
            and self.size == other.size
            and self.labelled == other.labelled
        )


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

    def labelled(self):
        return False

    def atoms(self):
        return set()


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

    def labelled(self):
        return any((arg.labelled() for arg in self.args))

    def atoms(self):
        return reduce(lambda x, y: x | y, (arg.atoms() for arg in self.args))


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

    def labelled(self):
        return any((arg.labelled() for arg in self.args))

    def atoms(self):
        return reduce(lambda x, y: x | y, (arg.atoms() for arg in self.args))


class Seq(Rule):
    """Sequence (with order) of rule.
    """

    def __init__(self, arg):
        """

        INPUT:

        - ``arg`` -- a rule or the name of a grammar symbol (string)
        """
        # TODO: add
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
        return "Seq({})".format(self.arg)

    def combsys(self):
        return SR(1 / (1 - self.arg.combsys()))

    def labelled(self):
        return self.arg.labelled()

    def atoms(self):
        return self.arg.atoms()


class Grammar(SageObject):
    """Context free grammars."""

    def __init__(self, rules=None):
        r"""Create a grammar.

        INPUT:

        - ``rules`` (default: None) -- dictionary mapping strings (non-terminal
          names) to Rules

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

    def labelled(self):
        return any((expr.labelled() for expr in self.rules.values()))

    def atoms(self):
        """Return all the atoms (terminals) appearing in the grammar.

        EXAMPLE:

            sage: z = Atom("z")
            sage: e = Atom("e", size=0)
            sage: g = Grammar(rules={"B": Union(e, Product(z, "B", "B"))})
            sage: g.atoms()
            {e, z}
        """
        return reduce(
            lambda x, y: x | y, (expr.atoms() for expr in self.rules.values())
        )
