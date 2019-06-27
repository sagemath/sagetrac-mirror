# coding: utf-8
r"""Boltzmann generator for Context-free grammars.

This module provides functions for generating combinatorial objects (i.e.
objects described by a combinatorial specification, see
:mod:`sage.combinat.boltzmann_sampling.grammar`) according to the Boltzmann
distribution.

Given an unlabelled combinatorial class A, the Boltzmann distribution of
parameter `x` is such that an object of size n is drawn with the probability
`\frac{x^n}{A(x)}` where `A(x)` denotes the ordinary generating function of A.
For labelled classes, this probability is set to `\frac{x^n}{n!A(x)}` where
`A(x)` denotes the exponential generating function of A. See [DuFlLoSc04] for
details.

By default, the objects produced by the generator are nested tuples of strings
(the atoms). For instance ``('z', ('z', 'e', 'e'), ('z', 'e', 'e'))`` is a
balanced binary tree with 3 internal nodes (z) and 4 leaves (e). To alter this
behaviour and generate other types of objects, you can specify a builder
function for each type of symbol in the grammar. The behaviour of the builders
is that they are applied in bottom up order to the structure "on the fly" during
generation. For instance, in order to generate Dyck words using the grammar for
binary trees, one can use a builder that return ``""`` for each leaf ``"(" +
left child + ")" + right child`` for each node. The builders will receive a
tuple for each product, a string for each atom and builders for unions should be
computed using the :func:`UnionBuilder` helper. See the example below for the
case of Dyck words.

EXAMPLES::

    sage: leaf = Atom("e", size=0)
    sage: z = Atom("z")
    sage: grammar = Grammar(rules={"B": Union(leaf, Product(z, "B", "B"))})
    sage: generator = Generator(grammar)

    sage: def leaf_builder(_):
    ....:     return ""

    sage: def node_builder(tuple):
    ....:     _, left, right = tuple
    ....:     return "(" + left + ")" + right

    sage: generator.set_builder("B", UnionBuilder(leaf_builder, node_builder))
    sage: dyck_word, _ = generator.gen("B", (10, 20))
    sage: dyck_word  # random
    "(()((()())))((())(()((()))))"
    sage: len(dyck_word) in range(20, 41)  # 2n in [20; 40]
    True

    sage: def leaf_builder(_):
    ....:     return BinaryTree()

    sage: def node_builder(tuple):
    ....:     _, left, right = tuple
    ....:     return BinaryTree([left, right])

    sage: generator.set_builder("B", UnionBuilder(leaf_builder, node_builder))
    sage: tree, _ = generator.gen("B", (10, 20))
    sage: tree.node_number() in range(10, 21)
    True

Note that the builders' mechanism can also be used to compute parameters on the
structure on the fly without building the whole structure such as the height of
the tree.

EXAMPLES::

    sage: leaf = Atom("e", size=0)
    sage: z = Atom("z")
    sage: grammar = Grammar(rules={"B": Union(leaf, Product(z, "B", "B"))})
    sage: generator = Generator(grammar)

    sage: def leaf_height(_):
    ....:     return 0

    sage: def node_height(tuple):
    ....:     _, left, right = tuple
    ....:     return 1 + max(left, right)

    sage: generator.set_builder("B", UnionBuilder(leaf_height, node_height))
    sage: height, _ = generator.gen("B", (10, 20))
    sage: height  # random
    6

A last example with labelled grammars:

EXAMPLE::

    sage: leaf = Atom("leaf")
    sage: u, z = Atom("u", labelled=True), Atom("z", labelled=True)
    sage: g = Grammar(rules={
    ....:     "T": Union(leaf, Product(u, "T"), Product(z, "T", "T"))
    ....: })
    sage: T(leaf, u, z) = ((1-u) - sqrt((1-u)^2 - 4 * z * leaf)) / (2 * z)
    sage: o = oracle({"T": T})
    sage: o.set_singularity({"leaf": 1/4, "u": 1/2, "z": 1/4})
    sage: generator = Generator(g, oracle=o)
    sage: tree, stats = generator.gen("T", (10, 20))
    sage: tree
    (('z', [1]),
     'leaf',
     (('u', [2]),
      (('u', [7]),
       (('u', [3]),
        (('z', [2]),
         (('u', [8]),
          (('u', [6]),
           (('u', [1]),
            (('z', [0]),
             (('u', [5]),
              (('u', [0]),
               'leaf')),
             'leaf')))),
        (('u', [4]),
         'leaf'))))))

REFERENCES:

.. [DuFlLoSc04] Philippe Duchon, Philippe Flajolet, Guy Louchard, and Gilles
   Schaeffer. 2004. Boltzmann Samplers for the Random Generation of
   Combinatorial Structures. Comb. Probab. Comput. 13, 4-5 (July 2004), 577-625.
   DOI=http://dx.doi.org/10.1017/S0963548304006315
"""

from functools import reduce

from sage.libs.gmp.random cimport gmp_randinit_set, gmp_randinit_default
from sage.libs.gmp.types cimport gmp_randstate_t
from sage.misc.randstate cimport randstate, current_randstate, SAGE_RAND_MAX

from .grammar import Atom, Product, Ref, Union, Seq
from .oracle import SimpleOracle

from libc.math cimport log, ceil

ctypedef enum options:
    # Grammar constructions
    REF,
    ATOM,
    UNION,
    PRODUCT,
    SEQ,
    # Special instructions used by the random generators
    TUPLE,
    LIST,
    FUNCTION,
    WRAP_CHOICE,
    WRAPPED_CHOICE

# ------------------------------------------------------- #
# Grammar preprocessing
# ------------------------------------------------------- #

# For performance reasons, we use integers rather that strings to identify
# symbols during generation and we drop the python classes in favour of tuples
# nested tuples. The following functions implement this transformation.
#
# All of this should remain hidden from the end user.

cdef _map_all_names_to_ids_expr(name_to_id, weights, expr):
    """Recursively transform an expression into a tuple of the form
    (RULE_TYPE, weight, args) where:
    - RULE_TYPE encodes the type of the expression [ATOM|REF|PRODUCT|...]
    - weight is the value of the generating function of this expression
    - args is auxilliary information (the name of an atom, the components of an
      union, ...)
    """
    if isinstance(expr, Ref):
        return (REF, weights[expr.name], name_to_id[expr.name])
    elif isinstance(expr, Atom):
        weight = 1. if expr.size == 0 else weights[expr.name] ** expr.size
        return (ATOM, weight, (expr.name, expr.size))
    elif isinstance(expr, Union):
        args = tuple((_map_all_names_to_ids_expr(name_to_id, weights, arg) for arg in expr.args))
        total_weight = sum([w for (__, w, __) in args])
        return (UNION, total_weight, args)
    elif isinstance(expr, Product):
        args = tuple((_map_all_names_to_ids_expr(name_to_id, weights, arg) for arg in expr.args))
        total_weight = reduce(
            lambda x, y: x * y,
            [w for (__, w, __) in args],
            1.
        )
        return (PRODUCT, total_weight, args)
    elif isinstance(expr, Seq):
        rule = _map_all_names_to_ids_expr(name_to_id, weights, expr.arg)
        return (SEQ, 1./(1. - rule[1]), rule)


cdef _map_all_names_to_ids_system(name_to_id, id_to_name, weights, rules):
    return [
        _map_all_names_to_ids_expr(name_to_id, weights, rules[id_to_name[i]])
        for i in range(len(name_to_id))
    ]

cdef _map_all_names_to_ids(rules):
    """Assign an integer (identifier) to each symbol in the grammar and compute
    two dictionaries:
    - one that maps the names to their identifiers
    - one that maps the identifier to the original names
    """
    name_to_id = {}
    id_to_name = {}
    for i, name in enumerate(rules.keys()):
        name_to_id[name] = i
        id_to_name[i] = name
    return name_to_id, id_to_name


# ------------------------------------------------------- #
# Generic builders
# ------------------------------------------------------- #

# When it chooses one possible derivation for a union rule, the random
# generator wraps the generated object in a tuple of the form ``(choice_id,
# object)`` so that the builder that will be called on this value has the
# information of which derivation was chosen. This helper function hides this
# machinery to the end user, allowing her to compose builders in a "high-level"
# manner.

def UnionBuilder(*builders):
    """Factory for generating builders for union rules.

    The recommended way to write a builder for a union rule is to use this
    helper: write a auxilliary builder for each component of the union and
    compose them with ``UnionBuilder``.

    Assume the symbol ``D`` is defined by ``Union("A", "B", "C")``, then
    defining a builder for D would look like:

    EXAMPLE::

        sage: def build_A(args):
        ....:     # Do something
        ....:     pass

        sage: def build_B(args):
        ....:     # Do something
        ....:     pass

        sage: def build_C(args):
        ....:     # Do something
        ....:     pass

        sage: build_D = UnionBuilder(build_A, build_B, build_C)

    For instance, for binary trees defined by ``B = Union(leaf, Product(z, B,
    B))``, this could be:

    EXAMPLE::

        sage: def build_leaf(_):
        ....:     return BinaryTree()

        sage: def build_node(args):
        ....:     z, left, right = args
        ....:     return BinaryTree([left, right])

        sage: build_binarytree = UnionBuilder(build_leaf, build_node)
    """
    def build(obj):
        index, content = obj
        builder = builders[index]
        return builder(content)
    return build

# The following functions generate the default builders for labelled and
# unlabelled structures. These builders produces nested tuples. Atom sizes and
# choices ids are omitted for readability.

cdef inline identity(x):
    return x

cdef inline first(x):
    a, __ = x
    return a

cdef inline ProductBuilder(builders):
    def build(terms):
        return tuple(builders[i](terms[i]) for i in range(len(terms)))
    return build

cdef inline SeqBuilder(builder):
    def build(terms):
        return list(builder(terms[i]) for i in range(len(terms)))
    return build

cdef make_default_builder(rule):
    """Generate the default builders for a rule.

    For use with Boltzmann samplers :mod:`sage.combinat.boltzmann_sampling`"""
    if isinstance(rule, Ref):
        return identity
    elif isinstance(rule, Atom):
        if rule.labelled():
            return identity
        else:
            return first
    elif isinstance(rule, Product):
        subbuilders = [make_default_builder(component) for component in rule.args]
        return ProductBuilder(subbuilders)
    elif isinstance(rule, Seq):
        subbuilder = make_default_builder(rule.arg)
        return SeqBuilder(subbuilder)
    elif isinstance(rule, Union):
        subbuilders = [make_default_builder(component) for component in rule.args]
        return UnionBuilder(*subbuilders)

# ------------------------------------------------------- #
# Auxilliary builders with size and choice annotations
# ------------------------------------------------------- #

# These builders are very verbose and are not meant to be used by end users.
# They are applied in the first phase of the generation of labelled structures,
# before the labelling phase and before the other builders (default or
# user-defined) are called.
#
# The produced trees are annotated with:
# - the size of every subtree with respect to each labelled atom, for later use
#   by the labelling algorithm
# - the choices made on the union rules, to pass the information to the other
#   builders

cdef size_ref_builder(id):
    def build(x):
        __, __, sizes = x
        return (REF, (id, x), sizes)
    return build

cdef size_atom_builder(labelled_atoms):
    def build(x):
        atom_name, atom_size = x
        sizes = [atom_size if name == atom_name else 0 for name in labelled_atoms]
        return (ATOM, (atom_name, atom_size), sizes)
    return build

cdef size_product_builder(builders, labelled_atoms):
    def build(terms):
        t = tuple(builders[i](terms[i]) for i in range(len(terms)))
        sizes = [0 for __ in labelled_atoms]
        for __, __, arg_sizes in t:
            for i in range(len(arg_sizes)):
                sizes[i] += arg_sizes[i]
        return (PRODUCT, t, sizes)
    return build

cdef size_union_builder(builders):
    def build(obj):
        index, content = obj
        builder = builders[index]
        content = builder(content)
        return (WRAPPED_CHOICE, (index, content), content[2])
    return build

cdef size_builder(name_to_id, labelled_atoms, rule):
    if isinstance(rule, Ref):
        return size_ref_builder(name_to_id[rule.name])
    elif isinstance(rule, Atom):
        return size_atom_builder(labelled_atoms)
    elif isinstance(rule, Union):
        subbuilders = [size_builder(name_to_id, labelled_atoms, component) for component in rule.args]
        return size_union_builder(subbuilders)
    elif isinstance(rule, Product):
        subbuilders = [size_builder(name_to_id, labelled_atoms, component) for component in rule.args]
        return size_product_builder(subbuilders, labelled_atoms)


# ------------------------------------------------------- #
# Random generation
# ------------------------------------------------------- #

# A free Boltzmann generator, but instead of generating a tree it only compute
# its size.

cdef int c_simulate(first_rule, int size_max, flat_rules, randstate rstate):
    cdef int size = 0
    cdef list todo = [first_rule]
    cdef double r = 0.

    while todo:
        type, weight, args = todo.pop()
        if type == REF:
            symbol = args
            todo.append(flat_rules[symbol])
        elif type == ATOM:
            __, atom_size = args
            size += atom_size
            if size > size_max:
                return size
        elif type == UNION:
            r = rstate.c_rand_double() * weight
            for arg in args:
                __, arg_weight, __ = arg
                r -= arg_weight
                if r <= 0:
                    todo.append(arg)
                    break
        elif type == PRODUCT:
            todo += args[::-1]
        elif type == SEQ:
            __, arg_weight, __ = args
            # k ~ Geom(arg_weight)
            k = int(ceil(log(rstate.c_rand_double())/log(arg_weight)-1))
            for __ in range(k):
                todo.append(args)
    return size

# Free Boltzmann sampler (actual generation)
# Does not handle the labelling of labelled structures

cdef c_generate(first_rule, rules, builders, randstate rstate):
    cdef list generated = []
    cdef list todo = [first_rule]
    cdef double r = 0.

    while todo:
        type, weight, args = todo.pop()
        if type == REF:
            symbol = args
            todo.append((FUNCTION, weight, symbol))
            todo.append(rules[symbol])
        elif type == ATOM:
            (name, size) = args
            generated.append((name, size))
        elif type == UNION:
            r = rstate.c_rand_double() * weight
            for i in range(len(args)):
                arg = args[i]
                __, arg_weight, __ = arg
                r -= arg_weight
                if r <= 0:
                    todo.append((WRAP_CHOICE, arg_weight, i))
                    todo.append(arg)
                    break
        elif type == PRODUCT:
            nargs = len(args)
            todo.append((TUPLE, weight, nargs))
            todo += args[::-1]
        elif type == SEQ:
            __, arg_weight, __ = args
            # k ~ Geom(arg_weight)
            k = int(ceil(log(rstate.c_rand_double())/log(arg_weight)-1))
            todo.append((LIST, weight, k))
            for __ in range(k):
                todo.append(args)
        elif type == TUPLE:
            nargs = args
            t = tuple(generated[-nargs:])
            generated = generated[:-nargs]
            generated.append(t)
        elif type == LIST:
            nargs = args
            if nargs == 0:
                generated.append([])
            else:
                t = list(generated[-nargs:])
                generated = generated[:-nargs]
                generated.append(t)
        elif type == FUNCTION:
            func = builders[args]
            x = generated.pop()
            generated.append(func(x))
        elif type == WRAP_CHOICE:
            choice = generated.pop()
            choice_number = args
            generated.append((choice_number, choice))

    obj, = generated
    return obj


# ------------------------------------------------------- #
# Labellings
# ------------------------------------------------------- #

# A few helpers first

cdef int position(str s, array):
    for i in range(len(array)):
        if s == array[i]:
            return i
    return -1

cdef int rand_int(int bound, randstate rstate):
    """Generate a uniform random integer in [0; bound[."""
    cdef int r = rstate.c_random()
    cdef int v = r % bound
    if r - v > SAGE_RAND_MAX - bound + 1:
        return rand_int(bound, rstate)
    else:
        return v

cdef int select(int r, sizes):
    """Compute the smallest i s.t. `r < sum_{j <= i} sizes[j]`."""
    cdef int t = r
    for i in range(len(sizes)):
        t -= sizes[i]
        if t < 0:
            return i

cdef gen_one_partition(bag, sizes, total_size, rstate):
    """Return a uniform partition of the list ``bag`` where the sizes of the
    elements of the partition are specified in ``sizes``."""
    total_size = total_size
    bags = [[] for __ in sizes]
    for total_size in range(total_size, 0, -1):
        r = rand_int(total_size, rstate)
        i = select(r, sizes)
        sizes[i] -= 1
        label = bag.pop()
        bags[i].append(label)
    return bags

cdef gen_partitions(bags, args_sizes, total_sizes, rstate):
    transposed = [
        gen_one_partition(
            bags[i],
            [arg_size[i] for arg_size in args_sizes],
            total_sizes[i],
            rstate
        )
        for i in range(len(bags))
    ]
    return [
        [transposed[i][j] for i in range(len(bags))]
        for j in range(len(args_sizes))
    ]

cdef shuffle_array(array, randstate rstate):
    cdef int l = len(array)
    cdef int tmp
    for i in range(l, 0, -1):
        j = rand_int(i, rstate)
        tmp = array[j]
        array[j] = array[i - 1]
        array[i - 1] = tmp

cdef tree_sizes(tree):
    __, __, sizes = tree
    return sizes

# The labelling function
# Sample a uniform labelling or a structure and apply the (default or
# user-defined) builders to the resulting tree.

cdef labelling(tree, builders, labelled_atoms, randstate rstate):
    sizes = tree_sizes(tree)
    bags = [list(range(size)) for size in sizes]
    generated = []
    todo = [(tree, bags)]

    while todo:
        (type, content, sizes), bags = todo.pop()
        if type == REF:
            rule_id, children = content
            todo.append(((FUNCTION, rule_id, sizes), None))
            todo.append((children, bags))
        elif type == ATOM:
            (name, size) = content
            i = position(name, labelled_atoms)
            if i >= 0:
                # Label!
                bag = bags[i]
                assert len(bag) == size
                if len(bag) > 1:
                    shuffle_array(bag, rstate)
                generated.append((name, bag))
            else:
                generated.append((name, size))
        elif type == PRODUCT:
            args_sizes = [tree_sizes(t) for t in content]
            args_bags = gen_partitions(bags, args_sizes, sizes, rstate)
            todo.append(((TUPLE, len(content), size), None))
            for i in range(len(content) - 1, -1, -1):
                todo.append((content[i], args_bags[i]))
        elif type == TUPLE:
            nargs = content
            t = tuple(generated[-nargs:])
            generated = generated[:-nargs]
            generated.append(t)
        elif type == FUNCTION:
            func = builders[content]
            x = generated.pop()
            generated.append(func(x))
        elif type == WRAPPED_CHOICE:
            (index, tree) = content
            todo.append(((WRAP_CHOICE, index, size), None))
            todo.append((tree, bags))
        elif type == WRAP_CHOICE:
            index = content
            x = generated.pop()
            generated.append((index, x))

    tree, = generated
    return tree


# ------------------------------------------------------- #
# High level interface
# ------------------------------------------------------- #

cdef c_gen(first_rule, rules, int size_min, int size_max, int max_retry, builders, labelled_atoms, size_builders=None):
    """Search for a tree in a given size window."""
    cdef int nb_rejections = 0
    cdef int cumulative_rejected_size = 0
    cdef int size = -1
    # A handle on the random generator
    cdef randstate rstate = current_randstate()
    # Allocate a gmp_randstate_t
    cdef gmp_randstate_t gmp_state
    gmp_randinit_default(gmp_state)

    while nb_rejections < max_retry:
        # save the random generator's state
        gmp_randinit_set(gmp_state, rstate.gmp_state)
        size = c_simulate(first_rule, size_max, rules, rstate)
        if size <= size_max and size >= size_min:
            break
        else:
            cumulative_rejected_size += size
            nb_rejections += 1

    statistics = {
        "size": size,
        "nb_rejections": nb_rejections,
        "cumulative_rejected_size": cumulative_rejected_size,
    }
    if not(size <= size_max and size >= size_min):
        return statistics, None

    # Reset the random generator to the state it was just before the simulation
    gmp_randinit_set(rstate.gmp_state, gmp_state)
    if size_builders is not None:
        assert len(labelled_atoms) > 0
        # Labelled grammar
        obj = c_generate(first_rule, rules, size_builders, rstate)
        __, __, id = first_rule
        __, __, sizes = obj
        tree = (REF, (id, obj), sizes)
        obj = labelling(tree, builders, labelled_atoms, rstate)
    else:
        obj = c_generate(first_rule, rules, builders, rstate)
    return statistics, obj

class Generator:
    """High level interface for Boltzmann samplers."""

    def __init__(self, grammar, oracle=None):
        """Make a Generator out of a grammar.

        INPUT:

        - ``grammar`` -- a combinatorial grammar

        - ``oracle`` (default: None) -- an oracle for the grammar. If not
          supplied, a default generic oracle is automatically generated.
        """
        # Load the default oracle if none is supplied
        if oracle is None:
            oracle = SimpleOracle(grammar)
        self.oracle = oracle
        # Store grammar-related information
        self.grammar = grammar
        self.atoms = sorted([
            atom.name
            for atom in grammar.atoms()
            if atom.size != 0
        ])
        self.labelled_atoms = [
            atom.name
            for atom in grammar.atoms()
            if atom.size != 0 and atom.labelled()
        ]
        # Replace all symbols in the grammar by an integer identifier
        # Use arrays rather than dictionaries
        name_to_id, id_to_name = _map_all_names_to_ids(grammar.rules)
        self.name_to_id = name_to_id
        self.id_to_name = id_to_name
        # Cached generating function values
        self.oracle_cache = {}
        self.singularity = None
        # init builders
        self.builders = [
            make_default_builder(self.grammar.rules[self.id_to_name[id]])
            for id in range(len(self.id_to_name))
        ]
        if self.grammar.labelled():
            self.size_builders = [
                size_builder(
                    self.name_to_id,
                    self.labelled_atoms,
                    self.grammar.rules[self.id_to_name[id]]
                )
                for id in range(len(self.id_to_name))
            ]

    def set_builder(self, non_terminal, func):
        """Set the builder for a non-terminal symbol.

        INPUT:

        - ``non_terminal`` -- string, the name of the non-terminal symbol

        - ``func`` -- function, the builder
        """
        symbol_id = self.name_to_id[non_terminal]
        self.builders[symbol_id] = func

    def get_builder(self, non_terminal):
        """Retrieve the current builder for a non-terminal symbol.

        INPUT:

        - ``non_terminal`` -- string, the name of the non-terminal symbol
        """
        symbol_id = self.name_to_id[non_terminal]
        return self.builders[symbol_id]

    def _precompute_oracle_values(self, z):
        variables = {self.atoms[i]: z[i] for i in range(len(z))}
        return [
            self.oracle.eval_rule(self.id_to_name[rule_id], variables)
            for rule_id in range(len(self.id_to_name))
        ]

    def gen(self, name, window, max_retry=2000, singular=True):
        """Generate a term of the grammar in a given size window.

        INPUT:

        - ``name`` -- string, the name of the symbol of the grammar you want to
          generate

        - ``window`` -- pair of integers, the size of the generated object will
          be greater than the first component of the window and lower than the
          second component

        - ``max_retry`` (default: 2000) -- integer, maximum number of attempts.
          If no object in the size window is found in less that ``max_retry``
          attempts, the generator returns None

        - ``singular`` (default: True) -- boolean, whether the generator should
          do singular sampling (the parameter z is the singularity of the
          generating function) or guess the best value for the parameter z in
          order to generate objects in the given size window.
        """
        size_min, size_max = window
        # Precompute the generating functions of each non-terminal in the
        # system
        z = None
        if singular:
            if self.singularity is None:
                values = self.oracle.find_singularity()
                self.singularity = tuple(values[name] for name in self.atoms)
                z = self.singularity
                self.oracle_cache[z] = values
            z = self.singularity
        else:
            raise NotImplementedError("Non-singular generation")
        if z not in self.oracle_cache:
            self.oracle_cache[z] = self._precompute_oracle_values(z)
        # Convert the grammar into an array of tuples (for performance)
        flat_rules = _map_all_names_to_ids_system(
            self.name_to_id,
            self.id_to_name,
            self.oracle_cache[z],
            self.grammar.rules
        )
        # Generate
        first_rule = (REF, self.oracle_cache[z][name], self.name_to_id[name])
        statistics, obj = c_gen(
            first_rule,
            flat_rules,
            size_min,
            size_max,
            max_retry,
            self.builders,
            self.labelled_atoms,
            size_builders=(self.size_builders if self.grammar.labelled() else None)
        )
        return obj, statistics
