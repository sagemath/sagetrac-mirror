# coding: utf-8
# cython: profile=True

from sage.libs.gmp.random cimport gmp_randinit_set, gmp_randinit_default
from sage.libs.gmp.types cimport gmp_randstate_t
from sage.misc.randstate cimport randstate, current_randstate
from .grammar import Atom, Product, Ref, Union

# ---
# Preprocessing
# ---

ctypedef enum options:
    REF,
    ATOM,
    UNION,
    PRODUCT,
    TUPLE,
    FUNCTION,
    WRAP_CHOICE

cdef inline wrap_choice(float weight, int id):
    return (WRAP_CHOICE, weight, id)

cdef map_all_names_to_ids_expr(name_to_id, expr):
    if isinstance(expr, Ref):
        return (REF, expr.weight, name_to_id[expr.name])
    elif isinstance(expr, Atom):
        return (ATOM, expr.weight, (expr.name, expr.size))
    elif isinstance(expr, Union):
        args = tuple((map_all_names_to_ids_expr(name_to_id, arg) for arg in expr.args))
        return (UNION, expr.weight, args)
    elif isinstance(expr, Product):
        args = tuple((map_all_names_to_ids_expr(name_to_id, arg) for arg in expr.args))
        return (PRODUCT, expr.weight, args)

cdef map_all_names_to_ids(rules):
    name_to_id = {}
    id_to_name = {}
    for i, name in enumerate(rules.keys()):
        name_to_id[name] = i
        id_to_name[i] = name
    rules = [
        map_all_names_to_ids_expr(name_to_id, rules[id_to_name[i]])
        for i in range(len(name_to_id))
    ]
    return name_to_id, id_to_name, rules

# ---
# Simulation phase
# ---

cdef int c_simulate(int id, float weight, int size_max, flat_rules, randstate rstate):
    cdef int size = 0
    cdef list todo = [(REF, weight, id)]
    cdef float r = 0.

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

    return size

# ---
# Actual generation
# ---

cdef c_generate(int id, float weight, flat_rules, builders, randstate rstate):
    generated = []
    cdef list todo = [(REF, weight, id)]
    cdef float r = 0.
    while todo:
        type, weight, args = todo.pop()
        if type == REF:
            symbol = args
            todo.append((FUNCTION, weight, symbol))
            todo.append(flat_rules[symbol])
        elif type == ATOM:
            atom_name, __ = args
            generated.append(atom_name)
        elif type == UNION:
            r = rstate.c_rand_double() * weight
            for i in range(len(args)):
                arg = args[i]
                __, arg_weight, __ = arg
                r -= arg_weight
                if r <= 0:
                    todo.append(wrap_choice(arg_weight, i))
                    todo.append(arg)
                    break
        elif type == PRODUCT:
            nargs = len(args)
            todo.append((TUPLE, weight, nargs))
            todo += args[::-1]
        elif type == TUPLE:
            nargs = args
            t = tuple(generated[-nargs:])
            generated = generated[:-nargs]
            generated.append(t)
        elif type == FUNCTION:
            func = builders[args]
            x = generated.pop()
            generated.append(func(x))
        elif type == WRAP_CHOICE:
            choice = generated.pop()
            generated.append((args, choice))

    obj, = generated
    return obj


cdef c_gen(int id, float weight, flat_rules, int size_min, int size_max, int max_retry, builders):
    cdef int nb_rejections = 0
    cdef int cumulative_rejected_size = 0
    cdef int size = -1
    cdef randstate rstate = current_randstate()
    cdef gmp_randstate_t gmp_state
    gmp_randinit_default(gmp_state)

    while nb_rejections < max_retry:
        gmp_randinit_set(gmp_state, rstate.gmp_state)
        size = c_simulate(id, weight, size_max, flat_rules, rstate)
        if size <= size_max and size >= size_min:
            break
        else:
            cumulative_rejected_size += size
            nb_rejections += 1

    if not(size <= size_max and size >= size_min):
        return None

    gmp_randinit_set(rstate.gmp_state, gmp_state)
    obj = c_generate(id, weight, flat_rules, builders, rstate)
    statistics = {
        "size": size,
        "nb_rejections": nb_rejections,
        "cumulative_rejected_size": cumulative_rejected_size,
    }
    return statistics, obj


cdef inline identity(x):
    return x

def UnionBuilder(*builders):
    def build(obj):
        index, content = obj
        builder = builders[index]
        return builder(content)
    return build

cdef inline ProductBuilder(builders):
    def build(terms):
        return tuple(builders[i](terms[i]) for i in range(len(terms)))
    return build

cdef make_default_builder(rule):
    type, __, args = rule
    if type == REF:
        return identity
    elif type == ATOM:
        return identity
    elif type == UNION:
        subbuilders = [make_default_builder(component) for component in args]
        return UnionBuilder(*subbuilders)
    elif type == PRODUCT:
        subbuilders = [make_default_builder(component) for component in args]
        return ProductBuilder(subbuilders)


class Generator:
    def __init__(self, grammar, oracle=None):
        # Load the default oracle if none is supplied
        if oracle is None:
            from .oracle import OracleSimple
            oracle = OracleSimple(grammar, e1=1e-6, e2=1e-6)
        self.oracle = oracle
        # flatten the grammar for faster access to rules
        self.grammar = grammar
        self.grammar.annotate(oracle)
        name_to_id, id_to_name, flat_rules = map_all_names_to_ids(grammar.rules)
        self.name_to_id = name_to_id
        self.id_to_name = id_to_name
        self.flat_rules = flat_rules
        # init builders
        self.builders = [
            make_default_builder(rule)
            for rule in self.flat_rules
        ]

    def set_builder(self, non_terminal, func):
        symbol_id = self.name_to_id[non_terminal]
        self.builders[symbol_id] = func

    def get_builder(self, non_terminal):
        symbol_id = self.name_to_id[non_terminal]
        return self.builders[symbol_id]

    def gen(self, name, window, max_retry=2000):
        id = self.name_to_id[name]
        weight = self.grammar.rules[name].weight
        size_min, size_max = window
        statistics, obj = c_gen(
            id,
            weight,
            self.flat_rules,
            size_min,
            size_max,
            max_retry,
            self.builders,
        )
        return obj, statistics
