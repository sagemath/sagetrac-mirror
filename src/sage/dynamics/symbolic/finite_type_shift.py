# -*- coding: utf-8 -*-
r"""
One-dimensional Subshifts of Finite Type.

This ``.py`` file contains the two SAGE classes ``SFT`` and ``sfts``. The
first one is an implementation of one-dimensional shifts of finite type which
provides methods for computing invariants, producing different (conjugate)
representations etc. useful in the study of these symbolic dynamical systems.
The second class provides a collection of predefined examples and families of
well-known shifts of finite type which can be used immediately after calling
their constructor methods with all functionalities provided by the class SFT.

AUTHORS:

- Michael Schraudner (2011):

    - initial version of class SFT: constructor method, basic functionalities
      and methods,
    - initial version of class sfts = SFTGenerators, mschraudner@dim.uchile.cl

- Sebastian Donoso (2011):

    - added some methods to class SFT, some documentation

- Sebastian Barbieri (2012):

    - added some methods and documentation to class SFT added examples and
      documentation to class sfts

- Michael Schraudner (2012-05-XX):

    - modifications throughout large parts of the code, including support for
      class Words, overall clean up to include classes SFT and sfts into SAGE

EXAMPLES:

.. linkall

We demonstrate 3 ways to create an SFT from different input objects.
    #. A directed graph (where the SFT's alphabet is extracted - depending
       on the chosen format - either from the graph's edge or the graph's
       vertex labels)::

            sage: G = DiGraph({"0":{"0":["a", "b"], "1":["c"]}, "1":{"1":["d"]}}, multiedges=True)
            sage: X = SFT(G); X
            A 1-step SFT over ['a', 'b', 'c', 'd'].

       ::

            sage: G = DiGraph({0:{0:["a"], 1:["b"]}, 1:{0:["c"]}}, multiedges=True)
            sage: X = SFT(G, format='vertex'); X
            A 1-step SFT over [0, 1].

    #. A square matrix either in vertex or in edge representation (note
       that as the alphabet is not specified, a standard alphabet is
       created)::

            sage: M = matrix(4, 4, [1,1,1,1, 0,1,0,1, 1,0,1,1, 0,0,0,1])
            sage: X = SFT(M, format='vertex'); X
            A 1-step SFT over [0, 1, 2, 3].

       ::

            sage: M = matrix(2, 2, [1,2, 3,0])
            sage: X = SFT(M, format='edge', alphabet=[str(i) for i in range(6)]); X
            A 1-step SFT over ['0', '1', '2', '3', '4', '5'].

    #. A list of forbidden words (note that the word '010' is implicitly
       forbidden as well, so the order of the SFT is 2 and not 3)::

            sage: X = SFT(["0101", "100"], alphabet=["0", "1"]); X
            A 2-step SFT over ['0', '1'].

       ::

            sage: X = SFT([Word("0101"), Word(["1", "0", "0"])]); X
            A 2-step SFT over ['0', '1'].

    #. A list of forbidden words in list format (containing a forbidden
       symbol)::

            sage: X = SFT([['1', '1'], ['2']], alphabet=["0", "1", "2"], name='Fibonacci shift'); X
            The Fibonacci shift of order 1 over ['0', '1'].

    #. A list of forbidden words over an alphabet with symbols of distinct
       lengths using a simple symbol separator::

            sage: X = SFT(['1.1'], ['0', '1'], name='Fibonacci shift', symb_sep="."); X
            The Fibonacci shift of order 1 over ['0', '1'].

       ::

            sage: X = SFT(['1.11', '11.0.11.00'], symb_sep="."); X
            A 3-step SFT over ['0', '00', '1', '11'].
"""

#*****************************************************************************
#     Copyright (C) 2011, 2012 Michael Schraudner <mschraudner@dim.uchile.cl>
#     Copyright (C) 2011 Sebastian Donoso
#     Copyright (C) 2012 Sebastian Barbieri
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.calculus.var import var
from sage.combinat.words.word import Word
from sage.functions.log import log
from sage.functions.other import sqrt, floor
from sage.graphs.digraph import DiGraph
from sage.groups.additive_abelian.additive_abelian_group import AdditiveAbelianGroup
from sage.matrix.constructor import matrix
from sage.misc.prandom import randint
from sage.rings.all import ZZ, NN, RR, Infinity
from sage.arith.all import gcd, divisors, moebius
from sage.rings.polynomial.polynomial_element import Polynomial
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.sets.set import Set
from sage.structure.element import is_Matrix
from sage.structure.sage_object import SageObject
from sage.symbolic.expression import Expression
from sage.combinat.backtrack import SearchForest

import warnings
from copy import copy

def sft_warning_style(msg, category, filename, lineno, file=None, line=None):
    return "%s: %s\n" % (category.__name__, msg)

# warnings.formatwarning = sft_warning_style
warnings.filterwarnings(action="always",
                        module="sage.dynamics.symbolic.finite_type_shift")


class SFT(SageObject):
    r"""
    A class for one-dimensional Subshifts of Finite Type (SFTs).

    For more information on one-dimensional SFTs, see the Wikipedia article
    `<http://en.wikipedia.org/wiki/Subshift_of_finite_type>`_
    containing some background information and further references on SFTs or
    take a look at the books mentioned in the REFERENCES section below.

    A Subshift of Finite Type can be created from different input objects
    as specified below. However the simplest way to generate an (empty) SFT in
    SAGE is by typing::

        sage: X = SFT()

    By typing the object's name, one can get some basic information about it::

        sage: X
        The empty shift.

    This first example is not very interesting as it is by default the empty
    shift. In the EXAMPLES section below we will explain how to define general
    SFTs calling the class' generator method with different input. SAGE
    however also contains a library of several pre-defined families of SFT
    examples in its class ``sfts`` which can be accessed in the following way:

    - Within a SAGE Session, type ``sfts.`` (do not press 'Enter', and do not
      forget the final period '.')
    - Hit the 'Tab' key

    You will see a list of commands which will construct some named and
    well-known SFTs. For example::

        sage: F = sfts.Fibonacci()
        sage: F.draw_graph()

    In order to obtain more information about these SFT-generators, access
    their documentation.

    Once you have defined the SFT you want, you can begin to work with it by
    using the methods predefined in the SAGE class ``SFT``. As usual a list of
    all available methods con be obtained in the following way:

    - Within a SAGE session, type the SFT's name followed by a period or type
      ``SFT.`` (do not press 'Enter', and do not forget the final period '.')
    - Hit the 'Tab' key

    As usual some additional information about what these methods do can be
    obtained by typing the method's name followed by a question mark. E.g. if
    you want to know about the entropy()-method just type: SFT.entropy?

    In order to use this (and the other) method(s) on a given SFT - which we
    will assume is called X - type and execute X.entropy() in order to compute
    the SFTs (topological) entropy.

    INPUT:

    - ``init_obj`` - (default: None) can be any of the following:

        - a (finite) directed (multi)graph
        - a square matrix with non-negative integer entries
        - a list of forbidden words (over some alphabet) either given as:

            - a list of strings (with or without symbol separators), e.g.
              ['11', '123', '4'] (no separator) or ['1.1', '1.2.3', '4']
              ('.' used as separator)
            - a list of lists containing elements of the alphabet, e.g.
              [[1, 1], [1, 2, 3], [4]] or [['a', 'a'], ['a', 'b', 'b']]
            - a list of instances of the SAGE class Words (over the alphabet),
              e.g. [Word("11"), Word(["1", "2", "3"]), Word("4")]

      If this parameter is not provided the full-shift over the specified
      alphabet will be created. If neither this parameter nor an alphabet
      is given the empty shift will be created.

    - ``alphabet`` - (default: None) a list of symbols (the alphabet) used in
      the definition of the SFT. If this parameter is not provided a standard
      alphabet will be created (whenever possible) from the ``init_obj``. The
      elements of the alphabet can be general SAGE objects which do not
      contain the substring ``symb_sep`` in their string representation.
      Nevertheless it is preferable to use symbols which are either
      non-negative integers or short strings (e.g. single letters).

      The alphabet is not allowed to contain repetitions (literal copies as
      well as symbols whose string representations do coincide) and whenever
      its elements have distinct lengths the use of ``symb_sep`` is strongly
      advised. We also highly recommend the use of unambiguous alphabets!

    - ``format`` - (default: 'edge') this flag selects one of two possible
      representations for the SFT:

        - 'vertex' - a representation of the SFT as a vertex shift
        - 'edge' - a representation of the SFT as an edge shift

      If this parameter is not provided in the constructor the 'edge'
      representation is chosen by default.

      Remember that the 'vertex' representation requires the transition matrix
      to be 0/1 resp. the digraph not to contain multiple edges.

    - ``symb_sep``- (default: None) a string (or None) which is used to
      separate distinct symbols in the forbidden words (given as strings). The
      use of this separator is necessary in case the alphabet contains symbols
      of distinct lengths and it will then also be used as a separator in the
      graph presentations (vertex and edge labels).

    - ``name`` - (default: None) a string specifying a particular name for the
      SFT. If this parameter is not provided the subshift will have no
      specific name.

    OUTPUT:

    - a Subshift of Finite Type

    REFERENCES:

    - 'D. Lind, B. Marcus: An introduction to symbolic dynamics and coding.'
    - 'B. Kitchens: Symbolic dynamics. One-sided, two-sided and countable state Markov shifts.'

    EXAMPLES:

    We demonstrate 3 ways to create an SFT from different input objects.
        #. A directed graph (where the SFT's alphabet is extracted - depending
           on the chosen format - either from the graph's edge or the graph's
           vertex labels)::

                sage: G = DiGraph({"0":{"0":["a", "b"], "1":["c"]}, "1":{"1":["d"]}}, multiedges=True)
                sage: X = SFT(G); X
                A 1-step SFT over ['a', 'b', 'c', 'd'].

           ::

                sage: G = DiGraph({0:{0:["a"], 1:["b"]}, 1:{0:["c"]}}, multiedges=True)
                sage: X = SFT(G, format='vertex'); X
                A 1-step SFT over [0, 1].

           Note that for a vertex representation the graph may not contain any
           multi-edges.

        #. A square matrix either in vertex or in edge representation (note
           that as the alphabet is not specified, a standard alphabet is
           created)::

                sage: M = matrix(4, 4, [1,1,1,1, 0,1,0,1, 1,0,1,1, 0,0,0,1])
                sage: X = SFT(M, format='vertex'); X
                A 1-step SFT over [0, 1, 2, 3].

           ::

                sage: M = matrix(2, 2, [1,2, 3,0])
                sage: X = SFT(M, format='edge', alphabet=[str(i) for i in range(6)]); X
                A 1-step SFT over ['0', '1', '2', '3', '4', '5'].

           Note that to produce a vertex representation the matrix has to be
           defined over 0 and 1 only while for an edge representation it can
           contain non-negative integer entries.

        #. A list of forbidden words either in string or in Word format (note
           that the word "010" is implicitly forbidden as well, hence the
           order of this example SFT is 2 and not 3)::

                sage: X = SFT(["0101", "100"], alphabet=["0", "1"]); X
                A 2-step SFT over ['0', '1'].

           ::

                sage: X = SFT([Word("0101"), Word(["1", "0", "0"])]); X
                A 2-step SFT over ['0', '1'].

        #. A list of forbidden words in list format (containing a forbidden
           symbol); this is an unusual way to define the well known Fibonacci
           shift from scratch (i.e. without using its predefined version in
           class ``sfts``)::

                sage: X = SFT([['1', '1'], ['2']], alphabet=["0", "1", "2"], name='Fibonacci shift'); X
                The Fibonacci shift of order 1 over ['0', '1'].

        #. A list of forbidden words over an alphabet with symbols of distinct
           lengths using a simple symbol separator::

                sage: X = SFT(['1.1'], ['0', '1'], name='Fibonacci shift', symb_sep="."); X
                The Fibonacci shift of order 1 over ['0', '1'].

           ::

                sage: X = SFT(['1:11', '11:0:11:00'], symb_sep=":"); X
                A 3-step SFT over ['0', '00', '1', '11'].

           ::

                sage: X = SFT(['one two', 'two zero two three'], symb_sep=" "); X
                A 3-step SFT over ['one', 'three', 'two', 'zero'].

        #. To avoid confusion we strongly recommend to define SFTs over
           ambiguous alphabets only by using either forbidden words being
           lists or Words or by using an (otherwise not appearing) symbol
           separator::

                sage: X = SFT(['12'], ['1', '2', '12']); X
                Traceback (most recent call last):
                ...
                ValueError: need a valid symbol separator symb_sep

                sage: X = SFT(['1.2'], ['1', '2', '12'], symb_sep="."); X
                A 1-step SFT over ['1', '12', '2'].

    AUTHORS:

    - Michael Schraudner (2011): initial version of class SFT: constructor
      method, basic functionalities and methods, mschraudner@dim.uchile.cl

    - Sebastian Donoso (2011): added methods to class SFT, basic documentation

    - Sebastian Barbieri (2012): added methods and documentation to class SFT

    - Michael Schraudner (2012-05-XX): modifications throughout large parts of
      the code, including support for class Words, overall clean up to include
      class SFT into SAGE
    """

    def __init__(self, init_obj=None, alphabet=None, format='edge',
                 symb_sep=None, name=None):
        r"""
        The __init__ method for creating an SFT.

        For instructions on how to use this method type: SFT?

        Private attributes defined inside the class:

        - self._alph:       alphabet of the SFT (list)
        - self._format:     format of the SFT (i.e. 'vertex' or 'edge')
        - self._fwords:     forbidden words of the SFT (list)
        - self._graph:      DiGraph of the SFT (if the SFT has order N>1 this
                           object will hold the smallest possible higher
                           block presentation, i.e. the N resp. N-1 higher
                           graph depending on the format being 'vertex' resp.
                           'edge')
        - self._is_empty:   'True' iff the SFT is empty (boolean)
        - self._matrix:     matrix defining the SFT according to the format;
                           either a 0/1 matrix ('vertex') or a matrix over the
                           non-negative integers ('edge')
        - self._name:       name of the SFT (string)
        - self._obj:        initial object the SFT was generated from
        - self._symb_sep:   separator of symbols in forbidden words (i.e. if
                           the alphabet has symbols of distinct lengths the
                           value of _symb_sep is used to separate symbols in
                           the forbidden words given as strings) and in the
                           graph presentation (string)

        AUTHORS:

        - Michael Schraudner (2011): initial version, mschraudner@dim.uchile.cl
        - Michael Schraudner (2012): implemented support for SAGE class Words

        TESTS::

            sage: M = matrix(2, 2, [2,1, 1,0])
            sage: SFT(M, format='vertex')
            doctest:...: UserWarning: matrix must be 0/1 to create a vertex shift. Will try to create an edge shift instead.
            A 1-step SFT over [0, 1, 2, 3].
            sage: G = DiGraph({0:{0:['a'], 1:['b','c']}, 1:{0:['d']}}, multiedges=True)
            sage: SFT(G, format='vertex')
            doctest:...: UserWarning: digraph must not have multiple edges to create a vertex shift. Will try to create an edge shift instead.
            A 1-step SFT over ['a', 'b', 'c', 'd'].

        ::

            sage: SFT([Word("1.0")], symb_sep=".")
            Traceback (most recent call last):
            ...
            ValueError: forbidden words contain the symbol separator symb_sep

        ::

            sage: SFT([[1, 0, 1]], alphabet=[2, 0], symb_sep='.')
            doctest:...: UserWarning: forbidden words contain symbols not in the given alphabet. Will enlarge alphabet.
            A 2-step SFT over [0, 1, 2].

        ::

            sage: SFT([Word("11"),[1,1]])
            Traceback (most recent call last):
            ...
            ValueError: forbidden words contain ambiguous symbols
        """
        if name is None or isinstance(name, basestring):
            self._name = name
        else:
            raise TypeError("name must be a string or None")
        if symb_sep is None:
            self._symb_sep = ""
        elif isinstance(symb_sep, basestring):
            self._symb_sep = symb_sep
        else:
            raise TypeError("symb_sep must be a string or None")
        if format == 'edge' or format == 'vertex':
            self._format = format
        else:
            raise ValueError("format not valid; must be 'vertex' or 'edge'")
        if alphabet is None:
            self._alph = []
        else: # force it to be a list (+ print a warning)
            try:
                self._alph = list(alphabet)
            except:
                warnings.warn("alphabet must be a list. Will create a " +
                              "standard alphabet.")
                self._alph = []
        self._check_alph(0)
        self._obj = init_obj # keep initial object for further reference

        ## initialize with a matrix ##
        if is_Matrix(init_obj):
            if not init_obj.is_square():
                raise TypeError("matrix must be square")
            if any([x not in NN for x in init_obj.list()]):
                raise ValueError("matrix must have non-negative integer " +
                                 "entries")
            if all([x == 0 for x in init_obj.eigenvalues()]):
                self._create_empty() # create empty shift if matrix is nilpot
            else:
                self._is_empty = False
                self._matrix = init_obj
                ## vertex shift ##
                if self._format == 'vertex':
                    if any([x not in [0,1] for x in init_obj.list()]):
                        warnings.warn("matrix must be 0/1 to create a " +
                                      "vertex shift. Will try to create an " +
                                      "edge shift instead.")
                        self._format = 'edge'
                    else:
                        if self._alph == []: # create standard alphabet
                            self._alph = range(init_obj.ncols())
                            if init_obj.ncols() > 10:
                                self._symb_sep = '.'
                        elif not len(self._alph) == init_obj.ncols():
                            warnings.warn("alphabet has wrong size. Will " +
                                          "create standard alphabet.")
                            self._alph = range(init_obj.ncols())
                            if init_obj.ncols() > 10:
                                self._symb_sep = '.'
                        self._create_vertex_label_DiGraph()
                ## edge shift ##
                if self._format == 'edge': # HAVE to use new if (not an else)
                    n = sum(sum(init_obj))
                    if self._alph == []: # create standard alphabet
                        self._alph = range(n)
                        if n > 10:
                            self._symb_sep = '.'
                    elif not len(self._alph) == n:
                        warnings.warn("alphabet has wrong size. Will create" +
                                      " standard alphabet.")
                        self._alph = range(n)
                        if n > 10:
                            self._symb_sep = '.'
                    self._create_edge_label_DiGraph()

        ## initialize with a digraph ##
        elif isinstance(init_obj, DiGraph):
            self._matrix = init_obj.adjacency_matrix()
            if all([x == 0 for x in self._matrix.eigenvalues()]):
                self._create_empty() # create empty shift if matrix is nilpot
            else:
                self._is_empty = False
                ## vertex shift ##
                if self._format == 'vertex':
                    if init_obj.has_multiple_edges():
                        warnings.warn("digraph must not have multiple edges" +
                                      " to create a vertex shift. Will try " +
                                      "to create an edge shift instead.")
                        self._format = 'edge'
                    else:
                        if self._alph == []: # if no alphabet given, create it
                            self._alph = init_obj.vertices()
                            self._check_alph(init_obj.num_verts())
                        elif not len(self._alph) == init_obj.num_verts():
                            warnings.warn("alphabet has wrong size. Will " +
                                          "create an alphabet from " +
                                          "vertices.")
                            self._alph = init_obj.vertices()
                            self._check_alph(init_obj.num_verts())
                        self._create_vertex_label_DiGraph()
                ## edge shift ##
                if self._format == 'edge': # HAVE to use new if (not an else)
                    if self._alph == []:
                        self._alph = init_obj.edge_labels()
                        self._check_alph(init_obj.num_edges())
                    elif not len(self._alph) == init_obj.num_edges():
                        warnings.warn("alphabet has wrong size. Will create" +
                                      " an alphabet from edges.")
                        self._alph = init_obj.edge_labels()
                        self._check_alph(init_obj.num_edges())
                    self._create_edge_label_DiGraph()

        ## initialize with a list of forbidden words ##
        elif isinstance(init_obj, (list, tuple, set, frozenset)):
            self._fwords = []
            alphabet = set([])
            for fword in init_obj:
                if isinstance(fword, basestring) and not self._symb_sep == "":
                    l = fword.split(self._symb_sep)
                else:
                    try:
                        l = list(fword)
                    except:
                        raise TypeError("forbidden words must be of type " +
                                        "string, list, tuple, Word etc.")
                self._fwords.append(l)
                alphabet.update(l)
            if self._alph == []:
                self._alph = list(alphabet)
            else:
                symbs = filter(lambda a: not a in self._alph, alphabet)
                if not symbs == []:
                    warnings.warn("forbidden words contain symbols not in" +
                                  " the given alphabet. Will enlarge " +
                                  "alphabet.")
                    self._alph.extend(symbs)
            self._alph = filter(lambda a: not [a] in self._fwords, self._alph)
            self._alph.sort()
            alph_str = [str(a) for a in self._alph]
            if not (self._symb_sep == "" or
                    all([not self._symb_sep in s for s in alph_str])):
                raise ValueError("forbidden words contain the symbol " +
                                 "separator symb_sep")
            if not len(self._alph) == len(set(alph_str)):
                raise ValueError("forbidden words contain ambiguous symbols")
            if self._symb_sep == "":
                lenth = [len(s) for s in alph_str]
                if not lenth == [] and not min(lenth) == max(lenth):
                    raise ValueError("need a valid symbol separator symb_sep")
            self._fwords = filter(lambda w: len(w) > 1, self._fwords)
            self._fwords = filter(lambda w: [a for a in w
                                  if not a in self._alph] == [], self._fwords)
            self._create_fwords_DiGraph()

        ## initialize the empty shift ##
        elif init_obj is None:
            self._create_empty()

        ## if init_obj is none of the above ##
        else:
            raise NotImplementedError("so far SFTs can only be created by a" +
                                      " list of forbidden words, a matrix " +
                                      "or a digraph.")

## __METHODS__ ##

    def _repr_(self):
        r"""
        A string representation of the SFT.

        OUTPUT:

        - the representation of the SFT as a string

        EXAMPLES:

        #. The empty shift::

            sage: X = SFT(); X
            The empty shift.

        #. An unnamed example SFT::

            sage: X = SFT(["0101", "100"], alphabet=["0", "1"]); X
            A 2-step SFT over ['0', '1'].

        #. A named well-known example SFT::

            sage: X = SFT(['1.1'], ['0', '1'], name='Fibonacci shift', symb_sep="."); X
            The Fibonacci shift of order 1 over ['0', '1'].

        AUTHORS:

        - Sebastian Donoso (2011): initial version
        """
        if self._is_empty:
            return "The empty shift."
        else:
            if self._name is None:
                return "A %d-step SFT over %s." % (self.order(), self._alph)
            else:
                return "The %s of order %d over %s." % (self._name, self.order(), self._alph)

    def __latex__(self):
        r"""
        A string representation of the SFT in LaTeX format.

        OUTPUT:

        - the LaTeX representation of the SFT

        EXAMPLES:

        #. The empty shift::

            sage: X = SFT(); X
            The empty shift.

        #. An unnamed example SFT::

            sage: X = SFT(["0101", "100"], alphabet=["0", "1"]); X
            A 2-step SFT over ['0', '1'].

        #. A named well-known example SFT::

            sage: X = SFT(['1.1'], ['0', '1'], name='Fibonacci shift', symb_sep="."); X
            The Fibonacci shift of order 1 over ['0', '1'].

        AUTHORS:

        - Sebastian Barbieri (2012): initial version
        """
        if self._is_empty:
            return "The empty shift."
        else:
            if self._name is None:
                return "A %d-step SFT over %s." % (self.order(), self._alph)
            else:
                return "The %s of order %d over %s." % (self._name, self.order(), self._alph)

## METHODS ##

    def alphabet(self):
        r"""
        Returns the alphabet associated to the SFT.

        OUTPUT:

        - a list of symbols representing the alphabet of the SFT

        EXAMPLES:

        #. Alphabets of some SFTs::

                sage: X = SFT(matrix(2, 2, [1,0, 1,1])); X.alphabet()
                [0, 1, 2]

           ::

                sage: X = sfts.Fibonacci(); X.alphabet()
                [0, 1]

           ::

                sage: X = SFT(['ab', 'bc'], ['a', 'b', 'c']); X.alphabet()
                ['a', 'b', 'c']

        AUTHORS:

        - Sebastian Donoso (2011): initial version
        """
        return self._alph

    def matrix(self):
        r"""
        Returns the transition matrix associated to the SFT.

        Note that two equivalent representations of an SFT may have different
        matrix representations. For example, the full-N-shift seen as an edge
        shift will be represented by a 1x1 matrix with the number of edges as
        its unique entry, while seen as a vertex shift the matrix will have
        size NxN and will be filled with 1s.

        OUTPUT:

        - a matrix over the non-negative integers representing the SFT

        EXAMPLES:

        #. The matrices associated to some SFTs::

                sage: X = SFT(matrix(2, 2, [2,1, 2,0])); X.matrix()
                [2 1]
                [2 0]

           ::

                sage: X = sfts.Fibonacci(); X.matrix()
                [1 1]
                [1 0]

           ::

                sage: X = SFT(['ab', 'bc'], ['a', 'b', 'c']); X.matrix()
                [1 0 1]
                [1 1 0]
                [1 1 1]

        AUTHORS:

        - Sebastian Donoso (2011): initial version
        """
        return self._matrix

    def graph(self):
        r"""
        Returns a graph presentation associated to the SFT.

        If the SFT was defined by a transition matrix (or a graph) this method
        will return the corresponding graph. If the SFT was defined by a list
        of forbidden words it will return the higher graph presentation in
        which the vertices correspond to words of length the order of the SFT
        and the edges correspond to words of length order + 1 (i.e. the length
        of the longest forbidden word used in the definition of the SFT).

        OUTPUT:

        - a digraph representing the SFT

        EXAMPLES:

        #. Graphs associated to some SFTs::

                sage: X = SFT(matrix(2, 2, [2,1, 2,0])); X.graph()
                Looped multi-digraph on 2 vertices

           ::

                sage: X = sfts.Fibonacci(); X.graph()
                Looped multi-digraph on 2 vertices

           ::

                sage: X = SFT(['ab', 'bc'], ['a', 'b', 'c']); X.graph()
                Looped multi-digraph on 3 vertices

        AUTHORS:

        - Sebastian Donoso (2011): initial version
        """
        return self._graph

    def draw_graph(self, size=8, **kwds):
        r"""
        Draws the graph of the SFT with its irreducible components colored
        distinctly.

        INPUT:

        - ``size`` - (default: 8) a natural number specifying the size
          (figsize) of the graph presentation drawn

        OUTPUT:

        - a plot of the graph which uses different colors for distinct
          irreducible components

          Depending on the format ('vertex' resp. 'edge') of the SFT the drawn
          graph contains the labels of either its vertices or its edges. Other
          parameters can be specified via the ``**kwds`` option in the call which
          are then just passed through to the DiGraph.show() method.

        EXAMPLES:

        #. Drawing the graphs associated with some SFTs::

                sage: X = SFT(matrix(2, 2, [1,2, 1,1])); X.draw_graph()

           ::

                sage: X = sfts.Fibonacci(); X.draw_graph()

           ::

                sage: X = SFT(['ab', 'bc'], ['a', 'b', 'c']); X.draw_graph(10)

        #. Drawing the graph of an SFT with two non-trivial irreducible
           components::

            sage: M = matrix(5, 5, [1,1,0,0,0, 1,0,1,0,0, 0,0,0,1,0, 0,0,0,0,1, 0,0,1,1,0])
            sage: X = SFT(M, format='vertex')
            sage: X.draw_graph()

        AUTHORS:

        - Sebastian Donoso (2011): initial version
        """
        if not size in NN:
            raise ValueError("size has to be a natural number")
        irr_comps = self._graph.strongly_connected_components()
        if self._format == 'vertex':
            self._graph.show(vertex_labels=True, edge_labels=False,
                             partition=irr_comps, graph_border=True,
                             figsize=size, **kwds)
        else:
            self._graph.show(vertex_labels=False, edge_labels=True,
                             partition=irr_comps, graph_border=True,
                             figsize=size, **kwds)

    def forbidden_words(self, format='string'):
        r"""
        Returns the forbidden words of the SFT in the specified format.

        If the SFT was defined by a transition matrix (or a graph) this method
        will return the list of forbidden words of length 2 corresponding to
        forbidden transitions. If the SFT was defined by a list of forbidden
        words it will return the corresponding list.

        INPUT:

        - ``format`` - (default: 'string') a string defining the format in
          which the forbidden words are output; it can be one of the
          following:

            - 'string': the forbidden words are returned as a list of strings.
              If the alphabet has symbols of distinct lengths and a symbol
              separator was given when defining the SFT, symbols in the
              forbidden word strings will be separated by this separator
            - 'list': the forbidden words are returned as a list of lists
              containing alphabet symbols
            - 'Word': the forbidden words are returned as a list of SAGE Words

          If ``format`` is not specified, the format 'string' is chosen by
          default.

        OUTPUT:

        - a list of forbidden words in the output format specified by the
          parameter ``format``

        EXAMPLES:

        #. Forbidden words of some SFTs::

                sage: X = SFT(matrix(2, 2, [1,0, 1,1]), alphabet=["a", "b", "c"])
                sage: X.forbidden_words()
                ['ab', 'ac', 'bb', 'bc', 'ca']

           ::

                sage: X = sfts.Fibonacci()
                sage: X.forbidden_words()
                ['11']

           ::

                sage: X = SFT(['ab', 'bc'], ['a', 'b', 'c'])
                sage: X.forbidden_words()
                ['ab', 'bc']
                sage: X.forbidden_words(format='list')
                [['a', 'b'], ['b', 'c']]
                sage: X.forbidden_words(format='Word')
                [word: ab, word: bc]

           ::

                sage: X = SFT([[1,11], [11,2,11,22]], symb_sep=".")
                sage: X.forbidden_words()
                ['1.11', '11.2.11.22']
                sage: X.forbidden_words('Word')
                [word: 1,11, word: 11,2,11,22]

        AUTHORS:

        - Sebastian Donoso (2011): initial version
        - Sebastian Barbieri (2012): included presentation choices
        - Michael Schraudner (2012): included Word support,
          mschraudner@dim.uchile.cl
        """
        if format == 'list':
            return self._fwords
        elif format == 'Word':
            return [Word(w) for w in self._fwords]
        elif format == 'string':
            return [self._symb_sep.join([str(a) for a in w])
                                                for w in self._fwords]
        else:
            raise ValueError("format must be 'string', 'list' or 'Word'")

    def forbidden_words_iterator(self, format='string'):
        r"""
        Returns an iterator yielding the forbidden words of the SFT in the
        specified format.

        If the SFT was defined by a transition matrix (or a graph) this method
        will return an iterator over the forbidden words of length 2
        corresponding to forbidden transitions. If the SFT was defined by a
        list of forbidden words it will return the corresponding iterator.

        INPUT:

        - ``format`` - (default: 'string') a string defining the format in
          which the forbidden words are output; it can be one of the
          following:

            - 'string': the forbidden words are returned as a list of strings.
              If the alphabet has symbols of distinct lengths and a symbol
              separator was given when defining the SFT, symbols in the
              forbidden words string will be separated by this separator
            - 'list': the forbidden words are returned as a list of lists
              containing alphabet symbols
            - 'Word': the forbidden words are returned as a list of SAGE Words

          If ``format`` is not specified, the format 'string' is chosen by
          default.

        OUTPUT:

        - an iterator yielding the forbidden words of the SFT in the output
          format specified by the parameter ``format``

        EXAMPLES:

        #. Forbidden words iterator of the S-Gap shift::

            sage: X = sfts.S_Gap([1,2,4,5,7,8]); X
            The [1, 2, 4, 5, 7, 8]-Gap shift of order 9 over [0, 1].
            sage: T = X.forbidden_words_iterator()
            sage: for t in range(6): next(T)
            '101'
            '1001'
            '100001'
            '1000001'
            '100000001'
            '1000000001'
            sage: T = X.forbidden_words_iterator('Word')
            sage: for t in range(6): next(T)
            word: 101
            word: 1001
            word: 100001
            word: 1000001
            word: 100000001
            word: 1000000001

        AUTHORS:

        - Sebastian Barbieri (2012): initial version
        - Michael Schraudner (2012): different formats support,
          mschraudner@dim.uchile.cl
        """
        if format == 'list':
            for word in self._fwords:
                yield word
        elif format == 'Word':
            for word in [Word(w) for w in self._fwords]:
                yield word
        elif format == 'string':
            for word in [self._symb_sep.join([str(a) for a in w])
                                                     for w in self._fwords]:
                yield word
        else:
            raise ValueError("format must be 'string', 'list' or 'Word'")

    def get_info(self):
        r"""
        Print all the information stored in the SFT.

        OUTPUT:

        - prints the SFT's name, its alphabet and format, the matrix
          representing the SFT, a list of forbidden words and the SFT's order

        EXAMPLES:

        #. Information on the Empty shift::

            sage: SFT().get_info()
            The empty shift.

        #. Information on the Fibonacci shift::

            sage: X = sfts.Fibonacci()
            sage: X.get_info()
            Name: Fibonacci shift
            Alphabet: [0, 1]
            Format: vertex
            A representing matrix:
            [1 1]
            [1 0]
            Forbidden words: [[1, 1]]
            Order: 1

        #. Information on the standard Whistle Free shift::

            sage: X = sfts.WhistleFree([("01",2), ("10",2)]); X.get_info()
            Name: generalized Whistle-Free shift
            Alphabet: [0, 1]
            Format: edge
            A representing matrix:
            [1 1 0 0 0 0 0 0]
            [0 0 1 1 0 0 0 0]
            [0 0 0 0 1 0 0 0]
            [0 0 0 0 0 0 1 1]
            [1 1 0 0 0 0 0 0]
            [0 0 0 1 0 0 0 0]
            [0 0 0 0 1 1 0 0]
            [0 0 0 0 0 0 1 1]
            Forbidden words: [[0, 1, 0, 1], [1, 0, 1, 0]]
            Order: 3

        AUTHORS:

        - Sebastian Donoso (2011): initial version
        """
        if self._is_empty == False:
            print("""Name: %s
Alphabet: %s
Format: %s
A representing matrix:
%s
Forbidden words: %s
Order: %d""" % (self._name, self._alph, self._format, self._matrix,
                self._fwords, self.order()))
        else:
            print("The empty shift.")

    def reduce(self):
        r"""
        Returns an essential presentation of the SFT (removing all sources and
        sinks from the graph).

        As stranded vertices in a graph presentation can not be visited in any
        bi-infinite walk on the graph these vertices do never show up in
        points of the SFT and thus can (should) be removed without any loss of
        generality.

        OUTPUT:

        - an SFT that has an essential graph representation and is equivalent
          to the input SFT

        EXAMPLES:

        #. Reducing a non-essential graph presentation::

            sage: X = SFT([[0,1,0], [0,1,1], [1,1,0]]); X
            A 2-step SFT over [0, 1].
            sage: X.graph()
            Looped multi-digraph on 4 vertices
            sage: X.reduce().graph()
            Looped multi-digraph on 2 vertices

        #. Reducing a big graph presentation to a single 3-cycle::

            sage: G = DiGraph({1:[2,7], 2:[3,8], 3:[4,9], 4:[5,9,10,11], 5:[6,11], 6:[12], 7:[8,13], 8:[14], 9:[8,15], 10:[16], 11:[16,17], 12:[11], 13:[14], 14:[15], 15:[10], 17:[12]})
            sage: X = SFT(G, format='vertex', symb_sep='.'); X
            A 1-step SFT over [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17].
            sage: X.reduce()
            A 1-step SFT over [11, 12, 17].

        AUTHORS:

        - Sebastian Donoso (2011): initial version
        - Michael Schraudner (2012): simplified method, cleanup,
          mschraudner@dim.uchile.cl
        """
        if self._is_empty:
            return self
        X = copy(self)
        G = copy(self._graph)
        X._fwords = copy(self._fwords)
        if self._symb_sep == "":
            lenth = len(str(self._alph[0]))
            stranded = [v for v in G.vertices() # find sources and sinks
                          if G.in_degree(v) == 0 or G.out_degree(v) == 0]
            while not stranded == []:
                G.delete_vertices(stranded)
                for v in stranded:
                    symbs = [v[i:i + lenth] for i in range(0, len(v), lenth)]
                    X._fwords.append([b for a in symbs for b in self._alph
                                                             if str(b) == a])
                stranded = [v for v in G.vertices()
                              if G.in_degree(v) == 0 or G.out_degree(v) == 0]
            symbs = []
            for s in G.edge_labels():
                symbs.extend([s[i:i + lenth] for i in range(0, len(s), lenth)])
        else:
            stranded = [v for v in G.vertices() # find sources and sinks
                          if G.in_degree(v) == 0 or G.out_degree(v) == 0]
            while not stranded == []:
                G.delete_vertices(stranded)
                for v in stranded:
                    symbs = v.split(self._symb_sep)
                    X._fwords.append([b for a in symbs for b in self._alph
                                                             if str(b) == a])
                stranded = [v for v in G.vertices()
                              if G.in_degree(v) == 0 or G.out_degree(v) == 0]
            symbs = []
            for s in G.edge_labels():
                symbs.extend(s.split(self._symb_sep))
        X._alph = filter(lambda a: str(a) in symbs, self._alph)
        X._fwords = [w for w in X._fwords
                       if [a for a in w if not a in X._alph] == []]
        X._fwords, X._alph = X.first_offenders('list')
        X._create_fwords_DiGraph()
        return X

    def first_offenders(self, format='string'):
        r"""
        Returns the pair ``(first_offend, alph_red)`` where ``first_offend``
        is the list of first offender words and ``alph_red`` is the reduced
        alphabet of the SFT corresponding to these first offenders.

        For a given SFT a word is said to belong to the set of first offenders
        if it does not belong to the language of the SFT, but every proper
        subword of it does.

        In particular, the list of first offenders is a (the) minimal set of
        forbidden words for the given SFT in the sense that it defines the
        same shift space as the original list of forbidden words, but in
        addition it has the least possible length (i.e. sum of the lengths of
        its elements) under all such lists of forbidden words.

        INPUT:

        - ``format`` - (default: 'string') a string defining the format in
          which the set of first offenders is output; it can be one of the
          following:

            - 'string': the first offenders are returned as a list of strings.
              If the alphabet has symbols of distinct lengths and a symbol
              separator was given when defining the SFT, symbols in the first
              offender words/strings will be separated by this separator
            - 'list': the first offenders are returned as a list of lists
              containing alphabet symbols
            - 'Word': the first offenders are returned as a list of SAGE Words

          If ``format`` is not specified, the format 'string' is chosen by
          default.

        OUTPUT:

        - a tuple ``(first_offend, alph_red)`` where:

          - ``first_offend`` is the list of first offenders whose format
            depends on the parameter ``format`` and
          - ``alph_red`` is the list of symbols in the reduced alphabet of the
            SFT corresponding to the list of first offenders

        REFERENCES:

        - 'D. Lind, B. Marcus: An introduction to symbolic dynamics and coding. Sections 1.3 and 2.1'

        EXAMPLES:

        #. Getting the list of First Offenders of some SFTs::

                sage: S = SFT(["0000", "11", "010101", "101010"], ["0", "1"])
                sage: S.first_offenders()
                (['11', '0000', '10101'], ['0', '1'])
                sage: S.first_offenders(format='Word')
                ([word: 11, word: 0000, word: 10101], ['0', '1'])

           ::

                sage: X = SFT(["000", "10001", "101", "110", "11111", "01100"], ["0", "1"])
                sage: X.first_offenders('list')
                ([['1', '1'], ['0', '0', '0'], ['1', '0', '1']], ['0', '1'])

           ::

                sage: Y = sfts.RLL(2,6,3)
                sage: Y.forbidden_words()
                ['0000000', '11', '101', '10001', '100001', '10000001']
                sage: Y.first_offenders()
                (['11', '101', '10001', '000000', '100001'], [0, 1])

        #. A non-obvious way to define a full-shift on a single symbol::

            sage: X = SFT(["10", "11"], ["0", "1"])
            sage: X.forbidden_words()
            ['10', '11']
            sage: X.is_fullshift()
            True
            sage: X.first_offenders()
            ([], ['0'])

        ALGORITHM:

        The following is a rough description of the algorithm used to find the
        set of first offenders and the reduced alphabet:

        #. If there are no forbidden words, then the shift is a full-shift
           over the original alphabet.

        #. If there are only symbols as forbidden words, then the shift is the
           full-shift over the reduced alphabet being the original alphabet
           minus the forbidden symbols.

        #. Remove all forbidden symbols from the alphabet.

        #. Remove the forbidden symbols and all words containing them from the
           list of forbidden words.

        #. Take max_len to be the size of the largest forbidden word and
           inductively create a sequence of max_len-1 lists where:

            - the first list is the reduced alphabet
            - the list m+1 is created concatenating the words of list m with
              the symbols of the reduced alphabet including only those words
              that do not contain any forbidden word

           The last list (called A) thus contains exactly the admissible words
           of length max_len-1 over the reduced alphabet.

        #. Create a list T of triples (w, L, R) where:

            - w is an element of the list A created in point 5
            - L is the list of elements of the list A created in point 5
              which overlap from the left with w and such that their
              overlapping concatenation is not a forbidden word
            - R is the same as L but with right overlapping

           (We say that u overlaps with w from the left if u[1:]==w[:len(u)-1]
           and its overlapping concatenation is u[0]+w. Similarly for overlaps
           from the right.)

        #. If for a triple (w, L, R), L or R are empty, remove this triple
           from the list T and remove the element w from all the lists L and R
           of the other triples. Repeat this until there are no triples with
           empty lists L or R left over.

        #. For l in {1,...,max_len}:
            Produce all subwords of length l (over the reduced alphabet)
            which appear in one of the forbidden words and store all of them
            (ordered by increasing length) as possible first offenders in a
            list pos_first_off.

        #. For every possible first offender u, see if it appears in the
           remaining triples (i.e if there exists some triple (w, L, R) with u
           a subword of w). If it does, then it can be continued to the left
           and right and thus it can not be a first offender. If u does not
           show up as subword of w in any of the triples, store u as a first
           offender and remove from the list of possible first offenders all
           the words that contain u.

        AUTHORS:

        - Michael Schraudner (2011): algorithm design, mschraudner@dim.uchile.cl
        - Sebastian Donoso (2011): initial implementation and documentation
        - Sebastian Barbieri (2012): adapted algorithm to support different
          forbidden words formats
        - Michael Schraudner (2012): final cleanup to include method into SAGE
        """
        if not format in ['string', 'list', 'Word']:
            raise ValueError("format must be either 'string' or 'list' or " +
                             "'Word'")
        ## point 1 ##
        if self._fwords == []:
            return [], self._alph
        else:
            lenth = [len(w) for w in self._fwords]
            max_len = max(lenth)
            ## point 2 ##
            if max_len == 1:
                return [], [a for a in self._alph
                              if not [a] in self._fwords]
            elif len(self._fwords) == 1: # only 1 forbidden word, not a symbol
                return self._fwords, self._alph # has to be a first offender
            else:
                ## points 3 and 4 ##
                min_len = min(lenth)
                if min_len == 1:
                    alph_red = [a for a in self._alph
                                  if not [a] in self._fwords]
                    fwords = [w for w in self._fwords
                                if [a for a in w if not a in alph_red] == []]
                else:
                    alph_red = self._alph
                    fwords = self._fwords
                ## point 5 ##
                A = [[a] for a in alph_red]
                for m in range(2, max_len):
                    A = [u + [a] for u in A for a in alph_red
                                 if [w for w in fwords
                                       if w == (u + [a])[-len(w):]] == []]
                ## point 6 ##
                T = []
                max_fwords = [w for w in fwords if len(w) == max_len]
                for w in A:
                    L = [u for u in A
                         if u[1:] == w[:-1] and not [u[0]] + w in max_fwords]
                    R = [u for u in A
                         if w[1:] == u[:-1] and not [w[0]] + u in max_fwords]
                    T.append((w, L, R))
                ## point 7 ##
                while True:
                    for t in T:
                        exist_empty_triple = False
                        if (t[1] == [] or t[2] == []):
                            exist_empty_triple = True
                            T.remove(t)
                            for s in T:
                                if t[0] in s[1]:
                                    s[1].remove(t[0])
                                if t[0] in s[2]:
                                    s[2].remove(t[0])
                    if not exist_empty_triple:
                        break
                ## point 8 ##
                pos_first_off = []
                for m in range(1, max_len + 1):
                    for w in fwords:
                        for k in range(len(w) - m + 1):
                            if not w[k:k + m] in pos_first_off:
                                pos_first_off.append(w[k:k + m])
                ## point 9 ##
                first_off = []
                while not pos_first_off == []:
                    w = pos_first_off.pop(0)
                    len_w = len(w)
                    if [t for t in T for j in range(len(t[0]) - len_w + 1)
                          if w == t[0][j:j + len_w]] == []:
                        first_off.append(w)
                        pos_first_off = filter(lambda u:
                                         [u for j in range(len(u) - len_w + 1)
                                            if w == u[j:j + len_w]] == [],
                                               pos_first_off)
                symbol_off = [w for w in first_off if len(w) == 1]
                if not symbol_off == []:
                    alph_red = [a for a in alph_red if not [a] in symbol_off]
                    first_off = [w for w in first_off if not w in symbol_off]
                ## output ##
                if format == 'list':
                    return first_off, alph_red
                elif format == 'Word':
                    return [Word(w) for w in first_off], alph_red
                elif format == 'string':
                    return [self._symb_sep.join([str(a) for a in w])
                                                 for w in first_off], alph_red

    def admissible_words_iterator(self, n, m=None, format='string'):
        r"""
        Returns an iterator generating the admissible words of the SFT of
        lengths between ``n`` and ``m`` sorted in length-lexicographic order.

        If the parameter ``m`` is not given, the method returns an iterator
        over the admissible words of length ``n`` only. If the parameter ``m``
        is ``+Infinity`` the method will return an iterator over all
        admissible words of length ``n`` and larger.

        INPUT:

        - ``n`` - a positive integer

        - ``m`` - (default: None) a positive integer (larger or equal than
          ``n``) or ``+Infinity``

        - ``format`` - (default: 'string') a string defining the format in
          which the admissible words are yielded; it can be one of the
          following:

            - 'string': the admissible words are yielded as strings. If the
              alphabet has symbols of distinct lengths and a symbol separator
              was given when defining the SFT, symbols in the admissible words
              will be separated by this separator
            - 'list': the admissible words are yielded as lists containing
              alphabet symbols
            - 'Word': the admissible words are yielded as SAGE Words

          If ``format`` is not specified, the format 'string' is chosen by
          default.

        OUTPUT:

        - an iterator yielding the admissible words of the SFT sorted in
          length-lexicographic order

        EXAMPLES:

        #. The first six admissible words of length 4 in the Fibonacci shift::

            sage: F = sfts.Fibonacci()
            sage: I = F.admissible_words_iterator(4, format='list')
            sage: for i in range(6): print(next(I))
            [0, 0, 0, 0]
            [0, 0, 0, 1]
            [0, 0, 1, 0]
            [0, 1, 0, 0]
            [0, 1, 0, 1]
            [1, 0, 0, 0]

        #. The admissible words of length between 1 and 5 not containing the
           subword '11' in a given SFT and the first six admissible words of
           length 60::

            sage: X = SFT([Word([0,1,0,1]), Word([1,0,0])], [0, 1])
            sage: [i for i in X.admissible_words_iterator(1, 5, 'Word') if not Word([1,1]).is_factor(i)]
            [word: 0, word: 1, word: 00, word: 01, word: 10, word: 000, word: 001, word: 101, word: 0000, word: 0001, word: 00000, word: 00001]
            sage: I = X.admissible_words_iterator(60)
            sage: for i in range(6): print(next(I))
            000000000000000000000000000000000000000000000000000000000000
            000000000000000000000000000000000000000000000000000000000001
            000000000000000000000000000000000000000000000000000000000011
            000000000000000000000000000000000000000000000000000000000110
            000000000000000000000000000000000000000000000000000000000111
            000000000000000000000000000000000000000000000000000000001101

        #. An iterator over all admissible words of a given SFT and the output
           of its first 20 elements::

            sage: X = SFT(matrix(3, 3, [0,1,0, 1,0,1, 0,1,1]), ['a', 'b', 'c', 'd', 'e'])
            sage: I = X.admissible_words_iterator(1, +Infinity)
            sage: w = " ".join(next(I) for t in range(20)); print(w)
            a b c d e ab ac ba cd ce db dc ed ee aba acd ace bab bac cdb

        #. An iterator for a run length limited SFT and its 1000th admissible
           word::

            sage: X = sfts.RLL(2,5)
            sage: I = X.admissible_words_iterator(1, +Infinity)
            sage: for i in range(1000): w = next(I)
            sage: print(next(I))
            1000010010010001

        AUTHORS:

        - Michael Schraudner (2012): initial version, mschraudner@dim.uchile.cl
        """
        def words_recursion(word, word_len, word_last):
            word_len += 1
            if word_len < current_len:
                for c in continuations[word_last]:
                    for w in words_recursion(word + [c[0]], word_len, c[1]):
                        yield w
            else:
                for c in continuations[word_last]:
                    yield word + [c[0]]

        if not (n in NN and n > 0):
            raise ValueError("n has to be a positive integer")
        if not (m is None or (m in NN and m >= n) or m == +Infinity):
            raise ValueError("m has to be a positive integer not less than " +
                             "n or +Infinity")
        if format == 'list':
            form = lambda w: w
        elif format == 'Word':
            form = lambda w: Word(w)
        elif format == 'string':
            form = lambda w: (self._symb_sep.join([str(a) for a in w]))
        else:
            raise ValueError("format has to be 'string', 'list' or 'Word'")
        if self._is_empty:
            return
        if m is None:
            m = n
        first_off, alph_red = self.first_offenders(format='list')
        if first_off:
            max_len = max([len(w) for w in first_off])
        else:
            max_len = 0
        admiss_words = [[]]
        for i in range(1, min(m, max_len) + 1):
            prefixes = admiss_words
            admiss_words = []
            for p in prefixes:
                for a in alph_red:
                    w = p + [a]
                    if all([u != w[-len(u):] for u in first_off]):
                        admiss_words.append(w)
            if n <= i and i <= m:
                for w in admiss_words:
                    yield form(w)
        if m > max_len:
            num_prefixes = len(prefixes)
            continuations = [0]*num_prefixes
            for i in range(num_prefixes):
                continuations[i] = [(w[-1], prefixes.index(w[1:]))
                                    for w in admiss_words
                                    if w[:-1] == prefixes[i]]
            current_len = max(n-1, max_len)
            while current_len < m:
                current_len += 1
                for i in range(num_prefixes):
                    for t in words_recursion(prefixes[i], max_len-1, i):
                        yield form(t)

    def admissible_words(self, n, m=None, format='string'):
        r"""
        Returns the list of all admissible words of the SFT of lengths between
        ``n`` and ``m`` sorted in length-lexicographic order.

        If the parameter ``m`` is not given, the method returns the list of
        admissible words of length ``n`` only.

        INPUT:

        - ``n`` - a positive integer

        - ``m`` - (default: None) a positive integer (larger or equal than
          ``n``)

        - ``format`` - (default: 'string') a string defining the format of
          the admissible words; it can be one of the following:

            - 'string': the admissible words are returned as a list of
              strings. If the alphabet has symbols of distinct lengths and a
              symbol separator was given when defining the SFT, symbols in the
              admissible words will be separated by this separator
            - 'list': the admissible words are returned as a list of lists
              containing alphabet symbols
            - 'Word': the admissible words are returned as a list of SAGE
              Words

          If ``format`` is not specified, the format 'string' is chosen by
          default.

        OUTPUT:

        - a list containing all of the SFT's admissible words of length equal
          to ``n``, or in the range ``n`` up to ``m`` (inclusive) sorted in
          length-lexicographic order

        EXAMPLES:

        #. Admissible words of different sizes (and different formats) for a
           given SFT::

            sage: X = SFT(["0101", "100", "111"], ["0", "1"])
            sage: X.admissible_words(4)
            ['0000', '0001', '0011', '0110', '1011', '1101']
            sage: X.admissible_words(1, 3, 'list')
            [['0'], ['1'], ['0', '0'], ['0', '1'], ['1', '0'], ['1', '1'], ['0', '0', '0'], ['0', '0', '1'], ['0', '1', '1'], ['1', '0', '1'], ['1', '1', '0']]
            sage: X.admissible_words(2, format='Word')
            [word: 00, word: 01, word: 10, word: 11]

        #. Admissible words of size between 2 and 4 for the Fibonacci shift
           and the number of admissible words of fixed lengths from 1 up to
           15::

            sage: F = sfts.Fibonacci(); F
            The Fibonacci shift of order 1 over [0, 1].
            sage: F.admissible_words(2, 4)
            ['00', '01', '10', '000', '001', '010', '100', '101', '0000', '0001', '0010', '0100', '0101', '1000', '1001', '1010']
            sage: [len(F.admissible_words(n)) for n in range(1, 16)]
            [2, 3, 5, 8, 13, 21, 34, 55, 89, 144, 233, 377, 610, 987, 1597]

        #. There are no admissible words in the empty shift::

            sage: X = SFT([[0, 0], [0, 1], [1, 1, 1, 1]], [0,1])
            sage: X.admissible_words(10)
            []

        #. If ``n`` or ``m`` are not positive integers satisfying ``n <= m``,
           an error will be raised::

            sage: X = sfts.Fibonacci()
            sage: X.admissible_words(-3)
            Traceback (most recent call last):
            ...
            ValueError: n has to be a positive integer
            sage: X.admissible_words(2, 0)
            Traceback (most recent call last):
            ...
            ValueError: m has to be a positive integer not less than n

        AUTHORS:

        - Sebastian Donoso (2011): initial version
        - Sebastian Barbieri (2012): added documentation, changed algorithm
          and fixed various problems
        - Michael Schraudner (2012): added parameter checks, variable cleanup,
          mschraudner@dim.uchile.cl
        """
        if not (n in NN and n > 0):
            raise ValueError("n has to be a positive integer")
        if not (m is None or (m in NN and m >= n)):
            raise ValueError("m has to be a positive integer not less than n")
        if self._is_empty:
            return []
        if m is None:
            m = n
        first_off, alph_red = self.first_offenders(format='list')
        words = [[]]
        for i in range(1,n):
            prefixes = words
            words = []
            for p in prefixes:
                for a in alph_red:
                    w = p + [a]
                    if all([u != w[-len(u):] for u in first_off]):
                        words.append(w)
        admissible_words = []
        for i in range(n,m+1):
            prefixes = words
            words = []
            for p in prefixes:
                for a in alph_red:
                    w = p + [a]
                    if all([u != w[-len(u):] for u in first_off]):
                        words.append(w)
            admissible_words.extend(words)
        if format == 'list':
            return admissible_words
        elif format == 'Word':
            return [Word(w) for w in admissible_words]
        elif format == 'string':
            return [self._symb_sep.join([str(a) for a in w])
                                                    for w in admissible_words]
        else:
            raise ValueError("format has to be 'string', 'list' or 'Word'")

    def num_admissible_words(self, n):
        r"""
        Returns the number of admissible words of length  ``n`` in the SFT.

        INPUT:

        - ``n`` - a positive integer

        OUTPUT:

        - a non-negative integer, the number of admissible words of length
          ``n``

        EXAMPLES:

        #. The number of admissible words in a full shift on three symbols::

            sage: X = sfts.Full(3)
            sage: for n in range(1,5): X.num_admissible_words(n)
            3
            9
            27
            81

        #. The cardinalities of admissible words in the Fibonacci shift are
           Fibonacci numbers::

            sage: F = sfts.Fibonacci()
            sage: for n in range(1,10): F.num_admissible_words(n)
            2
            3
            5
            8
            13
            21
            34
            55
            89

        AUTHORS:

        - Michael Schraudner (2012): initial version, mschraudner@dim.uchile.cl
        """
        if not (n in NN and n > 0):
            raise ValueError("n has to be a positive integer")
        if self._is_empty:
            return 0
        first_off, alph_red = self.first_offenders(format='list')
        words = [[]]
        for i in range(1,n+1):
            prefixes = words
            words = []
            for p in prefixes:
                for a in alph_red:
                    w = p + [a]
                    if all([u != w[-len(u):] for u in first_off]):
                        words.append(w)
        return len(words)


    def num_periodic_points(self, n, least=False):
        r"""
        Returns the number of periodic points of (least) period ``n`` in the
        SFT.

        If the parameter ``least`` is set to False (default), the method
        returns the number of (all) periodic point of period ``n``, if it is
        set to True the method returns the number of periodic points of least
        period ``n`` only.

        INPUT:

        - ``n`` - a positive integer

        - ``least`` - (default: False) a boolean

        OUTPUT:

        - a non-negative integer, the number of periodic points of (least)
          period ``n``

        EXAMPLES:

        #. The number of periodic points in a full shift on three symbols::

            sage: X = sfts.Full(3)
            sage: for n in range(1,5): X.num_periodic_points(n)
            3
            9
            27
            81
            sage: for n in range(1,5): X.num_periodic_points(n, least=True)
            3
            6
            24
            72

        #. The number of orbits of length 1 to 9 in the Fibonacci shift::

            sage: F = sfts.Fibonacci()
            sage: for n in range(1,10): F.num_periodic_points(n, least=True)/n
            1
            1
            1
            1
            2
            2
            4
            5
            8

        AUTHORS:

        - Michael Schraudner (2012): initial version, mschraudner@dim.uchile.cl
        """
        if not (n in NN and n > 0):
            raise ValueError("n has to be a positive integer")
        if self._is_empty:
            return 0
        if not least:
            return((self._matrix**n).trace())
        else:
            t = 0
            for d in divisors(n):
                t += moebius(n/d)*(self._matrix**d).trace()
            return(t)


    def random_element(self):
        r"""
        Returns a random element of the SFT.

        The element, being a (bi-)infinite sequence of symbols, is given as an
        iterator over the symbols in a random point of the SFT. The iterator
        thus actually yields consecutive symbols of such a point/sequence.

        OUTPUT:

        - yields an iterator providing consecutive symbols of a random point
          in the SFT

        EXAMPLES:

        #. A (random) point in the Fibonacci shift never contains the word
           "11"::

            sage: F = sfts.Fibonacci()
            sage: I = F.random_element()
            sage: w = "".join([str(next(I)) for i in range(10000)])
            sage: "11" in w
            False

        #. A (random) word of length 10 in the (2,7)-RLL shift for sure
           contains the subword "010" but does not contain the subword "101"::

            sage: X = sfts.RLL(2,7)
            sage: I = X.random_element()
            sage: w = "".join([str(next(I)) for i in range(10)])
            sage: "010" in w
            True
            sage: "101" in w
            False

        AUTHORS:

        - Michael Schraudner (2012): initial version, mschraudner@dim.uchile.cl
        """
        first_off, alph_red = self.first_offenders(format='list')
        if first_off:
            max_len = max([len(w) for w in first_off])
        else:
            max_len = 0
        admiss_words = [[a] for a in alph_red]
        for i in range(2, max_len + 1):
            prefixes = admiss_words
            admiss_words = []
            for p in prefixes:
                for a in alph_red:
                    w = p + [a]
                    if all([u != w[-len(u):] for u in first_off]):
                        admiss_words.append(w)
        num_prefixes = len(prefixes)
        continuations = [0]*num_prefixes
        for i in range(num_prefixes):
            continuations[i] = [(w[-1], prefixes.index(w[1:]))
                                for w in admiss_words
                                if w[:-1] == prefixes[i]]
        last = randint(0, num_prefixes-1)
        while True:
            c = continuations[last]
            i = randint(0, len(c)-1)
            last = c[i][1]
            yield c[i][0]

    def order(self):
        r"""
        Returns the order of the SFT.

        This method returns the length of the longest first offender
        (forbidden word). Hence the order of a full shift is 0 (no forbidden
        words at all), while all other SFTs have a positive order.

        OUTPUT:

        - a non-negative integer representing the order of the SFT

        EXAMPLES:

        #. Getting the order of some SFTs::

            sage: SFT().order()
            0
            sage: F = sfts.Fibonacci()
            sage: F.order()
            1
            sage: SFT(["011"], ["0", "1"]).order()
            2
            sage: sfts.RLL(2,6,3).order()
            5
            sage: X = SFT([[0,1,0], [0,1,1], [1,1,0]])
            sage: X.order()
            2
            sage: X.reduce().order()
            1

        AUTHORS:

        - Sebastian Donoso (2011): initial version
        """
        if self._is_empty:
            return 0
        first_off = self.first_offenders("list")[0]
        if first_off == []:
            return 0
        else:
            return max([len(w) for w in first_off])-1

    def period(self):
        r"""
        Returns the period of an SFT.

        This invariant is defined as the greatest common divisor of the
        periods of all periodic points in the SFT which is the same as the
        greatest common divisor of the lengths of all cycles in a (any) graph
        presentation of the SFT.

        If the SFT is not irreducible the period of every (non-trivial)
        irreducible component is computed and the greatest common divisor of
        those periods is returned.

        The period of the empty shift is returned as +Infinity.

        OUTPUT:

        - a non-negative integer representing the period of the SFT

        EXAMPLES:

        #. Period of some irreducible SFTs::

            sage: SFT([], ['a', 'b', 'c']).period()
            1
            sage: M = matrix(6, 6, [0,0,0,0,0,1, 0,0,0,0,1,0, 1,0,0,1,0,0, 0,1,0,0,0,0, 0,1,1,0,0,0, 1,0,0,0,1,0])
            sage: X = SFT(M)
            sage: X.period()
            2
            sage: M = matrix(5, 5, [0,2,0,0,0, 0,0,2,0,0, 2,0,0,1,0, 0,0,0,0,3, 0,0,1,0,0])
            sage: X = SFT(M)
            sage: X.period()
            3

        #. Period of some reducible or even not essential SFTs::

            sage: X = SFT(matrix(4, 4, [4,0,3,1, 1,2,1,0, 0,0,2,0, 1,0,1,0]))
            sage: X.period()
            1
            sage: X = SFT(["0000","0011","0101","0110","1001","1010","1111"], ["0","1"])
            sage: X.period()
            4
            sage: X = SFT(matrix(4, 4, [0,1,2,1, 0,0,2,1, 0,3,0,2, 0,0,0,0]))
            sage: X.period()
            2

        #. Period of the empty shift::

            sage: SFT().period()
            +Infinity

        AUTHORS:

        - Michael Schraudner (2012): initial version, mschraudner@dim.uchile.cl
        """
        if self._is_empty:
            return +Infinity
        irr_comps = self._graph.strongly_connected_components()
        if len(irr_comps) == 1:
            L = [set([i]) for i in range(self._matrix.ncols())]
            for k in range(self._matrix.ncols()):
            # each iteration reduces the size of L, until it has final size
                for l in L:
                    S = set([])
                    for i in l:
                        for j in range(self._matrix.ncols()):
                            if self._matrix[i, j] > 0:
                                for m in L:
                                    if j in m:
                                        L.remove(m)
                                        S.update(m)
                                        break
                    L.append(S)
            return len(L) # the remaining sets in L give the period
        else:
            vertices = self._graph.vertices()
            per = 0
            for comp in irr_comps:
                comp_ind = [i for c in comp
                             for i in range(self._graph.num_verts())
                             if c == vertices[i]]
                L = [set([i]) for i in comp_ind]
                for k in comp_ind:
                    for l in L:
                        S = set([])
                        for i in l:
                            for j in comp_ind:
                                if self._matrix[i, j] > 0:
                                    for m in L:
                                        if j in m:
                                            L.remove(m)
                                            S.update(m)
                                            break
                        if S: # if S was updated (i.e. is not empty)
                            L.append(S)
                        else:
                            L.remove(l)
                per = gcd(per, len(L))
            return per

    def is_fullshift(self):
        r"""
        Returns ``True`` if the SFT is a full-shift, and ``False`` otherwise.

        OUTPUT:

        - ``bool``: True or False

        EXAMPLES:

        #. A full-shift on 4 symbols and the Fibonacci shift::

            sage: sfts.Full(4).is_fullshift()
            True
            sage: sfts.Fibonacci().is_fullshift()
            False

        #. Some inputs that don't look like full-shifts but in fact are::

            sage: X = SFT(matrix(2, 2, [1,1, 1,1]), format='vertex')
            sage: X.is_fullshift()
            True
            sage: SFT(matrix(2, 2, [2,2, 0,0])).is_fullshift()
            True
            sage: SFT(matrix(3, 3, [1,0,0, 1,0,0, 1,1,0])).is_fullshift()
            True

        #. Some SFTs that aren't full-shifts::

            sage: X = SFT(matrix(2, 2, [1,1, 1,1]), format='edge')
            sage: X.is_fullshift()
            False
            sage: SFT(matrix(2, 2, [1,2, 1,0])).is_fullshift()
            False

        AUTHORS:

        - Sebastian Donoso (2011): initial version
        """
        return bool(self.first_offenders("list")[0] == [])

    def is_irreducible(self):
        r"""
        Returns ``True`` if the SFT is irreducible, and ``False`` otherwise.

        OUTPUT:

        - ``bool``: True or False

        EXAMPLES:

        #. Checking whether some SFTs are irreducible::

            sage: SFT(matrix(2, 2, [1,2, 1,0])).is_irreducible()
            True
            sage: SFT(["00", "11"]).is_irreducible()
            True
            sage: SFT(matrix(2, 2, [1,0, 2,1])).is_irreducible()
            False
            sage: SFT(["01", "10"]).is_irreducible()
            False

        #. The empty shift is considered irreducible::

            sage: SFT().is_irreducible()
            True

        AUTHORS:

        - Sebastian Donoso (2011): initial version
        """
        return self._graph.is_strongly_connected() # possible to speed this up?

    def irreducible_components(self, format='graph'):
        r"""
        Returns the list of irreducible components of the SFT in two possible
        formats, either as a list of graphs, or as a list of adjacency
        matrices.

        INPUT:

        - ``format`` - (default: 'graph') a string which can be either:

            - 'graph':  a list of the subgraphs representing the irreducible
              components of the SFT will be returned
            - 'matrix': a list of the adjacency matrices representing the
              irreducible components of the SFT will be returned

          If ``format`` is not given, 'graph' format is chosen by default.

        OUTPUT:

        - a list of the irreducible components either as subgraphs ('graph'
          format) or as adjacency matrices ('matrix' format)

        EXAMPLES:

        #. The irreducible componentes of an SFT in the 2 available formats::

            sage: M = matrix(5, 5, [1,1,0,0,0, 1,0,0,1,0, 1,0,0,0,1, 0,0,0,1,1, 0,0,0,1,0])
            sage: X = SFT(M, format='vertex')
            sage: X.irreducible_components('graph')
            [Subgraph of (): Looped multi-digraph on 2 vertices, Subgraph of (): Looped multi-digraph on 2 vertices, Subgraph of (): Looped multi-digraph on 1 vertex]
            sage: X.irreducible_components('matrix')
            [
            [1 1]  [1 1]
            [1 0], [1 0], [0]
            ]

        AUTHORS:

        - Sebastian Donoso (2011): initial version
        """
        if format == 'graph':
            return self._graph.strongly_connected_components_subgraphs()
        elif format == 'matrix':
            return [G.adjacency_matrix() for G in
                      self._graph.strongly_connected_components_subgraphs()]
        else:
            raise ValueError("format not valid, has to be either 'graph' or" +
                             " 'matrix'")

    def is_mixing(self):
        r"""
        Returns ``True`` if the SFT is (topologically) mixing, and ``False``
        otherwise.

        OUTPUT:

        - ``bool``: True or False

        EXAMPLES:

        #. Checking whether some irreducible SFTs are mixing::

            sage: SFT(matrix(2, 2, [1,2, 1,0])).is_mixing()
            True
            sage: SFT(["00", "11"]).is_mixing()
            False
            sage: X = SFT(matrix(5, 5, [0,1,0,0,0, 0,0,1,0,0, 0,0,0,1,0, 0,0,0,0,1, 1,0,0,1,0]))
            sage: X.is_mixing()
            True

        #. Non-irreducible SFTs clearly can not be mixing::

            sage: SFT(matrix(2, 2, [1,0, 1,1])).is_mixing()
            False
            sage: SFT(["01", "10"]).is_mixing()
            False

        #. The empty shift is considered mixing::

            sage: SFT().is_mixing()
            True

        AUTHORS:

        - Sebastian Donoso (2011): initial version
        - Sebastian Barbieri (2012): simplified using self.period() method
        """
        if self._is_empty:
            return True
        return bool(self.is_irreducible() and self.period() == 1)

    def entropy(self, logarit=True):
        r"""
        Returns the topological entropy (or its exponential) of the SFT.

        The topological entropy of an SFT is actually the logarithm of the
        Perron value (i.e. the largest non-negative real eigenvalue) of the
        (a) matrix representing the SFT.

        If the SFT is the empty shift, the method will return -Infinity (or
        its exponential being 0).

        INPUT:

        - ``logarit`` - (default: True) option to select the output value:

            - True: the logarithm of the Perron value (the actual entropy)
            - False: the actual Perron value (exponential of the entropy)

        OUTPUT:

        - a real number (or -Infinity) representing the topological entropy of
          the SFT or its exponential value

        EXAMPLES:

        #. Entropy of a full-shift::

            sage: sfts.Full(3).entropy()
            log(3)

        #. Entropy of the Fibonacci shift::

            sage: F = sfts.Fibonacci(); F.entropy()
            0.481211825059603
            sage: F.entropy(logarit=False)
            1.618033988749895?

        #. Entropy of a given SFT::

            sage: SFT(matrix(3, 3, [1,1,0, 1,0,2, 2,1,2])).entropy()
            1.19476321728711

        #. Entropy of the empty shift::

            sage: SFT().entropy()
            -Infinity
            sage: SFT().entropy(logarit=False)
            0

        AUTHORS:

        - Sebastian Donoso (2011): initial version
        - Sebastian Barbieri (2012): documentation and the empty shift case
        """
        if self._is_empty:
            if logarit:
                return -Infinity
            else:
                return 0
        else:
            lamda = max([x for x in self._matrix.eigenvalues() if x in RR])
            if logarit:
                return log(lamda)
            else:
                return lamda

    def zeta_function(self):
        r"""
        Returns the Zeta function of a given SFT.

        The (Artin-Mazur) zeta function of an SFT combines the combinatorial
        data about its periodic points. It is defined as the exponential of
        the generating function of the sequence of cardinalities (divided by
        n) of periodic points of period length n. For an SFT this always gives
        a rational function of the form 1 divided by a polynomial with integer
        coefficients.

        OUTPUT:

        - a rational function, in fact the inverse of a polynomial over the
          integers

        REFERENCES:

        - 'D. Lind, B. Marcus: An introduction to symbolic dynamics and coding. Section 6.4'

        EXAMPLES:

        #. Zeta function of the full-4-shift::

            sage: X = sfts.Full(4); X
            The full-shift on 4 symbols of order 0 over [0, 1, 2, 3].
            sage: X.zeta_function()
            -1/(4*t - 1)

        #. Zeta function of the Fibonacci shift::

            sage: X = sfts.Fibonacci(); X
            The Fibonacci shift of order 1 over [0, 1].
            sage: X.zeta_function()
            -1/(t^2 + t - 1)

        #. Zeta function of the RLL-(2,4) shift::

            sage: X = sfts.RLL(2, 4); X
            The (2,4)-RLL shift of order 4 over [0, 1].
            sage: X.zeta_function()
            -1/(t^5 + t^4 + t^3 - 1)

        #. Zeta function of the empty shift::

            sage: SFT().zeta_function()
            1

        AUTHORS:

        - Sebastian Barbieri (2012): initial version
        """
        t = var("t")
        f = (1 - t*(self._matrix)).determinant()
        f = f.expand() # every zeta function gets the same polynomial repr.
        return 1/f

    def nonzero_spectrum(self):
        r"""
        Returns the spectrum away from zero of an SFT.

        The spectrum away from zero of a matrix is defined as the set of its
        eigenvalues (excluding zero) counted with multiplicities. For an SFT
        given by some matrix presentation this set is a conjugacy invariant.

        OUTPUT:

        - a sorted list representing the spectrum away from zero of the SFT

        EXAMPLES:

        #. Spectrum away from zero of the Fibonacci shift::

            sage: X = sfts.Fibonacci(); X; X.nonzero_spectrum()
            The Fibonacci shift of order 1 over [0, 1].
            [-0.618033988749895?, 1.618033988749895?]

        #. Spectrum away from zero of a given SFT::

            sage: X = SFT(matrix(4, 4, [4,0,3,1, 1,2,1,0, 0,0,2,0, 1,0,1,0]))
            sage: X.nonzero_spectrum()
            [-0.2360679774997897?, 2, 2, 4.236067977499789?]

        #. Spectrum away from zero of two conjugate shifts::

            sage: X = SFT([], ['a', 'b', 'c']); X
            The full-shift on 3 symbols of order 0 over ['a', 'b', 'c'].
            sage: Y = SFT(matrix(3, 3, [1,1,1, 1,1,1, 1,1,1])); Y
            A 1-step SFT over [0, 1, 2, 3, 4, 5, 6, 7, 8].
            sage: X.nonzero_spectrum(); Y.nonzero_spectrum()
            [3]
            [3]

        AUTHORS:

        - Sebastian Barbieri (2012): initial version
        """
        eig_vals = [x for x in self._matrix.eigenvalues() if not x == 0]
        eig_vals.sort()
        return eig_vals

    def bowen_franks_group(self, poly=[1, -1], format='smith'):
        r"""
        Returns the (generalized) Bowen-Franks group of the SFT in the
        specified format.

        The (ordinary) Bowen-Franks group of an SFT defined by an NxN matrix A
        is the finitely generated additive abelian group given as the quotient
        of ``\ZZ^N`` by its image under the matrix ``(Id-A)``.

        The (generalized) Bowen-Franks group of an SFT is defined similarly
        replacing the matrix ``(Id-A)`` by a more general matrix given as
        ``p(A)`` where ``p(t)`` is a polynomial over the integers satisfying
        ``|p(0)|=1`` which can be specified via the parameter ``poly``. In
        this setting the ordinary Bowen-Franks group is recovered for the
        polynomial ``p(t) = 1 - t`` (the default value of the parameter
        ``poly``).

        Both definitions yield conjugacy invariants for the given SFT.

        The method returns the (generalized) Bowen-Franks group either as a
        finitely generated additive abelian group (i.e. a sum of cyclic
        groups) or as a matrix in Smith form. The correspondence between the
        two presentations is given by the theorem about finitely generated
        modules over a PID (or the structure theorem for finitely generated
        additive abelian groups) which implies that the sizes of the cyclic
        groups appearing in the sum are given by the elementary divisors
        appearing in the Smith form of the matrix.

        INPUT:

        - ``poly`` - (default: [1, -1]) an optional polynomial used in the
          definition of the generalized Bowen-Franks groups

          This parameter accepts either a polynomial over the integers, a
          symbolic expression representing a polynomial or a list of integers
          corresponding to the polynomial's coefficients ordered by ascending
          degree, for example ``poly = [1, 2, -4, 8]`` will be equivalent to
          the polynomial ``p(t) = 1 + 2*t - 4*t^2 + 8*t^3``.

          If ``poly`` is not specified, the polynomial ``p(t) = 1 - t``
          (corresponding to the ordinary Bowen-Franks group) is chosen by
          default.

        - ``format`` - (default: 'smith') a string specifying the output
          format; it can be either of the following:

            - 'smith': the Bowen-Franks group will be returned as a matrix in
              Smith form, containing the elementary divisors (in ascending
              order) in its diagonal
            - 'group': the Bowen-Franks group will be returned as a direct sum
              of cyclic groups

          If ``format`` is not specified, the format 'smith' is chosen by
          default.

        OUTPUT:

        - a non-negative square matrix over the integers in Smith form
          representing the (generalized) Bowen-Franks group
        - a finitely generated additive abelian group (i.e. a sum of cyclic
          groups) representing the (generalized) Bowen-Franks group

        REFERENCES:

        - 'D. Lind, B. Marcus: An introduction to symbolic dynamics and coding. Section 7.4'

        EXAMPLES:

        #. The ordinary Bowen-Franks group of a given shift::

            sage: X = SFT(matrix(2, 2, [4,1, 1,0])); X
            A 1-step SFT over [0, 1, 2, 3, 4, 5].
            sage: X.bowen_franks_group()
            [1 0]
            [0 4]

        #. The Bowen-Franks group of the full-N-shift is isomorphic to the
           cyclic group of order (N-1)::

            sage: X = sfts.Full(3); X
            The full-shift on 3 symbols of order 0 over [0, 1, 2].
            sage: X.bowen_franks_group()
            [2]
            sage: X.bowen_franks_group(format='group')
            Additive abelian group isomorphic to Z/2

        #. Some (generalized) Bowen-Franks groups of two strong shift
           equivalent (i.e. conjugate) SFTs::

            sage: X = SFT(matrix(2, 2, [1,3, 2,1])); X
            A 1-step SFT over [0, 1, 2, 3, 4, 5, 6].
            sage: Y = SFT(matrix(2, 2, [1,6, 1,1])); Y
            A 1-step SFT over [0, 1, 2, 3, 4, 5, 6, 7, 8].
            sage: X.bowen_franks_group(); Y.bowen_franks_group()
            [1 0]
            [0 6]
            [1 0]
            [0 6]
            sage: l = [1,5,8]    # p(t) = 1 + 5*t + 8*t^2
            sage: X.bowen_franks_group(l); Y.bowen_franks_group(l)
            [   1    0]
            [   0 1198]
            [   1    0]
            [   0 1198]

        #. Some other ways to compute generalized Bowen-Franks groups::

            sage: X = SFT(matrix(3, 3, [1,3,4, 2,1,0, 0,1,3]))
            sage: expr = x^4 - 4*x^2 - 1    # p(t) = t^4 - 4*t^2 - 1
            sage: R.<t> = ZZ['t']; p = t^4 - 4*t^2 - 1
            sage: l = [-1, 0, -4, 0, 1]
            sage: expr; p; l
            x^4 - 4*x^2 - 1
            t^4 - 4*t^2 - 1
            [-1, 0, -4, 0, 1]
            sage: X.bowen_franks_group(expr)
            [  4   0   0]
            [  0   4   0]
            [  0   0 356]
            sage: X.bowen_franks_group(p)
            [  4   0   0]
            [  0   4   0]
            [  0   0 356]
            sage: X.bowen_franks_group(l)
            [  4   0   0]
            [  0   4   0]
            [  0   0 356]

        #. If ``format='group'``, an additive abelian group will be returned::

            sage: X = SFT(matrix(3, 3, [1,3,4, 2,1,0, 0,1,3]))
            sage: X.bowen_franks_group(format='smith')
            [1 0 0]
            [0 2 0]
            [0 0 2]
            sage: X.bowen_franks_group(format='group')
            Additive abelian group isomorphic to Z/2 + Z/2
            sage: X.bowen_franks_group([-1, 3, 4], 'smith')
            [  1   0   0]
            [  0 178   0]
            [  0   0   0]
            sage: X.bowen_franks_group([-1, 3, 4], 'group')
            Additive abelian group isomorphic to Z/178 + Z

        #. An error will be returned if the polynomial p(t) is not over ZZ or
           does not satisfy ``|p(0)|=1``::

                sage: X = SFT(matrix(3, 3, [1,3,4, 2,1,0, 0,1,3]))
                sage: X.bowen_franks_group([1, 1/3])    # p(t) = 1 + t/3
                Traceback (most recent call last):
                ...
                ValueError: poly must have integer coefficients

           ::

                sage: R.<t> = QQ['t']; p = t^4 - (1/4)*t^2 - 1
                sage: X.bowen_franks_group(p)
                Traceback (most recent call last):
                ...
                ValueError: poly must be a polynomial with base ring ZZ

           ::

                sage: X.bowen_franks_group([2, 4])    # p(t) = 2 + 4*t
                Traceback (most recent call last):
                ...
                ValueError: poly does not satisfy |p(0)|=1

        AUTHORS:

        - Sebastian Barbieri (2012): initial version
        - Michael Schraudner (2012): extended documentation, cleanup,
          mschraudner@dim.uchile.cl
        """
        if not format in ['smith', 'group']:
            raise ValueError("format must be either 'smith' or 'group'")
        if isinstance(poly, list):
            if poly == []:
                raise ValueError("poly does not satisfy |p(0)|=1")
            else:
                if not [c for c in poly if not c in ZZ] == []:
                    raise ValueError("poly must have integer coefficients")
                if not abs(poly[0]) == 1:
                    raise ValueError("poly does not satisfy |p(0)|=1")
                A = sum([poly[k] * (self._matrix**k)
                                                 for k in range(len(poly))])
        else:
            if isinstance(poly, Expression):
                try:
                    p = poly.polynomial(ZZ)
                except:
                    raise TypeError("poly must be a polynomial with base " +
                                    "ring ZZ")
            elif isinstance(poly, Polynomial):
                try:
                    p = poly.change_ring(ZZ)
                except:
                    raise ValueError("poly must be a polynomial with base " +
                                    "ring ZZ")
            else:
                raise TypeError("poly must be a coefficients list, an " +
                                "integer polynomial or a symbolic " +
                                "expression representing such a polynomial")
            if not abs(p.coefficients()[0]) == 1:
                raise ValueError("poly does not satisfy |p(0)|=1")
            A = p(self._matrix)
        A = A.smith_form()[0]
        if format == 'smith':
            return A
        else:
            return AdditiveAbelianGroup(A.diagonal())

    def parry_measure(self):
        r"""
        Returns a tuple ``(p, P)`` representing the Parry measure of an
        irreducible SFT presented as a vertex shift.

        Here ``p`` is the stationary probability vector (initial distribution)
        and ``P`` is the stochastic transition matrix of a particular Markov
        measure called the Parry measure of the irreducible SFT given as a
        vertex shift. This measure is the unique shift-invariant probability
        measure of maximal (measure-theoretic) entropy and it also satisfies
        the Gibbs condition (making it a Gibbs measure).

        OUTPUT:

        - a vector ``p`` (the stationary probability vector)
        - a matrix ``P`` (the stochastic transition matrix)

        REFERENCES:

        - 'B. Kitchens: Symbolic dynamics. One-sided, two-sided and countable state Markov shifts.'
        - 'M. Denker, C. Grillenberger, K. Sigmund: Ergodic theory on compact spaces.'

        EXAMPLES:

        #. Parry measure of an irreducible SFT in vertex representation::

            sage: X = SFT(matrix(3, 3, [1,0,1, 1,1,1, 0,1,0]), format='vertex')
            sage: X.parry_measure()
            (
                           [1/2   0 1/2]
                           [1/4 1/2 1/4]
            [1/4 1/2 1/4], [  0   1   0]
            )

        #. The Parry measure of the Fibonacci shift::

            sage: X = sfts.Fibonacci()
            sage: X.parry_measure()
            (
            [0.7236067977499790? 0.2763932022500211?],
            [ 0.618033988749895? 0.3819660112501051?]
            [                  1                   0]
            )

        #. Parry measure of a full-shift (given as a vertex shift)::

                sage: X = SFT(matrix(2, 2, [1,1, 1,1]), format='vertex')
                sage: X.parry_measure()
                (
                           [1/2 1/2]
                [1/2 1/2], [1/2 1/2]
                )

           ::

                sage: X = sfts.Full(2); X
                The full-shift on 2 symbols of order 0 over [0, 1].
                sage: X.parry_measure()
                doctest:...: UserWarning: the SFT has to be a vertex shift; will convert it
                (
                           [1/2 1/2]
                [1/2 1/2], [1/2 1/2]
                )

        #. A non-irreducible SFT will produce a NonImplementedError with the
           option to compute the Parry measure for each irreducible component::

                sage: X = SFT(matrix(2, 2, [1,1, 0,1]), format='vertex')
                sage: X.parry_measure()
                Traceback (most recent call last):
                ...
                NotImplementedError: computing the Parry measure requires an irreducible
                SFT; try computing it for each irreducible component separately

        AUTHORS:

        - Sebastian Barbieri (2012): initial version
        """
        if self.is_irreducible() == False:
            raise NotImplementedError("computing the Parry measure requires" +
                                      " an irreducible SFT; try computing " +
                                      "it for each irreducible component " +
                                      "separately")
        if self._format == 'vertex':
            A = self._matrix
        else:
            warnings.warn("the SFT has to be a vertex shift; will convert it")
            # instead to create a new SFT, could we just use a higher graph?
            A = SFT(self._fwords, self._alph, format='vertex')._matrix
        vecs_left = A.left_eigenvectors()
        vecs_right = A.right_eigenvectors()
        Perron = max([x for x in A.eigenvalues() if x in RR])
        for u in vecs_left:
            if u[0] == Perron:
                U = u[1][0]
                break
        for v in vecs_right:
            if v[0] == Perron:
                V = v[1][0]
                break
        U = U / (U.inner_product(V))
        n = A.ncols()
        p = matrix(1, n, map(lambda x,y: x*y, U,V))
        P = matrix(n, n, [A[i,j] * V[j] / V[i] / Perron for i in range(n)
                                                        for j in range(n)])
        return p, P

    def intersection(self, other, name=None):
        r"""
        Returns the intersection of two SFTs ``self`` and ``other``.

        The intersection of two SFTs is another SFT over an alphabet which is
        (a subset of) the intersection of the two SFT's alphabets and whose
        set of forbidden words equals the union of the two SFT's finite lists
        of forbidden words.

        The name of the intersection SFT can be specified using the parameter
        ``name``. Note also that the intersection will have the format
        ('vertex' or 'edge') and the symbol separator of the SFT ``self``
        while the format and separator of ``other`` is discarded.

        INPUT:

        - ``other`` - a Subshift Of Finite Type

        - ``name`` - (default: None) a string specifying the name of the new
          SFT

        OUTPUT:

        - a Subshift of Finite Type being the intersection of ``self`` with
          ``other``

        EXAMPLES:

        #. Two examples of non-trivial intersections of two SFTs::

                sage: X = SFT(['0.0.0.0','1.1.0','0.1.1.1','1.0.1'], symb_sep='.')
                sage: Y = SFT(['1,1,1','1,1,0,1','0,1,1,1'], symb_sep=',')
                sage: Z = X.intersection(Y); Z
                A 3-step SFT over ['0', '1'].
                sage: Z.forbidden_words()
                ['0.0.0.0', '1.1.0', '0.1.1.1', '1.0.1', '1.1.1', '1.1.0.1']

           ::

                sage: X = SFT(matrix(3, 3, [1,1,0, 0,1,1, 0,1,0]), [0, 1, 2], format='vertex')
                sage: Y = SFT([[1, 1, 2], [2, 2], [2, 3, 2], [3, 1]])
                sage: Z = X.intersection(Y, 'SFT'); Z
                The SFT of order 2 over [1, 2].
                sage: Z.forbidden_words()
                ['22', '112']
                sage: Z.get_info()
                Name: SFT
                Alphabet: [1, 2]
                Format: vertex
                A representing matrix:
                [1 0 0]
                [0 0 1]
                [1 1 0]
                Forbidden words: [[2, 2], [1, 1, 2]]
                Order: 2

        #. An intersection resulting in a trivial shift consisting of a single
           point::

            sage: X = SFT(matrix(2, 2, [1,1, 0,1]), format='vertex')
            sage: F = sfts.Fibonacci()
            sage: X.intersection(F)
            A 0-step SFT over [0, 1].
            sage: X.intersection(F).reduce()
            The full-shift on 1 symbols of order 0 over [0].

        #. An intersection resulting in the empty shift::

            sage: X = SFT([], ['a', 'c'])
            sage: F = SFT(matrix(2, 2, [1,1, 1,0]))
            sage: Y = X.intersection(F); Y
            The empty shift.

        #. An intersection of a full-2-shift in edge presentation with a shift
           on 3 symbols::

            sage: F = SFT(matrix(2, 2, [1,1, 1,1]), ['a', 'b', 'c', 'd'])
            sage: X = SFT(['ab'], ['a', 'b', 'c'])
            sage: Y = X.intersection(F); Y
            A 1-step SFT over ['a', 'b', 'c'].
            sage: Y.forbidden_words()
            ['ab', 'ac', 'ba', 'bb', 'cc']
            sage: Y.admissible_words(1, 2)
            ['a', 'b', 'c', 'aa', 'bc', 'ca', 'cb']
            sage: X.admissible_words(1, 2)
            ['a', 'b', 'c', 'aa', 'ac', 'ba', 'bb', 'bc', 'ca', 'cb', 'cc']
            sage: F.admissible_words(1, 2)
            ['a', 'b', 'c', 'd', 'aa', 'ab', 'bc', 'bd', 'ca', 'cb', 'dc', 'dd']

        AUTHORS:

        - Sebastian Donoso (2011): initial version
        - Sebastian Barbieri (2012): documentation and minor improvements
        - Michael Schraudner (2012): included support for format, symb_sep and
          name, mayor cleanup, mschraudner@dim.uchile.cl
        """
        if not isinstance(other, SFT):
            raise TypeError("other must be an SFT")
        if not (name is None or isinstance(name, basestring)):
            raise TypeError("name must be a string or None")
        X = copy(self)
        X._name = name
        X._alph = filter(lambda a: a in other.alphabet(), self._alph)
        if X._alph == []:
            X._obj = []
            X._create_empty()
            return X
        X._fwords = copy(self._fwords)
        X._fwords.extend(filter(lambda w: not w in self._fwords,
                                            other.forbidden_words('list')))
        X._fwords = filter(lambda w: [a for a in w if not a in X._alph] == [],
                                                               X._fwords)
        X._obj = X._fwords
        X._create_fwords_DiGraph()
        return X

    def cartesian_product(self, other, name=None):
        r"""
        Returns the cartesian product of the two SFTs ``self`` and ``other``.

        The cartesian product of two SFTs is another SFT over an alphabet
        which is a subset of the cartesian product of the two SFT's alphabets
        and whose set of forbidden words equals the union of the two sets
        ``{(u,w)}`` where either ``u`` is a forbidden word of ``self`` and
        ``w`` is any word of the same length in ``other`` or where ``w`` is a
        forbidden word in ``other`` and ``u`` is any word of the same length
        in ``self``

        The name of the cartesian product SFT can be specified using the
        parameter ``name``. Note also that the product will have the format
        ('vertex' or 'edge') of the SFT ``self`` while the format of ``other``
        is discarded. Finally the method tries first to re-use the symbol
        separator ``symb_sep`` from ``self``. If this fails (or if ``self``
        does not have a symbol separator) it will try the symbol separator of
        ``other`` and if this also fails will raise an error.

        INPUT:

        - ``other`` - a Subshift Of Finite Type

        - ``name`` - (default: None) a string specifying the name of the new
          SFT

        OUTPUT:

        - a Subshift of Finite Type being the cartesian product of ``self``
          with ``other``

        EXAMPLES:

        #. Cartesian product of an empty shift with another shift::

            sage: X = SFT(); Y = SFT(matrix(2, 2, [1,1, 1,1])); X; Y
            The empty shift.
            A 1-step SFT over [0, 1, 2, 3].
            sage: Z = X.cartesian_product(Y); Z
            The empty shift.

        #. Note that the entropy of the cartesian product is the sum of the
           entropies of the two SFTs::

            sage: X = SFT(matrix(2, 2, [2,0, 1,1])); F = sfts.Fibonacci()
            sage: Y = X.cartesian_product(F, name='SFT'); Y
            The SFT of order 1 over [(0, 0), (0, 1), (1, 0), (1, 1), (2, 0), (2, 1), (3, 0), (3, 1)].
            sage: (F.entropy() + X.entropy()).n() == Y.entropy()
            True

        #. Another example, a full-shift with a 2-cycle::

            sage: X = SFT(matrix(2, 2, [0,1, 1,0])); F = sfts.Full(3)
            sage: Y = X.cartesian_product(F); Y
            A 1-step SFT over [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2)].
            sage: Y.entropy()
            log(3)
            sage: Y.matrix()
            [0 0 0 1 1 1]
            [0 0 0 1 1 1]
            [0 0 0 1 1 1]
            [1 1 1 0 0 0]
            [1 1 1 0 0 0]
            [1 1 1 0 0 0]

        #. Careful with using symbol separators::

            sage: X = SFT(matrix(2, 2, [2,0, 1,1]), symb_sep=',')
            sage: Y = SFT(matrix(2, 2, [1,1, 1,0]), format='vertex', symb_sep=' ')
            sage: X.cartesian_product(Y)
            A 1-step SFT over [(0, 0), (0, 1), (1, 0), (1, 1), (2, 0), (2, 1), (3, 0), (3, 1)].
            sage: X = SFT(matrix(2, 2, [2,0, 1,1]), [0, 1, 2, 33], symb_sep=',')
            sage: X.cartesian_product(Y)
            Traceback (most recent call last):
            ...
            ValueError: symbols contain both symbol separators

        AUTHORS:

        - Sebastian Donoso (2011): initial version
        - Sebastian Barbieri (2012): documentation and forbidden words changes
        - Michael Schraudner (2012): included support for format, symb_sep and
          name, mayor cleanup, mschraudner@dim.uchile.cl
        """
        if not isinstance(other, SFT):
            raise TypeError("other must be an SFT")
        if not (name is None or isinstance(name, basestring)):
            raise TypeError("name must be a string or None")
        if self._is_empty:
            return self
        if other._is_empty:
            return other
        alphabet = [(a, b) for a in self.alphabet() for b in other.alphabet()]
        alph_str = [str(a) for a in alphabet]
        if not (self._symb_sep == "" or
                filter(lambda s: self._symb_sep in s, alph_str) == []):
            if not (other._symb_sep == "" or
                    filter(lambda s: other._symb_sep in s, alph_str) == []):
                symb_sep = ""
            else:
                symb_sep = other._symb_sep
        else:
            symb_sep = self._symb_sep
        if symb_sep == "":
            lenth = [len(s) for s in alph_str]
            if not lenth == [] and not min(lenth) == max(lenth):
                raise ValueError("symbols contain both symbol separators")
        fwords = []
        for u in self._fwords:  # construct pairs (u,w) as above
            for w in other._allwords(len(u)):
                fwords.append([(u[i], w[i]) for i in range(len(u))])
        for w in other.forbidden_words('list'):
            for u in self._allwords(len(w)):
                fwords.append([(u[i], w[i]) for i in range(len(w))])
        # have to create a new SFT ???
        return SFT(fwords, alphabet, format=self._format, symb_sep=symb_sep,
                   name=name)

## AUXILIARY_METHODS ##

    def _allwords(self, n):
        r"""
        Helper function that generates a list containing all words of size
        ``n`` of the full-shift on alphabet ``self._alph``.

        INPUT:

        - ``n`` - a positive integer

        OUTPUT:

        - a list of words in list format

        AUTHORS:

        - Sebastian Barbieri (2012): initial version
        """
        if not n in NN:
            raise TypeError("n must be a positive integer")
        l = [[a] for a in self._alph]
        for i in range(n-1):
            l = [a + [b] for a in l for b in self._alph]
        return l

    def _check_alph(self, n):
        r"""
        Helper function to check certain features (no repetitions, no symbol
        separator in symbols etc.) of the alphabet.

        INPUT:

        - ``n`` - a non-negative integer specifying the size of the expected
          alphabet. If one of the checks fails this number is used to create
          the standard alphabet ``{0, 1, ..., n-1}``.

        AUTHORS:

        - Michael Schraudner (2012): initial version, mschraudner@dim.uchile.cl
        """
        alph_str = [str(a) for a in self._alph]
        if not (self._symb_sep == "" or
                filter(lambda s: self._symb_sep in s, alph_str) == []):
            warnings.warn("symbols contain the symbol separator symb_sep. " +
                          "Will create a standard alphabet.")
            self._alph = range(n)
            if n > 10:
                self._symb_sep = '.'
            return
        if not len(self._alph) == len(set(alph_str)):
            warnings.warn("alphabet contains repetitions. Will create a " +
                          "standard alphabet.")
            self._alph = range(n)
            if n > 10:
                self._symb_sep = '.'
            return
        if self._symb_sep == "":
            lenth = [len(s) for s in alph_str]
            if not lenth == [] and not min(lenth) == max(lenth):
                raise ValueError("need a valid symbol separator symb_sep")

    def _create_empty(self):
        r"""
        Helper function to create the empty shift. Sets all internal variables
        of class SFT to hold the data of the empty vertex shift.

        AUTHORS:

        - Sebastian Donoso (2011): initial version
        - Sebastian Barbieri (2012): changed matrix representation
        """
        self._name = "The empty shift."
        self._is_empty = True
        self._alph = []
        self._format = 'vertex'
        self._fwords = []
        self._graph = DiGraph()
        self._matrix = matrix()
        self._symb_sep = ""

    def _create_edge_label_DiGraph(self):
        r"""
        Helper function that uses ``self._matrix`` to create a digraph whose
        edges are labeled by the elements of ``self._alph`` and whose vertices
        are generically numbered by 0,1,2,...
        This function also fills in the forbidden words using the not allowed
        transitions specified by ``self._matrix``.

        AUTHORS:

        - Sebastian Donoso (2011): initial version
        - Michael Schraudner (2012): cleanup, added forbidden words part,
          mschraudner@dim.uchile.cl
        """
        dic = {}
        z = 0
        for i in range(self._matrix.nrows()):
            d = {}
            for j in range(self._matrix.ncols()):
                d[j] = [str(a) for a in self._alph[z:z + self._matrix[i, j]]]
                z += self._matrix[i, j]
            dic[i] = d
        self._graph = DiGraph(dic, loops=True, multiedges=True)
        row_sums = [sum(self._matrix.row(i))
                                      for i in range(self._matrix.nrows())]
        row_sums = [sum(row_sums[:i]) for i in range(self._matrix.nrows()+1)]
        z = 0
        self._fwords = []
        for i in range(self._matrix.nrows()):
            for j in range(self._matrix.ncols()):
                self._fwords.extend([[a, b]
                                for a in self._alph[z:z + self._matrix[i, j]]
                                for b in self._alph[:row_sums[j]]])
                self._fwords.extend([[a, b]
                                for a in self._alph[z:z + self._matrix[i, j]]
                                for b in self._alph[row_sums[j+1]:]])
                z += self._matrix[i, j]

    def _create_vertex_label_DiGraph(self):
        r"""
        Helper function that uses ``self._matrix`` to create a digraph whose
        vertices are labeled by the elements of ``self._alph`` and whose
        edges are labeled by the words formed as the concatenation
        ``"initial vertex symb_sep terminal vertex"``.
        This function also fills in the forbidden words using the not allowed
        transitions specified by ``self._matrix``.

        AUTHORS:

        - Sebastian Donoso (2011): initial version
        - Michael Schraudner (2012): cleanup, added forbidden words part,
          mschraudner@dim.uchile.cl
        """
        dic = {}
        for i in range(self._matrix.nrows()):
            d = {}
            for j in range(self._matrix.ncols()):
                if self._matrix[i, j] == 1:
                    d[str(self._alph[j])] = [str(self._alph[i]) +
                                             self._symb_sep +
                                             str(self._alph[j])]
            dic[str(self._alph[i])] = d
        self._graph = DiGraph(dic, loops=True)
        self._fwords = [[self._alph[i], self._alph[j]]
                                        for i in range(self._matrix.nrows())
                                        for j in range(self._matrix.ncols())
                                        if self._matrix[i, j] == 0]

    def _create_fwords_DiGraph(self):
        r"""
        Helper function that uses ``self._alph`` and ``self._fwords`` to
        create the digraph ``self._graph`` whose vertices are words of length
        ``N`` and whose edges are words of length ``N+1`` over the alphabet
        ``self._alph``, where ``N`` is the length of the longest element in
        ``self._fwords``.
        This function also defines ``self._matrix`` as the graph's adjacency
        matrix.

        AUTHORS:

        - Michael Schraudner (2012): initial version, mschraudner@dim.uchile.cl
        """
        if self._fwords == []:
            if self._alph == []:
               self._create_empty()
            else:
                ## full-shift as vertex or edge shift ##
                self._is_empty = False
                self._name = "full-shift on %d symbols" % (len(self._alph))
                if self._format == 'vertex':
                    self._matrix = matrix(len(self._alph), len(self._alph),
                                          [1] * len(self._alph)**2)
                    dic = {}
                    for i in range(len(self._alph)):
                        d = {}
                        for j in range(len(self._alph)):
                            d[str(self._alph[j])] = [str(self._alph[i]) +
                                                     self._symb_sep +
                                                     str(self._alph[j])]
                        dic[str(self._alph[i])] = d
                    self._graph = DiGraph(dic, loops=True)
                else:
                    self._matrix = matrix(1, 1, [len(self._alph)])
                    self._graph = DiGraph({0: {0: [str(i) for i in self._alph]
                                              }}, loops=True, multiedges=True)
        else:
            ## edge or vertex shift ##
            max_len = max([len(w) for w in self._fwords])-1
            vert = [[a] for a in self._alph]
            for m in range(1, max_len):
                vert = [w+[a] for w in vert for a in self._alph
                              if filter(lambda u: w[1-len(u):] == u[:-1] and
                                        a == u[-1], self._fwords) == []]
            dic = {} # create set of edges and labels
            for i in vert:
                d = {}
                for j in vert:
                    if i[1:] == j[:-1]:
                        w = i + j[-1:]
                        if not w in self._fwords:
                            d[self._symb_sep.join([str(a) for a in j])
                              ] = [self._symb_sep.join([str(a) for a in w])]
                dic[self._symb_sep.join([str(a) for a in i])] = d
            self._graph = DiGraph(dic, loops=True, multiedges=True)
            self._matrix = self._graph.adjacency_matrix()
            if all([x == 0 for x in self._matrix.eigenvalues()]):
                self._create_empty()
            else:
                self._is_empty = False


class SFTGenerators(SageObject):
    r"""
    A class providing constructors for several families of common (well known)
    Subshifts of Finite Type.

    A list of all SFTs in this database is available via tab completion.
    Within a SAGE session, type sfts. (do not press 'Enter' and do not forget
    the final period '.') and then hit the 'Tab' key to see which SFTs are
    available.

    The constructors currently included in this class are::

      "Full"-shifts
      "Fibonacci"-shift
      "RLL"-shifts (Run-Length-Limited-shifts) and multiple spaced RLL-shifts
      "WhistleFree"-shifts
      "S_Gap"-shifts (of finite type)

    The docstring of each constructor provides information about how to
    generate corresponding SFTs in each family.

    AUTHORS:

    - Michael Schraudner (2011): initial version of class sfts = SFTGenerators,
      mschraudner@dim.uchile.cl
    - Sebastian Barbieri (2012): documentation and some new families added
    - Michael Schraudner (2012): some modifications (Word support) and cleanup
      to include class sfts into SAGE
    """

    def Full(self, n):
        r"""
        The full-shift on ``n`` symbols. Its elements are all bi-infinite
        sequences over the alphabet ``{0,1,...,n-1}``. This shift can also be
        seen as representing the set of all functions from ``\ZZ`` into the
        alphabet ``{0,1,...,n-1}``.

        INPUT:

        - ``n`` - a non-negative integer specifying the size of the
          full-shift's alphabet

        OUTPUT:

        - a (full-)Shift of Finite Type

        EXAMPLES::

            sage: X = sfts.Full(3); X
            The full-shift on 3 symbols of order 0 over [0, 1, 2].
            sage: X.entropy()
            log(3)

        ::

            sage: X = sfts.Full(0); X
            The empty shift.

        AUTHORS:

        - Michael Schraudner (2011): initial version, mschraudner@dim.uchile.cl
        - Sebastian Barbieri (2012): documentation
        """
        if n in NN:
            return SFT([], range(n), name="full-shift on %d symbols" % (n))
        else:
            raise ValueError("n must be a non-negative integer")

    def Fibonacci(self):
        r"""
        The Fibonacci shift, also known as the golden mean shift. It consists
        of all bi-infinite sequences over the alphabet ``{0,1}`` that do not
        contain two consecutive ones, i.e. the word ``11`` is forbidden.

        This shift is one of the most basic examples of an SFT. It has many
        interesting properties which also explain where its names come from,
        such as the fact that its entropy is the logarithm of the golden mean,
        ``log((1+sqrt(5))/2)`` while the number of admissible words of length
        ``n`` is given by the ``(n+2)``-th Fibonacci number.

        OUTPUT:

        - a (particular) Subshift of Finite Type

        EXAMPLES::

            sage: X = sfts.Fibonacci(); X
            The Fibonacci shift of order 1 over [0, 1].
            sage: X.forbidden_words()
            ['11']
            sage: exp(X.entropy())
            1.61803398874989

        AUTHORS:

        - Michael Schraudner (2011): initial version, mschraudner@dim.uchile.cl
        - Sebastian Barbieri (2012): documentation
        """
        return SFT([[1,1]], [0, 1], format='vertex', name="Fibonacci shift")

    def RLL(self, d, k, s=1):
        r"""
        The ``(d,k,s)``-Run-Length-Limited shift. For the default value of
        ``s=1`` it consists of all bi-infinite sequences over the alphabet
        ``{0,1}`` such that every sequence contains infinitely many ``1``\ s
        in its left and its right half, and such that the number of
        consecutive ``0``\ s between two occurrences of the symbol ``1`` lies
        between ``d`` and ``k`` (inclusive).

        If the parameter ``s`` is specified to be a natural number other than
        its default value ``s=1``, the length of each block of consecutive
        ``0``\ s between two occurrences of the symbol ``1`` has to be an
        element in the arithmetic progression ``d``, ``d+s``, ``d+2s`` ...
        ``d+\floor{(k-d)/s}*s`` (i.e. ``d`` plus some non-negative multiple of
        ``s`` not exceeding ``k``). Those subshifts are known under the name
        multiple spaced RLLs.

        (Multiple spaced) RLL shifts appear in the theory of communication
        channels (both in telecommunication with limited bandwidth as well as
        in the storage of arbitrary data on hard drives and magnetic or
        optical disks) where certain technical conditions (interference and
        clock drift) impede consecutive occurrences of ``1``\ s with either
        too small or too large separation.

        INPUT:

        - ``d`` - a non-negative integer specifying the minimum amount of
          consecutive ``0``\ s allowed between two ``1``\ s. If ``d=0``, the
          shift contains the fixed point consisting only of ``1``\ s

        - ``k`` - a non-negative integer or ``+Infinity`` which has to be
          greater than or equal to ``d``. It specifies the maximal allowed
          length of a block of consecutive ``0``\ s without seeing any symbol
          ``1``

        - ``s`` - (default: 1) a positive integer. If specified, the blocks of
          ``0``\ s have to have lengths ``d``, ``d+s``, ``d+2s`` ... up to at
          most ``k``

        OUTPUT:

        - a Subshift of Finite Type (representing the (d,k)-RLL or the
          (d,k,s)-RLL)

        REFERENCES:

        - 'D. Lind, B. Marcus: An introduction to symbolic dynamics and coding. Chapter 1'
        - 'P. Funk: Run-Length-Limited Codes with Multiple Spacing, IEEE Transactions on Magnetics, 18/2 (1982)'

        EXAMPLES:

        #. Some run-lenth-limited shifts::

                sage: X = sfts.RLL(1, 3); X
                The (1,3)-RLL shift of order 3 over [0, 1].
                sage: X.forbidden_words()
                ['0000', '11']
                sage: X.entropy()
                0.382245085840036

           ::

                sage: X = sfts.RLL(2, +Infinity); X
                The (2,+Infinity)-RLL shift of order 2 over [0, 1].
                sage: X.forbidden_words()
                ['11', '101']

           ::

                sage: X = sfts.RLL(2, 6, 2); X
                The (2,6,2)-RLL shift of order 6 over [0, 1].
                sage: X.forbidden_words()
                ['0000000', '11', '101', '10001', '1000001']

        #. For some choices of parameters the actual order of the RLL shift is
           smaller than expected::

            sage: X = sfts.RLL(2, 5, 2); X
            The (2,5,2)-RLL shift of order 4 over [0, 1].
            sage: X.forbidden_words()
            ['000000', '11', '101', '10001', '1000001']
            sage: X.first_offenders()
            (['11', '101', '00000', '10001'], [0, 1])

        AUTHORS:

        - Michael Schraudner (2011): initial version, mschraudner@dim.uchile.cl
        - Sebastian Barbieri (2012): documentation, added multiple spaced RLLs
        """
        if d in NN and (k in NN or k == +Infinity) and d <= k:
            if s in NN and s > 0:
                if k < +Infinity:
                    fwords = [[0]*(k+1)]
                else:
                    fwords = []
                fwords.extend([[1] + [0]*i + [1] for i in range(d)])
                if s == 1:
                    return SFT(fwords, [0, 1],
                               name="(%s,%s)-RLL shift" % (d, k))
                else:
                    fwords.extend([[1] + [0]*(d+i) + [1]
                                   for i in range(k-d+1) if not (i%s == 0)])
                    return SFT(fwords, [0, 1],
                               name="(%s,%s,%s)-RLL shift" % (d, k, s))
            else:
                raise ValueError("s must be a positive integer")
        else:
            raise ValueError("parameters d, k must be non-negative integers" +
                             " (k can also be +Infinity) and k must not be " +
                             "less than d")

    def WhistleFree(self, L):
        r"""
        The ``L``-Whistle-Free shift. It consists of all bi-infinite sequences
        over the alphabet ``{0,1}`` that do not contain any of the words
        w_0^n_0, w_1^n_1, ..., w_k^n_k specified by the list ``L`` formed by
        the pairs ``(w_i, n_i)`` (for ``i = 0, 1, ..., k``).

        Those shifts appear in the theory of communication channels where
        certain technical conditions preclude the possibility to transmit
        signals which contain long words with (particular) small periods.

        INPUT:

        - ``L`` - a list of 2-tuples, i.e. ``L = [(w_i, n_i) for i=0,...,k]``
          where each ``w_i`` is a word over the alphabet ``{0,1}`` which may
          be given as a string, a list/tuple or a SAGE Word and ``n_i`` is a
          non-negative integer

          The list ``L`` determines the forbidden words of the corresponding
          general Whistle-Free shift.

        OUTPUT:

        - a Subshift of Finite Type (representing the general
          ``L``-Whistle-Free shift)

        EXAMPLES:

        #. The standard 2-Whistle-Free shift, forbidding the words ``0101``
           and ``1010``::

            sage: X = sfts.WhistleFree([('01',2), ('10',2)]); X
            The generalized Whistle-Free shift of order 3 over [0, 1].
            sage: X.forbidden_words()
            ['0101', '1010']
            sage: X.entropy()
            0.609377863436006

        #. A generalized Whistle-Free shift::

            sage: X = sfts.WhistleFree([(Word('1'),5), (Word('01'), 3), (Word('101'), 3)]); X
            The generalized Whistle-Free shift of order 8 over [0, 1].
            sage: X.forbidden_words()
            ['11111', '010101', '101101101']

        #. Testing::

            sage: X = sfts.WhistleFree([('0123',2)])
            Traceback (most recent call last):
            ...
            ValueError: forbidden words in L have to be over the alphabet [0,1]

        AUTHORS:

        - Michael Schraudner (2011): initial version, mschraudner@dim.uchile.cl
        - Sebastian Barbieri (2012): documentation, extended to general WFS
        - Michael Schraudner (2012): error handling, support for list and Word
        """
        if not isinstance(L, list) or any([not (isinstance(l, tuple) and
                                                len(l) == 2) for l in L]):
            raise TypeError("L must be a list of 2-tuples")
        for l in L:
            try:
                symbs = filter(lambda a: not str(a) in ["0","1"],l[0])
            except:
                raise TypeError("2-tuples in L must contain an iterable as " +
                                "their first component")
            if len(symbs):
                raise ValueError("forbidden words in L have to be over the " +
                                 "alphabet [0,1]")
        try:
            fwords = [map(lambda a: str(a) == "1" and 1 or 0, l[0]) * l[1]
                      for l in L]
        except:
            raise ValueError("2-tuples in L must be of the form (string, " +
                             "non-negative integer) or (list, non-negative " +
                             "integer) or (Word, non-negative integer)")
        return SFT(fwords, [0, 1], name="generalized Whistle-Free shift")

    def S_Gap(self, S):
        r"""
        The ``S``-Gap shift. It consists of all bi-infinite sequences over the
        alphabet ``{0,1}`` such that the words ``1 0^i 1`` with ``i`` in ``S``
        are forbidden.

        INPUT:

        - ``S`` - a finite list of non-negative integers

        OUTPUT:

        - a Subshift of Finite Type (representing the S-Gap shift)

        REFERENCES:

        - 'D. Lind, B. Marcus: An introduction to symbolic dynamics and coding. Chapter 1'

        EXAMPLES::

            sage: sfts.S_Gap([1, 0, 4])
            The [0, 1, 4]-Gap shift of order 5 over [0, 1].

        AUTHORS:

        - Michael Schraudner (2012): suggested, mschraudner@dim.uchile.cl
        - Sebastian Barbieri (2012): initial version, documentation
        """
        if not isinstance(S, list):
            raise ValueError("S must be a list")
        if any([not x in NN for x in S]):
            raise ValueError("S must be a list of non-negative integers")
        S = copy(S)
        S.sort()
        fwords = [[1] + [0]*i + [1] for i in S]
        return SFT(fwords, [0, 1], name="%s-Gap shift" % (S))

sfts=SFTGenerators()
