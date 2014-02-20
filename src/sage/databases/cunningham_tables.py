r"""
Cunningham tables of factors.

The Cunningham project is a long-running project to find factors of numbers of the form
`b^k \pm 1`, where `b` is an element of `\{2, 3, 5, 6, 7, 10, 11, 12\}`.  The project
is currently hosted at [SSW]_, and the tables in Sage are generated from the
machine-readable tables at [JDe]_.

In Sage, the Cunningham tables are accessed as follows::

    sage: CD = CunninghamDatabase(); CD
    Database of factors of numbers of the form b^k + 1 and b^k - 1
    for b in 2,3,5,6,7,10,11,12

There are two main interfaces for the factors in the tables.

List Interface
--------------

In some ways, the Cunningham tables behave as a sorted list::

    sage: len(CD)
    26242
    sage: CD[1000]
    21841
    sage: for d in CD:
    ...       if d < 2500:
    ...           print d
    2467
    2473
    2477
    sage: [i for i in range(len(CD)) if CD[i].divides(2^6010 - 1)]
    [118, 5371, 12351, 20515, 23617, 24937]

Here the database is serving as an interface to the underlying list of
all factors appearing in the Cunningham tables.  You can access the underlying list::

    sage: type(CD.list)
    <type 'list'>

In fact, not all factors in the Cunningham tables appear here: 2 certainly divides
every `3^k \pm 1`.  The Sage interface omits factors less than a certain bound::

    sage: CD.bound
    2441

These small factors can be accessed by :meth:`~CunninghamDatabase.list_primes`::

    sage: CD.list_primes()[0]
    2

The composites are also accessible.  If you're looking for an integer factoring challenge,
try factoring one of these::

    sage: CD.list_composites()[0].nbits()
    545

Dictionary Interface
--------------------

The list interface above gives no information about which expressions `b^k \pm 1` a given
factor divides.  So the Cunningham tables also have an interface that allows access
to this data.

    sage: v = CD[2, -140]; v
    (0, [86171, 122921, 7416361, 47392381])
    sage: factor(2^140 - 1)
    3 * 5^2 * 11 * 29 * 31 * 41 * 43 * 71 * 113 * 127 * 281 * 86171 * 122921 * 7416361 * 47392381

When passed a triple of integers ``(b, k, s)``, the database returns a pair ``(f, L)``.
`L` consists of a list of factors of `b^k + s`, and `f` is a flag whose bits indicate
which entries of `L` are composite.

    sage: v = CD[2, -3788]; v[0].bits()
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1]
    sage: [0 if d.is_pseudoprime() else 1 for d in v[1]]
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0]
    sage: len(v[1])
    20
    sage: v[1][16].nbits()
    701

You can access the underlying dictionary::

    sage: type(CD.dict)
    <type 'dict'>

WARNING:

You should not modify either the ``CunninghamDatabase``\ 's :attr:`list` or :attr:`dict`: they
are provided without copying for efficiency's sake.

REFERENCES:

.. [Wik] Cunningham Project, Wikipedia.  http://en.wikipedia.org/wiki/Cunningham_project

.. [SSW] The Cunningham Project.  http://homes.cerias.purdue.edu/~ssw/cun/index.html

.. [JDe] Machine-readable Cunningham Tables.  http://cage.ugent.be/~jdemeyer/cunningham/

AUTHORS:

- Yann Laigle-Chapuy (2009-10-18): initial version
- David Roe (2011-12-4): Updated data and format
"""


import os
from sage.misc.cachefunc import cached_method
from sage.rings.integer import Integer
from sage.structure.sage_object import load, SageObject
from sage.rings.fast_arith import prime_range
from sage.misc.superseded import deprecated_function_alias
from sage.structure.unique_representation import UniqueRepresentation
import bisect
from sage.misc.misc import SAGE_SHARE

# We do this here to avoid CunninghamDatabase() having an `a` attribute.
_base_list = [Integer(a) for a in [2,3,5,6,7,10,11,12]]

class CunninghamDatabase(UniqueRepresentation, SageObject):
    r"""
    An object storing prime factors of numbers of the form `b^k + 1`
    and `b^k - 1` for `b \in \{2,3,5,6,7,10,11,12\}`.

    Data is generated from http://cage.ugent.be/~jdemeyer/cunningham/

    EXAMPLES::

        sage: CD = CunninghamDatabase(); CD
        Database of factors of numbers of the form b^k + 1 and b^k - 1
        for b in 2,3,5,6,7,10,11,12

    :attr:`list` stores a sorted list of all of the factors stored::

        sage: len(CD.list)
        26242
        sage: CD.list[0]
        2467
        sage: CD.list[-1].nbits()
        1872

    :attr:`dict` stores a dictionary with information on which factors occur
    for what bases and exponents.  The keys of the dictionary are triples
    ``(b, k, s)`` where

    - `b` is a Cunningham base: 2, 3, 5, 6, 7, 10, 11 or 12
    - `s` is 1 or -1
    - `k` is an exponent so that the factorization of `b^k + s` is
      found in the Cunningham tables at
      http://cage.ugent.be/~jdemeyer/cunningham/

    The value attached to ``(b, k, s)`` is a pair ``(f, L)`` where

    - `L` is a list of relatively prime divisors of `b^k + s`, so that
      `b^k + s` is the product of the entries of `L` together with
      primes up to :attr:`bound`.
    - `f` is a flag giving which entries of `L` are BPSW pseudoprimes.  If
    `f = 0`, then all entries of `L` are pseudoprimes.  More generally,
    the positions of the 1s in the binary expansion of `f` correspond
    to the composite values in `L`.

    ::

        sage: flag, plist = CD.dict[2,544]
        sage: (2^544 + 1) / prod(plist)
        641
        sage: flag
        0
        sage: all([p.is_pseudoprime() for p in plist])
        True

    The factorizations are not necessarily into primes::

        sage: flag, plist = CD.dict[10,-926]
        sage: flag
        48
        sage: plist[4].is_pseudoprime()
        False
        sage: plist[5].is_pseudoprime()
        False
        sage: plist[4].nbits()
        1456

    Note that small primes (below :attr:`bound`) are omitted from :attr:`list`
    and :attr:`dict`

        sage: 31 in CD.list
        False

    TESTS:

    We check that the only prime factors not included in the dict are
    smaller than the cunningham bound.::

        sage: from sage.rings.factorint import factor_trial_division
        sage: for ky in CD.dict.keys():
        ...       b, ks = ky
        ...       k, s = ks.abs(), ks.sign()
        ...       n = b^k + s
        ...       rem = ZZ(n/prod(CD.dict[ky][1]))
        ...       F = factor_trial_division(rem, CD.bound)
        ...       assert(F.prod() == rem)

    We check that no factor appears twice in one row of the cunningham table:
    this assumption is used in :meth:`~sage.rings.factorint.factor_cunningham`::

        sage: for ky in CD.dict.keys():
        ...       for i in range(len(CD.dict[ky][1]) - 1):
        ...           assert(CD.dict[ky][1][i] < CD.dict[ky][1][i+1])
    """

    bound = Integer(load(os.path.join(SAGE_SHARE, "cunningham_tables", "cunningham_bound.sobj")))
    base_list = _base_list

    def __init__(self):
        """
        Initialization.

        TESTS::

            sage: TestSuite(CunninghamDatabase()).run()
        """
        self.list = map(Integer, load(os.path.join(SAGE_SHARE, "cunningham_tables", "cunningham_list.sobj")))
        D = load(os.path.join(SAGE_SHARE, "cunningham_tables", "cunningham_dict.sobj"))
        self.dict = {}
        for ky, val in D.iteritems():
            self.dict[Integer(ky[0]),Integer(ky[1])] = Integer(val[0]), [self.list[a] for a in val[1]]

    def _repr_(self):
        """
        String representation.

        EXAMPLES::

            sage: CunninghamDatabase() # indirect doctest
            Database of factors of numbers of the form b^k + 1 and b^k - 1
            for b in 2,3,5,6,7,10,11,12
        """
        return "Database of factors of numbers of the form b^k + 1 and b^k - 1\nfor b in 2,3,5,6,7,10,11,12"

    def __getitem__(self, v):
        """
        If `v` is an integer or of length 1, returns ``self.list[v]``
        Otherwise, returns ``self.dict[v]``.

        EXAMPLES::

            sage: CD = CunninghamDatabase()
            sage: CD[2]
            2477
            sage: CD[(2,)]
            2477
            sage: CD[2, 71]
            (0, [56409643, 13952598148481])
        """
        if not isinstance(v, (list, tuple)):
            return self.list[v]
        elif len(v) == 1:
            return self.list[v[0]]
        elif len(v) == 2:
            return self.dict[v]
        else:
            raise ValueError("v must be an integer or tuple of length 1 or 2")

    def __len__(self):
        """
        Length of the list of Cunningham factors.

        EXAMPLES::

            sage: CD = CunninghamDatabase()
            sage: len(CD)
            26242
        """
        return len(self.list)

    def __contains__(self, x):
        """
        Containment testing.

        Should behave identically to ``x in self.list``.

        EXAMPLES::

            sage: CD = CunninghamDatabase()
            sage: 2441 in CD
            False
            sage: 2477 in CD
            True
        """
        i = bisect.bisect_left(self.list, x)
        return (i < len(self.list)) and x == self.list[i]

    def count(self, x):
        """
        Count.

        Should behave identically to ``self.list.count(x)``.

        EXAMPLES::

            sage: CD = CunninghamDatabase()
            sage: CD.count(2441)
            0
            sage: CD.count(2477)
            1
        """
        i = bisect.bisect_left(self.list, x)
        return 1 if ((i < len(self.list)) and x == self.list[i]) else 0

    def index(self, x, start=None, stop=None):
        """
        Index in list.

        Should behave identically to ``self.list.index(x, start, stop)``.

        EXAMPLES::

            sage: CD = CunninghamDatabase()
            sage: CD.index(2441)
            Traceback (most recent call last):
            ...
            ValueError: x not in list
            sage: CD.index(2477)
            2
            sage: CD.index(293871)
            Traceback (most recent call last):
            ...
            ValueError: x not in list
            sage: CD.index(CD[57])
            57
        """
        if start is None:
            start = 0
        if stop is None:
            stop = len(self.list)
        i = bisect.bisect_left(self.list, x, start, stop)
        if (i < len(self.list)) and x == self.list[i]:
            return i
        raise ValueError("x not in list")

    def __reversed__(self):
        """
        Provide a reversed iterator over ``self.list``.

        EXAMPLES::

            sage: CD = CunninghamDatabase()
            sage: for d in reversed(CD): #indirect doctest
            ...       if d.divides(2^140 + 1):
            ...           print CD.index(d)
            15591
            4723
            1618
        """
        return self.list.__reversed__()

    def items(self):
        """
        List of ``self.dict``\ 's (key, value) pairs, as 2-tuples

        EXAMPLES::

            sage: CD = CunninghamDatabase()
            sage: for ky, v in CD.items():
            ...       if ky[1].abs() == 5:
            ...           print prod(v[1])
            3221
            9091
            2801
            13421
            22621
            19141
        """
        return self.dict.items()

    def keys(self):
        """
        List of the keys ``(b, ks)`` that have divisors in the database above the bound.

        EXAMPLES::

            sage: len(CunninghamDatabase().keys())
            11249
        """
        return self.dict.keys()

    def values(self):
        """
        List of the values ``(f, L)`` in the dictionary.

        EXAMPLES::

            sage: len(CunninghamDatabase().values())
            11249
        """
        return self.dict.values()

    def iteritems(self):
        """
        An iterator over the items in the dictionary.

        EXAMPLES::

            sage: CD = CunninghamDatabase()
            sage: for ky, v in CD.iteritems():
            ...       if ky[1] == 5:
            ...           print prod(v[1])
            9091
            13421
            19141
        """
        return self.dict.iteritems()

    def iterkeys(self):
        """
        An iterator over the items in the dictionary.

        EXAMPLES::

            sage: CD = CunninghamDatabase()
            sage: for ky in CD.iterkeys():
            ...       if ky[1] == 5:
            ...           print ky[0]**(ky[1].abs()) + ky[1].sign()
            100001
            161052
            248833
        """
        return self.dict.iterkeys()

    def itervalues(self):
        """
        An iterator over the values in the dictionary.

        EXAMPLES::

            sage: CD = CunninghamDatabase()
            sage: L = []
            sage: for v in CD.itervalues():
            ...       if sum(v[0].bits()) == 3:
            ...           L.append(prod(v[1]).nbits())
            sage: L
            [4345, 4433, 4645, 4601, 4261, 4621, 4129, 4033, 3785, 4391, 4596, 4385, 3961, 4528, 4609, 4252, 4369, 4816]
        """
        return self.dict.itervalues()

    def get(self, ky, v=None):
        """
        If ``ky`` is in the dictionary, return the associated value.  Otherwise, return `v`, or
        ``(1, [ky[0]**ky[1] + ky[2]])`` if `v` is ``None``.

        EXAMPLES::

            sage: CD = CunninghamDatabase()
            sage: CD.get((2, 10))
            (1, [1025])
            sage: CD.dict.get((2, -91))
            (0, [8191, 112901153, 23140471537])
        """
        if v is None and isinstance(ky, tuple) and len(ky) == 2 and all([isinstance(a, (Integer, int)) for a in ky]):
            v = (Integer(1), [ky[0]**ky[1].abs() + ky[1].sign()])
        return self.dict.get(ky, v)

    @cached_method
    def bdict(self, b):
        """
        Returns the subdictionary of :attr:`dict` corresponding to the base `b`.

        Note that `b` is omitted from the keys in the resulting dictionary.

        INPUT:

        - `b` -- one of 2, 3, 5, 6, 7, 10, 11, 12

        OUTPUT:

        - A dictionary with keys ``ks`` and values ``(f, L)``, where

          - `k` is the exponent
          - `s` is 1 or -1
          - `L` is a list of relatively prime factors of `b^k + s`
          - `f` is a flag: the `i^{th}` bit of `f` is set iff the `i^{th}` element
            of `L` is composite.

        EXAMPLES::

            sage: CD = CunninghamDatabase()
            sage: D = CD.bdict(2); D.has_key(71)
            True
            sage: D[71]
            (0, [56409643, 13952598148481])
        """
        bkeys = self.bkeys(b)
        cdict = self.dict
        ans = {}
        for ky in bkeys:
            ans[ky[1]] = cdict[ky]
        return ans

    @cached_method
    def bkeys(self, b):
        """
        Returns a sorted list of the keys of ``self.dict`` corresponding to
        the base `b`.

        INPUT:

        - `b` -- one of 2, 3, 5, 6, 7, 10, 11, 12

        OUTPUT:

        - A sorted list of keys ``ks`` corresponding to factorizations in
          the Cunningham tables with base `b`.

        EXAMPLES::

            sage: CD = CunninghamDatabase()
            sage: CD.bkeys(5)[33]
            (5, -845)
        """
        ans = []
        for ky in self.dict.iterkeys():
            if ky[0] == b:
                ans.append(ky)
        ans.sort()
        return ans

    @cached_method
    def list_primes(self):
        """
        Sorted list of the primes appearing in the Cunningham table,
        including the primes below :attr:`bound` that are ommitted from :attr:`list`.

        EXAMPLES::

            sage: L = CunninghamDatabase().list_primes()
            sage: len(L)
            26186
            sage: L[0]
            2
            sage: L[-1].nbits()
            1872

        The old ``cunningham_prime_factors`` has been deprecated::

            sage: from sage.databases.cunningham_tables import cunningham_prime_factors
            sage: cunningham_prime_factors()[0]
            doctest:...: DeprecationWarning: cunningham_prime_factors is deprecated. Please use sage.databases.cunningham_tables.CunninghamDatabase().list_primes instead.
            See http://trac.sagemath.org/7240 for details.
            2
        """
        L = []
        for v in self.dict.itervalues():
            bits = v[0].bits()
            factors = v[1]
            for i in range(len(factors)):
                if i >= len(bits) or not bits[i]:
                    L.append(factors[i])
        return prime_range(self.bound + 1) + sorted(list(set(L)))

    @cached_method
    def list_composites(self):
        """
        Sorted list of composite numbers occuring in the Cunningham tables.

        If you succeed in factoring one of these composites,
        you should write to Sam Wagstaff <ssw@cerias.purdue.edu>.

        EXAMPLES::

            sage: L = CunninghamDatabase().list_composites()
            sage: len(L)
            418
            sage: L[0].nbits()
            545
            sage: L[-1].nbits()
            1580
        """
        L = []
        for v in self.dict.itervalues():
            bits = v[0].bits()
            factors = v[1]
            for i in range(len(bits)):
                if bits[i]:
                    L.append(factors[i])
        return sorted(list(set(L)))

    def covering_keys(self, b, ks0):
        """
        Returns an unsorted list of keys in the Cunningham tables that are relevant
        to factoring a number whose exponent lies outside the scope of the tables.

        INPUT:

        - ``b`` -- integer, the desired base

        - ``ks0`` -- integer, the product of the desired exponent
          ``k0`` and the sign ``s0``: factors of `b^{k0} + s0` are
          sought.

        OUTPUT:

        - An unsorted list of keys ``(b, ks)`` so that any factor in the Cunningham
          tables appearing in the factorization of `b^{k0} + s0` (for an algebraic reason)
          will appear the factors associated to some key `(b, ks)`.

        EXAMPLES:

        If ``(b, ks0)`` is already a key for the dictionary, this function will just return
        one value.

            sage: CD = CunninghamDatabase()
            sage: CD.covering_keys(2, -2310)
            [(2, -2310)]
            sage: CD.dict.has_key((2, -2310))
            True

        Of course, if we lie outside the scope of the tables we don't usually catch all the factors::

            sage: L = CD.covering_keys(2, 3*5*7*11*13); L
            [(2, 1155), (2, 1001), (2, 715), (2, 455), (2, 429), (2, 273), (2, 195)]
            sage: N = 2^(3*5*7*11*13) + 1
            sage: N.nbits()
            15016
            sage: for ky in L:
            ...       for d in CD.dict[ky][1]:
            ...           N = N // (N.gcd(d))
            sage: N.nbits()
            11719
            sage: N.is_pseudoprime()
            False

        TESTS::

            sage: CD.covering_keys(2, 2*3*5*7*11*13)
            [(2, 2310), (2, 2002), (2, 1430), (2, 910), (2, 858), (2, 546), (2, 390)]
            sage: CD.covering_keys(2, -2*3*5*7*11*13)
            [(2, 1155), (2, 1001), (2, 715), (2, 455), (2, 429), (2, 273), (2, 195), (2, -1155), (2, -1001), (2, -715), (2, -455), (2, -429), (2, -273), (2, -195)]
            sage: CD.covering_keys(2, -3*5*7*11*13)
            [(2, -1155), (2, -1001), (2, -715), (2, -455), (2, -429), (2, -273), (2, -195)]
            sage: CD.covering_keys(2, -12400)
            [(2, 248), (2, 200), (2, 620), (2, 100), (2, 1550), (2, 775), (2, -775)]
        """
        b = Integer(b)
        ks0 = Integer(ks0)
        cdict = self.dict
        if cdict.has_key((b, ks0)):
            return [(b, ks0)]
        bkeys = self.bkeys(b)
        k0 = ks0.abs()
        L = []
        val, u = k0.val_unit(2)
        if ks0 < 0:
            # we can factor b^k - 1 formally if k factors
            if k0.is_prime():
                return L
            k = k0
            # b^(2^val * u) - 1 = (b^(2^(val-1) * u) + 1) * ... * (b^u + 1) * (b^u - 1)
            for w in range(val):
                k = k >> 1
                L.extend(self.covering_keys(b, k))
            # now k == u
            # For each divisor d of u,
            # b^u - 1 = (b^d - 1) * (b^(d * (u/d-1)) + ... + b^(d*2) + 1)
            found_divisors = []
            for d in reversed(u.divisors()[1:]):
                ky = (b, -d)
                if cdict.has_key(ky):
                    for D in found_divisors:
                        if d.divides(D):
                            break
                    else:
                        L.append(ky)
                        found_divisors.append(d)
        else:
            # For d odd, b^(xd) + 1 = (b^x + 1) * (b^(x(d-1)) - b^(x(d-2)) + ... + 1)
            found_divisors = []
            for d in u.divisors()[1:]:
                ky = (b, k0 // d)
                if cdict.has_key(ky):
                    for D in found_divisors:
                        if D.divides(d):
                            break
                    else:
                        L.append(ky)
                        found_divisors.append(d)
        return L

def _dep_list_primes():
    """
    Auxilliary function so that the deprecation warning prints correctly.

    EXAMPLES::

        sage: from sage.databases.cunningham_tables import _dep_list_primes
        sage: len(_dep_list_primes())
        26186
    """
    return CunninghamDatabase().list_primes()
_dep_list_primes.__name__ = 'CunninghamDatabase().list_primes'

cunningham_prime_factors = deprecated_function_alias(7240, _dep_list_primes)
cunningham_prime_factors.__name__ = 'cunningham_prime_factors'


