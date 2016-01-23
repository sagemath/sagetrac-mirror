"""
Necklaces

The algorithm used in this file comes from

- Sawada, Joe.  "A fast algorithm to generate necklaces with fixed content", Source
  Theoretical Computer Science archive Volume 301 , Issue 1-3 (May
  2003)
"""
#*****************************************************************************
#       Copyright (C) 2007 Mike Hansen <mhansen@gmail.com>,
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.lazy_attribute import lazy_attribute
from sage.combinat.composition import Composition
from sage.combinat.combinat import CombinatorialClass
from sage.arith.all import euler_phi,factorial, divisors, gcd
from sage.rings.integer import Integer
from sage.misc.all import prod
from sage.combinat.misc import DoublyLinkedList

def Necklaces(content):
    r"""
    Return the combinatorial class of necklaces with evaluation ``content``.

    A necklace is a list of integers that such that the list is
    the smallest lexicographic representative of all the cyclic shifts
    of the list.

    .. SEEALSO::

        :class:`LyndonWords`

    INPUT:

    - ``content`` -- a list of non-negative integers

    EXAMPLES::

        sage: Necklaces([2,1,1])
        Necklaces with evaluation [2, 1, 1]
        sage: Necklaces([2,1,1]).cardinality()
        3
        sage: Necklaces([2,1,1]).first()
        [1, 1, 2, 3]
        sage: Necklaces([2,1,1]).last()
        [1, 2, 1, 3]
        sage: Necklaces([2,1,1]).list()
        [[1, 1, 2, 3], [1, 1, 3, 2], [1, 2, 1, 3]]
        sage: Necklaces([0,2,1,1]).list()
        [[2, 2, 3, 4], [2, 2, 4, 3], [2, 3, 2, 4]]
        sage: Necklaces([2,0,1,1]).list()
        [[1, 1, 3, 4], [1, 1, 4, 3], [1, 3, 1, 4]]
    """
    return Necklaces_evaluation(content)

class Necklaces_evaluation(CombinatorialClass):
    """
    Necklaces with a fixed evaluation (content).

    INPUT:

    - ``content`` -- a list of non-negative integers
    """
    def __init__(self, content):
        r"""
        Initialize ``self``.

        TESTS::

            sage: N = Necklaces([2,2,2])
            sage: N == loads(dumps(N))
            True
        """
        if isinstance(content, Composition):
            self._content = content
        else:
            self._content = Composition(content)

    @lazy_attribute
    def e(self):
        """
        Deprecated in :trac:`17436`. Use :meth:`content` instead.

        TESTS::

            sage: N = Necklaces([2,2,2])
            sage: N.e
            doctest:...: DeprecationWarning: e attribute is deprecated. Use the content method instead
            See http://trac.sagemath.org/17436 for details.
            [2, 2, 2]
        """
        from sage.misc.superseded import deprecation
        deprecation(17436, 'e attribute is deprecated. Use the content method instead')
        return self._content

    def content(self):
        """
        Return the content (or evaluation) of the necklaces.

        TESTS::

            sage: N = Necklaces([2,2,2])
            sage: N.content()
            [2, 2, 2]
        """
        return self._content

    def __repr__(self):
        r"""
        TESTS::

            sage: repr(Necklaces([2,1,1]))
            'Necklaces with evaluation [2, 1, 1]'
        """
        return "Necklaces with evaluation %s"%self._content

    def __contains__(self, x):
        r"""
        Return ``True`` if ``x`` is the smallest word of all its cyclic shifts
        and the content vector of ``x`` is equal to ``content``.

        INPUT:

        - ``x`` -- a list of integers

        EXAMPLES::

            sage: [2,1,2,1] in Necklaces([2,2])
            False
            sage: [1,2,1,2] in Necklaces([2,2])
            True
            sage: [1,1,2,2] in Necklaces([2,2])
            True
            sage: [1,2,2,2] in Necklaces([2,2])
            False
            sage: all([ n in Necklaces([2,1,3,1]) for n in Necklaces([2,1,3,1])])
            True
            sage: all([ n in Necklaces([0,1,2,3]) for n in Necklaces([0,1,2,3])])
            True
        """
        xl = list(x)
        e = [0]*len(self._content)
        if len(xl) != sum(self._content):
            return False

        #Check to make sure xl is a list of integers
        for i in xl:
            if not isinstance(i, (int, Integer)):
                return False
            if i <= 0:
                return False
            if i > len(self._content):
                return False
            e[i-1] += 1

        #Check to make sure the evaluation is the same
        if e != self._content:
            return False

        #Check to make sure that x is lexicographically less
        #than all of its cyclic shifts
        cyclic_shift = xl[:]
        for i in range(len(xl) - 1):
            cyclic_shift = cyclic_shift[1:] + cyclic_shift[:1]
            if cyclic_shift < xl:
                return False

        return True

    def cardinality(self):
        r"""
        Return the number of integer necklaces with the evaluation ``content``.

        The formula for the number of necklaces of content `\alpha`
        a composition of `n` is:

        .. MATH::

            \sum_{d|gcd(\alpha)} \phi(d)
            \binom{n/d}{\alpha_1/d, \ldots, \alpha_\ell/d},

        where `\phi(d)` is the Euler `\phi` function.

        EXAMPLES::

            sage: Necklaces([]).cardinality()
            0
            sage: Necklaces([2,2]).cardinality()
            2
            sage: Necklaces([2,3,2]).cardinality()
            30
            sage: Necklaces([0,3,2]).cardinality()
            2

        Check to make sure that the count matches up with the number of
        necklace words generated.

        ::

            sage: comps = [[],[2,2],[3,2,7],[4,2],[0,4,2],[2,0,4]]+Compositions(4).list()
            sage: ns = [ Necklaces(comp) for comp in comps]
            sage: all( [ n.cardinality() == len(n.list()) for n in ns] )
            True
        """
        evaluation = self._content
        le = list(evaluation)
        if not le:
            return 0

        n = sum(le)

        return sum(euler_phi(j)*factorial(n/j) / prod(factorial(ni/j)
                    for ni in evaluation) for j in divisors(gcd(le))) / n

    def __iter__(self):
        r"""
        An iterator for the integer necklaces with evaluation ``content``.

        EXAMPLES::

            sage: Necklaces([]).list()    #indirect test
            []
            sage: Necklaces([1]).list()   #indirect test
            [[1]]
            sage: Necklaces([2]).list()   #indirect test
            [[1, 1]]
            sage: Necklaces([3]).list()   #indirect test
            [[1, 1, 1]]
            sage: Necklaces([3,3]).list() #indirect test
            [[1, 1, 1, 2, 2, 2],
             [1, 1, 2, 1, 2, 2],
             [1, 1, 2, 2, 1, 2],
             [1, 2, 1, 2, 1, 2]]
            sage: Necklaces([2,1,3]).list() #indirect test
            [[1, 1, 2, 3, 3, 3],
             [1, 1, 3, 2, 3, 3],
             [1, 1, 3, 3, 2, 3],
             [1, 1, 3, 3, 3, 2],
             [1, 2, 1, 3, 3, 3],
             [1, 2, 3, 1, 3, 3],
             [1, 2, 3, 3, 1, 3],
             [1, 3, 1, 3, 2, 3],
             [1, 3, 1, 3, 3, 2],
             [1, 3, 2, 1, 3, 3]]
        """
        if not self._content:
            return
        k = 0
        while not self._content[k]: # == 0
            k = k+1
        for z in _sfc(self._content[k:]):
            yield [x+1+k for x in z]


##############################
#Fast Fixed Content Algorithm#
##############################
def _ffc(content, equality=False):
    """
    EXAMPLES::

        sage: from sage.combinat.necklace import _ffc
        sage: list(_ffc([3,3])) #necklaces
        [[0, 1, 0, 1, 0, 1],
         [0, 0, 1, 1, 0, 1],
         [0, 0, 1, 0, 1, 1],
         [0, 0, 0, 1, 1, 1]]
        sage: list(_ffc([3,3], equality=True)) #Lyndon words
        [[0, 0, 1, 1, 0, 1], [0, 0, 1, 0, 1, 1], [0, 0, 0, 1, 1, 1]]
    """
    e = list(content)
    a = [len(e)-1]*sum(e)
    r = [0] * sum(e)
    a[0] = 0
    e[0] -= 1
    k = len(e)

    rng_k = range(k)
    rng_k.reverse()
    dll = DoublyLinkedList(rng_k)
    if not e[0]: # == 0
        dll.hide(0)

    for x in _fast_fixed_content(a, e, 2, 1, k, r, 2, dll, equality=equality):
        yield x

def _fast_fixed_content(a, content, t, p, k, r, s, dll, equality=False):
    """
    EXAMPLES::

        sage: from sage.combinat.necklace import _fast_fixed_content
        sage: from sage.combinat.misc import DoublyLinkedList
        sage: e = [3,3]
        sage: a = [len(e)-1]*sum(e)
        sage: r = [0]*sum(e)
        sage: a[0] = 0
        sage: e[0] -= 1
        sage: k = len(e)
        sage: dll = DoublyLinkedList(list(reversed(range(k))))
        sage: if e[0] == 0: dll.hide(0)
        sage: list(_fast_fixed_content(a,e,2,1,k,r,2,dll))
        [[0, 1, 0, 1, 0, 1],
         [0, 0, 1, 1, 0, 1],
         [0, 0, 1, 0, 1, 1],
         [0, 0, 0, 1, 1, 1]]
        sage: list(_fast_fixed_content(a,e,2,1,k,r,2,dll,True))
        [[0, 0, 1, 1, 0, 1], [0, 0, 1, 0, 1, 1], [0, 0, 0, 1, 1, 1]]
    """
    n = len(a)
    if content[k-1] == n - t + 1:
        if content[k-1] == r[t-p-1]:
            if equality:
                if n == p:
                    yield a
            else:
                if not n % p: # == 0
                    yield a
        elif content[k-1] > r[t-p-1]:
            yield a
    elif content[0] != n-t+1:
        j = dll.head()
        sp = s
        while j != 'end' and j >= a[t-p-1]:
            #print s, j
            r[s-1] = t-s
            a[t-1] = j
            content[j] -= 1

            if not content[j]: # == 0
                dll.hide(j)

            if j != k-1:
                sp = t+1

            if j == a[t-p-1]:
                for x in _fast_fixed_content(a[:], content, t+1, p+0, k, r, sp, dll, equality=equality):
                    yield x
            else:
                for x in _fast_fixed_content(a[:], content, t+1, t+0, k, r, sp, dll, equality=equality):
                    yield x

            if not content[j]: # == 0
                dll.unhide(j)

            content[j] += 1
            j = dll.next(j)
        a[t-1] = k-1
    return


################################
# List Fixed Content Algorithm #
################################
def _lfc(content, equality=False):
    """
    EXAMPLES::

        sage: from sage.combinat.necklace import _lfc
        sage: list(_lfc([3,3])) #necklaces
        [[0, 1, 0, 1, 0, 1],
         [0, 0, 1, 1, 0, 1],
         [0, 0, 1, 0, 1, 1],
         [0, 0, 0, 1, 1, 1]]
        sage: list(_lfc([3,3], equality=True)) #Lyndon words
        [[0, 0, 1, 1, 0, 1], [0, 0, 1, 0, 1, 1], [0, 0, 0, 1, 1, 1]]
    """
    content = list(content)
    a = [0]*sum(content)
    content[0] -= 1
    k = len(content)

    rng_k = range(k)
    rng_k.reverse()
    dll = DoublyLinkedList(rng_k)

    if not content[0]: # == 0
        dll.hide(0)

    for z in _list_fixed_content(a, content, 2, 1, k, dll, equality=equality):
        yield z

def _list_fixed_content(a, content, t, p, k, dll, equality=False):
    """
    EXAMPLES::

        sage: from sage.combinat.necklace import _list_fixed_content
        sage: from sage.combinat.misc import DoublyLinkedList
        sage: e = [3,3]
        sage: a = [0]*sum(e)
        sage: e[0] -= 1
        sage: k = len(e)
        sage: dll = DoublyLinkedList(list(reversed(range(k))))
        sage: if e[0] == 0: dll.hide(0)
        sage: list(_list_fixed_content(a,e,2,1,k,dll))
        [[0, 1, 0, 1, 0, 1],
         [0, 0, 1, 1, 0, 1],
         [0, 0, 1, 0, 1, 1],
         [0, 0, 0, 1, 1, 1]]
        sage: list(_list_fixed_content(a,e,2,1,k,dll,True))
        [[0, 0, 1, 1, 0, 1], [0, 0, 1, 0, 1, 1], [0, 0, 0, 1, 1, 1]]
    """
    n = len(a)
    if t > n:
        if equality:
            if n == p:
                yield a
        else:
            if not n % p: # == 0
                yield a
    else:
        j = dll.head()
        while j != 'end' and j >= a[t-p-1]:
            a[t-1] = j
            content[j] -= 1

            if not content[j]: # == 0
                dll.hide(j)

            if j == a[t-p-1]:
                for z in _list_fixed_content(a[:], content[:], t+1, p+0, k, dll, equality=equality):
                    yield z
            else:
                for z in _list_fixed_content(a[:], content[:], t+1, t+0, k, dll, equality=equality):
                    yield z

            if not content[j]: # == 0
                dll.unhide(j)

            content[j] += 1
            j = dll.next(j)



################################
#Simple Fixed Content Algorithm#
################################
def _sfc(content, equality=False):
    """
    This wrapper function calls :meth:`sage.combinat.necklace._simple_fixed_content`.
    If ``equality`` is ``True`` the function returns Lyndon words with content
    vector equal to ``content``, otherwise it returns necklaces.

    INPUT:

    - ``content`` -- a list of non-negative integers with no leading 0s
    - ``equality`` -- boolean (optional, default: ``True``)

    .. WARNING::

        You will get incorrect results if there are leading 0's in ``content``.
        See :trac:`12997` and :trac:`17436`.

    EXAMPLES::

        sage: from sage.combinat.necklace import _sfc
        sage: list(_sfc([3,3])) #necklaces
        [[0, 0, 0, 1, 1, 1],
         [0, 0, 1, 0, 1, 1],
         [0, 0, 1, 1, 0, 1],
         [0, 1, 0, 1, 0, 1]]
        sage: list(_sfc([3,3], equality=True)) #Lyndon words
        [[0, 0, 0, 1, 1, 1], [0, 0, 1, 0, 1, 1], [0, 0, 1, 1, 0, 1]]
    """
    content = list(content)
    a = [0]*sum(content)
    content[0] -= 1
    k = len(content)
    return _simple_fixed_content(a, content, 2, 1, k, equality=equality)

def _simple_fixed_content(a, content, t, p, k, equality=False):
    """
    EXAMPLES::

        sage: from sage.combinat.necklace import _simple_fixed_content
        sage: content = [3,3]
        sage: a = [0]*sum(content)
        sage: content[0] -= 1
        sage: k = len(content); k
        2
        sage: list(_simple_fixed_content(a, content, 2, 1, k))
        [[0, 0, 0, 1, 1, 1],
         [0, 0, 1, 0, 1, 1],
         [0, 0, 1, 1, 0, 1],
         [0, 1, 0, 1, 0, 1]]
        sage: list(_simple_fixed_content(a, content, 2, 1, k, True))
        [[0, 0, 0, 1, 1, 1], [0, 0, 1, 0, 1, 1], [0, 0, 1, 1, 0, 1]]
    """
    n = len(a)
    if t > n:
        if equality:
            if n == p:
                yield a
        else:
            if not n % p: # == 0
                yield a
    else:
        r = range(a[t-p-1],k)
        for j in r:
            if content[j] > 0:
                a[t-1] = j
                content[j] -= 1
                if j == a[t-p-1]:
                    for z in _simple_fixed_content(a[:], content, t+1, p+0, k, equality=equality):
                        yield z
                else:
                    for z in _simple_fixed_content(a[:], content, t+1, t+0, k, equality=equality):
                        yield z
                content[j] += 1


def _lyn(w):
    """
    Returns the length of the longest prefix of ``w`` that is a Lyndon word.

    EXAMPLES::

        sage: import sage.combinat.necklace as necklace
        sage: necklace._lyn([0,1,1,0,0,1,2])
        3
        sage: necklace._lyn([0,0,0,1])
        4
        sage: necklace._lyn([2,1,0,0,2,2,1])
        1
    """
    p = 1
    k = max(w)+1
    for i in range(1, len(w)):
        b = w[i]
        a = w[:i]
        if b < a[i-p] or b > k-1:
            return p
        elif b == a[i-p]:
            pass
        else:
            p = i+1
    return p
