
"""
Bracelets

The algorithm used in this file comes from

- Sawada, Joe.  "Generating Bracelets in constant amortized time", SIAM
  Journal on Computing archive Volume 31 , Issue 1 (2001)
"""
#*****************************************************************************
#       Copyright (C) 2011 Daniel Recoskie <danielrecoskie@gmail.com>,
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from sage.structure.parent import Parent
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.rings.arith import euler_phi, divisors
from sage.rings.integer import Integer


def Bracelets(n, k):
    """
    Returns the combinatorial class of bracelets.

    EXAMPLES::

        sage: Bracelets(3, 2)
        Bracelets with n = 3 and k = 2
        sage: Bracelets(3, 2).cardinality()
        4
        sage: Bracelets(3, 2).first()
        [0, 0, 0]
        sage: Bracelets(3, 2).last()
        [1, 1, 1]
        sage: Bracelets(3, 2).list()
        [[0, 0, 0], [0, 0, 1], [0, 1, 1], [1, 1, 1]]
    """
    return Bracelets_NK(n, k)


class Bracelets_NK(Parent):
    def __init__(self, n, k):
        """
        TESTS::

            sage: B = Bracelets(4, 2)
            sage: B == loads(dumps(B))
            True
        """
        super(Bracelets_NK, self).__init__(self,
                                           category=FiniteEnumeratedSets())

        self.n = n
        self.k = k
        self.a = [0] * (n + 1)

    def __repr__(self):
        """
        TESTS::

            sage: repr(Bracelets(2,1))
            'Bracelets with n = 2 and k = 1'
        """
        return 'Bracelets with n = {0} and k = {1}'.format(self.n, self.k)

    def __cmp__(self, other):
        """
        """
        return cmp(repr(self), repr(other))

    def __contains__(self, x):
        """
        EXAMPLES::

            sage: [0,0,2,2,1] in Bracelets(5,3)
            False
            sage: [1,1,2,1,2] in Bracelets(5,3)
            True
            sage: all([ n in Bracelets(5,3) for n in Bracelets(5,3)])
            True

        TESTS::

            sage: [0,0,2,2,1] in Bracelets(5,3)
            False
            sage: [1,1,2,1,2] in Bracelets(5,3)
            True
            sage: all([ n in Bracelets(5,3) for n in Bracelets(5,3)])
            True
        """
        l = list(x)

        if len(l) != self.n:
            return False

        #Check to make sure l is a list of valid integers
        if not all([isinstance(i, (int, Integer))
                    and i >= 0 and i < self.k for i in l]):
            return False

        #Check to make sure that l is lexicographically less
        #than all of its cyclic shifts
        cyclic_shift = l[:]
        for i in range(self.n - 1):
            cyclic_shift = cyclic_shift[1:] + cyclic_shift[:1]
            if cyclic_shift < l:
                return False

        #Also check its reverse
        cyclic_shift = l[:]
        cyclic_shift.reverse()
        for i in range(self.n - 1):
            cyclic_shift = cyclic_shift[1:] + cyclic_shift[:1]
            if cyclic_shift < l:
                return False

        return True

    def cardinality(self):
        """
        Returns the number of bracelets with given n and k.

        TESTS::

            sage: Bracelets(7, 4).cardinality() == len(Bracelets(7, 4).list())
            True
        """
        nk = sum([euler_phi(d) * self.k ** (self.n / d)
                  for d in divisors(self.n)]) / self.n

        if self.n % 2 == 0:
            return (nk + (self.k + 1) / 2 * self.k ** (self.n / 2)) / 2
        else:
            return (nk + self.k ** ((self.n + 1) / 2)) / 2

    def __iter__(self):
        """
        An iterator for bracelets.

        """
        for i in self._gen_b(1, 1, 1, -1, 0, False):
            yield i

    def _gen_b(self, t, p, r, u, v, rs):
        """
        Algorithm for generating bracelets as described in Sawada,
            Joe.  "Generating Bracelets in constant amortized time",
            SIAM Journal on Computing archive Volume 31 , Issue 1
            (2001)
        """
        if t - 1 > (self.n - r)/2 + r:
            if self.a[t - 1] > self.a[self.n - t + 2 + r]:
                rs = False
            elif self.a[t - 1] < self.a[self.n - t + 2 + r]:
                rs = True

        if t > self.n:
            if rs == False and self.n % p == 0:
                yield self.a[1:]
        else:
            self.a[t] = self.a[t - p]
            if self.a[t] == self.a[1]:
                v += 1
            else:
                v = 0

            if u == -1 and self.a[t-1] != self.a[1]:
                u = t-2
                r = t-2

            if u != -1 and t == self.n and self.a[self.n] == self.a[1]:
                pass
            elif u == v:
                rev = self.check_rev(t, u)
                if rev == 0:
                    for perm in self._gen_b(t+1, p, r, u, v, rs):
                        yield perm
                if rev == 1:
                    for perm in self._gen_b(t+1, p, t, u, v, False):
                        yield perm
            else:
                for perm in self._gen_b(t+1, p, r, u, v, rs):
                    yield perm

            for j in range(self.a[t-p]+1, self.k):
                self.a[t] = j
                for perm in self._gen_b(t+1, t, r, u, 0, rs):
                    yield perm

    def check_rev(self, t, i):
        """
        """
        for j in range(i+1, (t+1)/2 + 1):
            if self.a[j] < self.a[t-j+1]:
                return 0
            elif self.a[j] > self.a[t-j+1]:
                return -1

        return 1
