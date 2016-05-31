r"""
Dual groups of Finite Multiplicative Abelian Groups

The basic idea is very simple. Let G be an abelian group and `G^*` its
dual (i.e., the group of homomorphisms from G to `\CC^\times`). Let
`g_j`, `j=1,..,n`, denote generators of `G` - say `g_j` is of order
`m_j>1`. There are generators `X_j`, `j=1,..,n`, of `G^*` for which
`X_j(g_j)=\exp(2\pi i/m_j)` and `X_i(g_j)=1` if `i\not= j`. These are
used to construct `G^*`.

Sage supports multiplicative abelian groups on any prescribed finite
number `n > 0` of generators. Use
:func:`~sage.groups.abelian_gps.abelian_group.AbelianGroup` function
to create an abelian group, the
:meth:`~sage.groups.abelian_gps.abelian_group.AbelianGroup_class.dual_group`
method to create its dual, and then the :meth:`gen` and :meth:`gens`
methods to obtain the corresponding generators. You can print the
generators as arbitrary strings using the optional ``names`` argument
to the
:meth:`~sage.groups.abelian_gps.abelian_group.AbelianGroup_class.dual_group`
method.

EXAMPLES::

    sage: F = AbelianGroup(5, [2,5,7,8,9], names='abcde')
    sage: (a, b, c, d, e) = F.gens()

    sage: Fd = F.dual_group(names='ABCDE')
    sage: Fd.base_ring()
    Cyclotomic Field of order 2520 and degree 576
    sage: A,B,C,D,E = Fd.gens()
    sage: A(a)
    -1
    sage: A(b), A(c), A(d), A(e)
    (1, 1, 1, 1)

    sage: Fd = F.dual_group(names='ABCDE', base_ring=CC)
    sage: A,B,C,D,E = Fd.gens()
    sage: A(a)    # abs tol 1e-8
    -1.00000000000000 + 0.00000000000000*I
    sage: A(b); A(c); A(d); A(e)
    1.00000000000000
    1.00000000000000
    1.00000000000000
    1.00000000000000

AUTHORS:

- David Joyner (2006-08) (based on abelian_groups)

- David Joyner (2006-10) modifications suggested by William Stein

- Volker Braun (2012-11) port to new Parent base. Use tuples for immutables.
  Default to cyclotomic base ring.
"""

###########################################################################
#  Copyright (C) 2006 William Stein <wstein@gmail.com>
#  Copyright (C) 2006 David Joyner  <wdjoyner@gmail.com>
#  Copyright (C) 2012 Volker Braun  <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
###########################################################################
from __future__ import print_function

from sage.rings.infinity import infinity
from sage.structure.category_object import normalize_names
from sage.structure.unique_representation import UniqueRepresentation
from sage.groups.abelian_gps.dual_abelian_group_element import (
    DualAbelianGroupElement, is_DualAbelianGroupElement )
from sage.misc.mrange import mrange
from sage.misc.cachefunc import cached_method
from sage.groups.group import AbelianGroup as AbelianGroupBase


def is_DualAbelianGroup(x):
    """
    Return True if `x` is the dual group of an abelian group.

    EXAMPLES::

        sage: from sage.groups.abelian_gps.dual_abelian_group import is_DualAbelianGroup
        sage: F = AbelianGroup(5,[3,5,7,8,9], names=list("abcde"))
        sage: Fd = F.dual_group()
        sage: is_DualAbelianGroup(Fd)
        True
        sage: F = AbelianGroup(3,[1,2,3], names='a')
        sage: Fd = F.dual_group()
        sage: Fd.gens()
        (1, X1, X2)
        sage: F.gens()
        (1, a1, a2)
    """
    return isinstance(x, DualAbelianGroup_class)


class DualAbelianGroup_class(UniqueRepresentation, AbelianGroupBase):
    """
    Dual of abelian group.

    EXAMPLES::

        sage: F = AbelianGroup(5,[3,5,7,8,9], names="abcde")
        sage: F.dual_group()
        Dual of Abelian Group isomorphic to Z/3Z x Z/5Z x Z/7Z x Z/8Z x Z/9Z
        over Cyclotomic Field of order 2520 and degree 576
        sage: F = AbelianGroup(4,[15,7,8,9], names="abcd")
        sage: F.dual_group(base_ring=CC)
        Dual of Abelian Group isomorphic to Z/15Z x Z/7Z x Z/8Z x Z/9Z
        over Complex Field with 53 bits of precision
    """
    Element = DualAbelianGroupElement

    def __init__(self, G, names, base_ring):
        """
        The Python constructor

        EXAMPLES::

            sage: F = AbelianGroup(5,[3,5,7,8,9], names="abcde")
            sage: F.dual_group()
            Dual of Abelian Group isomorphic to Z/3Z x Z/5Z x Z/7Z x Z/8Z x Z/9Z
            over Cyclotomic Field of order 2520 and degree 576
       """
        self._base_ring = base_ring
        self._group = G
        names = normalize_names(G.ngens(), names)
        self._assign_names(names)
        AbelianGroupBase.__init__(self) # TODO: category=CommutativeGroups()

    def group(self):
        """
        Return the group that ``self`` is the dual of.

        EXAMPLES::

            sage: F = AbelianGroup(3,[5,64,729], names=list("abc"))
            sage: Fd = F.dual_group(base_ring=CC)
            sage: Fd.group() is F
            True
        """
        return self._group

    def base_ring(self):
        """
        Return the scalars over which the group is dualized.

        EXAMPLES::

            sage: F = AbelianGroup(3,[5,64,729], names=list("abc"))
            sage: Fd = F.dual_group(base_ring=CC)
            sage: Fd.base_ring()
            Complex Field with 53 bits of precision
        """
        return self._base_ring

    def __str__(self):
        """
        Print method.

        EXAMPLES::

            sage: F = AbelianGroup(3,[5,64,729], names=list("abc"))
            sage: Fd = F.dual_group(base_ring=CC)
            sage: print(Fd)
            DualAbelianGroup( AbelianGroup ( 3, (5, 64, 729) ) )
        """
        s = "DualAbelianGroup( AbelianGroup ( %s, %s ) )"%(self.ngens(), self.gens_orders())
        return s

    def _repr_(self):
        """
        Return a string representation.

        EXAMPLES::

            sage: F = AbelianGroup(5, [2,5,7,8,9], names='abcde')
            sage: Fd = F.dual_group(names='ABCDE', base_ring=CyclotomicField(2*5*7*8*9))
            sage: Fd   # indirect doctest
            Dual of Abelian Group isomorphic to Z/2Z x Z/5Z x Z/7Z x Z/8Z x Z/9Z
            over Cyclotomic Field of order 5040 and degree 1152
            sage: Fd = F.dual_group(names='ABCDE', base_ring=CC)
            sage: Fd
            Dual of Abelian Group isomorphic to Z/2Z x Z/5Z x Z/7Z x Z/8Z x Z/9Z
            over Complex Field with 53 bits of precision
        """
        G = self.group()
        eldv = G.gens_orders()
        gp = ""
        for x in eldv:
            if x!=0:
                gp = gp + "Z/%sZ x "%x
            if x==0:
                gp = gp + "Z x "
        gp = gp[:-2].strip()
        s = 'Dual of Abelian Group isomorphic to ' + gp + ' over ' + str(self.base_ring())
        return s

    def _latex_(self):
        r"""
        Return latex representation of this group.

        EXAMPLES::

            sage: F = AbelianGroup(3, [2]*3)
            sage: Fd = F.dual_group()
            sage: Fd._latex_()
            '$\\mathrm{DualAbelianGroup}( AbelianGroup ( 3, (2, 2, 2) ) )$'
        """
        s = "$\mathrm{DualAbelianGroup}( AbelianGroup ( %s, %s ) )$"%(self.ngens(), self.gens_orders())
        return s

    def random_element(self):
        """
        Return a random element of this dual group.

        EXAMPLES::

            sage: G = AbelianGroup([2,3,9])
            sage: Gd = G.dual_group(base_ring=CC)
            sage: Gd.random_element()
            X1^2

            sage: N = 43^2-1
            sage: G = AbelianGroup([N],names="a")
            sage: Gd = G.dual_group(names="A", base_ring=CC)
            sage: a, = G.gens()
            sage: A, = Gd.gens()
            sage: x = a^(N/4); y = a^(N/3); z = a^(N/14)
            sage: X = A*Gd.random_element(); X
            A^615
            sage: len([a for a in [x,y,z] if abs(X(a)-1)>10^(-8)])
            2
        """
        from sage.misc.prandom import randint
        result = self.one()
        for g in self.gens():
            order = g.order()
            result *= g**(randint(0,order))
        return result

    def gen(self, i=0):
        """
        The `i`-th generator of the abelian group.

        EXAMPLES::

            sage: F = AbelianGroup(3,[1,2,3],names='a')
            sage: Fd = F.dual_group(names="A")
            sage: Fd.0
            1
            sage: Fd.1
            A1
            sage: Fd.gens_orders()
            (1, 2, 3)
        """
        n = self.group().ngens()
        if i < 0 or i >= n:
            raise IndexError("Argument i (= %s) must be between 0 and %s."%(i, n-1))
        x = [0]*n
        if self.gens_orders()[i] != 1:
            x[i] = 1
        return self.element_class(self, x)

    def gens(self):
        """
        Return the generators for the group.

        OUTPUT:

        A tuple of group elements generating the group.

        EXAMPLES::

            sage: F = AbelianGroup([7,11]).dual_group()
            sage: F.gens()
            (X0, X1)
        """
        n = self.group().ngens()
        return tuple(self.gen(i) for i in range(n))

    def ngens(self):
        """
        The number of generators of the dual group.

        EXAMPLES::

            sage: F = AbelianGroup([7]*100)
            sage: Fd = F.dual_group()
            sage: Fd.ngens()
            100
        """
        return self.group().ngens()

    def gens_orders(self):
        """
        The orders of the generators of the dual group.

        OUTPUT:

        A tuple of integers.

        EXAMPLES::

            sage: F = AbelianGroup([5]*1000)
            sage: Fd = F.dual_group()
            sage: invs = Fd.gens_orders(); len(invs)
            1000
        """
        return self.group().gens_orders()

    def invariants(self):
        """
        The invariants of the dual group.

        You should use :meth:`gens_orders` instead.

        EXAMPLES::

            sage: F = AbelianGroup([5]*1000)
            sage: Fd = F.dual_group()
            sage: invs = Fd.gens_orders(); len(invs)
            1000
        """
        # TODO: deprecate
        return self.group().gens_orders()

    def __contains__(self,X):
        """
        Implements "in".

        EXAMPLES::

            sage: F = AbelianGroup(5,[2, 3, 5, 7, 8], names="abcde")
            sage: a,b,c,d,e = F.gens()
            sage: Fd = F.dual_group(names = "ABCDE")
            sage: A,B,C,D,E = Fd.gens()
            sage: A*B^2*D^7 in Fd
            True
        """
        return X.parent() == self and is_DualAbelianGroupElement(X)

    def order(self):
        """
        Return the order of this group.

        EXAMPLES::

            sage: G = AbelianGroup([2,3,9])
            sage: Gd = G.dual_group()
            sage: Gd.order()
            54
        """
        G = self.group()
        return G.order()

    def is_commutative(self):
        """
        Return True since this group is commutative.

        EXAMPLES::

            sage: G = AbelianGroup([2,3,9])
            sage: Gd = G.dual_group()
            sage: Gd.is_commutative()
            True
            sage: Gd.is_abelian()
            True
        """
        return True

    @cached_method
    def list(self):
        """
        Return tuple of all elements of this group.

        EXAMPLES::

            sage: G = AbelianGroup([2,3], names="ab")
            sage: Gd = G.dual_group(names="AB")
            sage: Gd.list()
            (1, B, B^2, A, A*B, A*B^2)
        """
        if not(self.is_finite()):
           raise NotImplementedError("Group must be finite")
        invs = self.gens_orders()
        T = mrange(invs)
        n = self.order()
        L = tuple( self(t) for t in T )
        return L

    def __iter__(self):
        """
        Return an iterator over the elements of this group.

        EXAMPLES::

            sage: G = AbelianGroup([2,3], names="ab")
            sage: Gd = G.dual_group(names="AB")
            sage: [X for X in Gd]
            [1, B, B^2, A, A*B, A*B^2]
            sage: N = 43^2-1
            sage: G = AbelianGroup([N],names="a")
            sage: Gd = G.dual_group(names="A", base_ring=CC)
            sage: a, = G.gens()
            sage: A, = Gd.gens()
            sage: x = a^(N/4)
            sage: y = a^(N/3)
            sage: z = a^(N/14)
            sage: len([X for X in Gd if abs(X(x)-1)>0.01 and abs(X(y)-1)>0.01 and abs(X(z)-1)>0.01])
            880
        """
        for g in self.list():
            yield g



