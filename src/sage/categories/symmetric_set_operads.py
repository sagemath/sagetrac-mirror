r"""
Symmetric Set Operads
"""
#*****************************************************************************
#  Copyright (C) 2011 Floent Hivert (CNRS) <Florent.Hivert@lri.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from category_types import Category
from sage.categories.all import SetOperads
from sage.misc.cachefunc import cached_method
from sage.categories.cartesian_product import cartesian_product
from sage.misc.abstract_method import abstract_method

class SymmetricSetOperads(Category):
    """
    The category of symmetric operads

    EXAMPLES::

      sage: SymmetricSetOperads()
      Category of symmetric set operads
      sage: SymmetricSetOperads().super_categories()
      [Category of set operads]

    TESTS::

        sage: C = SymmetricSetOperads()
        sage: TestSuite(C).run()
    """

    @cached_method
    def super_categories(self):
        """
        EXAMPLES::

            sage: SymmetricSetOperads().super_categories()
            [Category of set operads]
        """
        return [SetOperads()]

    class ParentMethods:

        @abstract_method(optional = True)
        def symmetric_group_action(self, elem, perm):
            """
            returns the action of ``perm`` over ``elem``
            """

    @cached_method
    def operad_orbits(self, n, generators):
        """
        Return the orbits in degree ``n``

        EXAMPLES::
        """
        if n == 2: return generators
        perms = Permutations(n)
        levn1 = list(operad_orbits(generators, n-1))
        total = set()
        res = dict()
        for j, x in enumerate(levn1):
            print "Progress: (%s/%s) #total=%s #orbits=%s"%(
                j,len(levn1), len(total), len(res))
            for g in generators:
                for i in range(1, n):
                    new = x.compose(g, i)
                    if new in total: continue
                    #orb = set(new.permute(p) for p in perms)
                    orb = orbit(new)
                    total.update(orb)
                    res[new] = list(orb)
        print "Done: (%s/%s) #total=%s  #orbits=%s"%(
            len(levn1), len(levn1), len(total), len(res))
        return res

    class ElementMethods:

        def permute(self, perm):
            """
            EXAMPLES::
            """
            return self.parent().symmetric_group_action(self, perm)

        def orbit(self):
            """
            Return the orbit of ``self`` under the group action

            EXAMPLES::
            """
            # importing this earlier leads to segfault !!
            from sage.groups.perm_gps.permgroup_named import SymmetricGroup
            n = self.degree()
            if n == 1:
                return [self]
            GR = SymmetricGroup(n)
            trans = GR([2,1]+range(3, n+1))
            self_permut = self.permute(trans)
            if n == 2:
                if self_permut == self:
                    return [self]
                else:
                    return [self, self_permut]
            cycle = GR([n]+range(1, n))
            done = set()
            todo = set([self, self_permut])
            while todo:
                new = todo.pop()
                done.add(new)
                for i in range(n):
                    new = new.permute(cycle)
                    if new not in done:
                        todo.add(new)
                    newp = new.permute(trans)
                    if newp not in done:
                        todo.add(newp)
            return list(done)


