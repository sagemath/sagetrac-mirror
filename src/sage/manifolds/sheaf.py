r"""
Presheaves and Sheaves over Manifolds

AUTHORS:

- Michael Jung (2021): initial version

"""

#******************************************************************************
#       Copyright (C) 2021 Michael Jung <m.jung at vu.nl>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.categories.functor import Functor
from sage.structure.sage_object import SageObject
from sage.misc.abstract_method import abstract_method
from .manifold import TopologicalManifold
from sage.categories.morphism import SetMorphism
from sage.categories.homset import Hom
from sage.categories.sets_cat import Sets

class Presheaf(Functor):
    r"""

    """
    def __init__(self, manifold, section_set_constructor):
        r"""

        EXAMPLES::

            sage: from sage.manifolds.sheaf import Presheaf
            sage: M = Manifold(2, 'M')
            sage: F = Presheaf(M, lambda d: d.scalar_field_algebra()); F
            Presheaf of sections constructed by <function <lambda> at 0x7f4d179cbaf0>
             on the 2-dimensional differentiable manifold M
            sage: U = M.open_subset('U')
            sage: F(U)
            Algebra of differentiable scalar fields on the Open subset U of the
             2-dimensional differentiable manifold M

        """
        if not isinstance(manifold, TopologicalManifold):
            raise ValueError(f'{manifold} must be an instance of {TopologicalManifold}')
        self._manifold = manifold
        self._section_set_constructor = section_set_constructor
        # TODO: assuming the category of open subsets is given by OpenSubsets(M)
        super(self).__init__(OpenSubsets(manifold), Sets())

    def _repr_name_(self):
        r"""

        """
        return 'Presheaf'

    def _repr_(self):
        r"""
        Return the string representation of ``self``.

        """
        repr = f'{self._repr_name_()} of sections constructed '
        repr += f'by {self._section_set_constructor} on the {self._manifold}'
        return repr

    def _apply_functor(self, open_subset):
        r"""
        Return the set of sections of ``self`` w.r.t. ``open_subset``.

        """
        return self._section_set_constructor(open_subset)

    def section_set(self, open_subset):
        r"""
        Return the set of sections of ``self`` w.r.t. ``open_subset``.

        """
        return self(open_subset)

    def restriction_morphism(self, from_open_subset, to_open_subset):
        r"""
        Return the restriction morphism from ``from_open_subset`` to
        ``to_open_subset``.

        EXAMPLES::

            sage: from sage.manifolds.sheaf import Presheaf
            sage: M = Manifold(2, 'M')
            sage: F = Presheaf(M, lambda d: d.scalar_field_algebra())
            sage: U = M.open_subset('U')
            sage: F.restriction_morphism(M, U)
            Generic morphism:
              From: Algebra of differentiable scalar fields on the 2-dimensional differentiable manifold M
              To:   Algebra of differentiable scalar fields on the Open subset U of the 2-dimensional differentiable manifold M

        """
        if open_subset not in self.domain():
            raise TypeError(f"{open_subset} must be an open subset of {self._manifold}")
        # TODO: assuming that ``inclusion`` yields the inclusion morphism between open subsets
        inclusion = from_open_subset.inclusion(to_open_subset)
        return self(inclusion)

    def _apply_functor_to_morphism(self, inclusion):
        r"""
        Apply ``self`` to an inclusion morphism from one open subset to
        another.

        EXAMPLES::



        """
        return SetMorphism(Hom(inclusion.domain(), inclusion.codomain(), Sets()),
                           lambda x: x.restrict(self(inclusion.codomain())))

class Sheaf(Presheaf):
    r"""

    """
    def concatenation(self, *args):
        r"""
        Return the concatenation of a list of sections.

        """
        if isinstance(args[0], (list, tuple)):
            sections = list(args[0])
        else:
            sections = args
        if len(sections) < 2:
            raise ValueError('input must contain at least two sheaf sections')
        res = sections[0].concatenate(sections[1])
        for sec in sections[2:]:
            res = res.concatenate(sec)
        return res

    def _repr_name_(self):
        r"""

        """
        return 'Sheaf'

class PresheafSection(SageObject):
    r"""

    """
    def __init__(self, domain):
        r"""

        """
        self._domain = domain
        self._restrictions = {} # dict. of restrictions of self on subsets
                                # of self._domain, with the subsets as keys
        self._extensions_graph = {self._domain: self}
                    # dict. of known extensions of self on bigger domains,
                    # including self, with domains as keys. Its elements can be
                    # seen as incoming edges on a graph.
        self._restrictions_graph = {self._domain: self}
                    # dict. of known restrictions of self on smaller domains,
                    # including self, with domains as keys. Its elements can be
                    # seen as outgoing edges on a graph.

    def _del_restrictions(self):
        r"""
        Delete the restrictions defined on ``self``.

        TESTS::



        """
        self._restrictions.clear()
        self._extensions_graph = {self._domain: self}
        self._restrictions_graph = {self._domain: self}

    @abstract_method
    def set_restriction(self, rst):
        r"""
        Set ``rst`` to a restriction of ``self``.

        """

    def restrict(self, subdomain):
        r"""
        Return the restriction of ``self`` to ``subdomain``.

        """
        if subdomain == self._domain:
            return self
        if subdomain not in self._restrictions:
            if not subdomain.is_subset(self._domain):
                raise ValueError("the provided domain is not a subset of " +
                                 "the presheaf section's domain")
            # First one tries to get the restriction from a tighter domain:
            for dom, rst in self._restrictions.items():
                if subdomain.is_subset(dom) and subdomain in rst._restrictions:
                    res = rst._restrictions[subdomain]
                    self._restrictions[subdomain] = res
                    self._restrictions_graph[subdomain] = res
                    res._extensions_graph.update(self._extensions_graph)
                    for ext in self._extensions_graph.values():
                        ext._restrictions[subdomain] = res
                        ext._restrictions_graph[subdomain] = res
                    return self._restrictions[subdomain]

            for dom, rst in self._restrictions.items():
                if subdomain.is_subset(dom) and dom is not self._domain:
                    self._restrictions[subdomain] = rst.restrict(subdomain)
                    self._restrictions_graph[subdomain] = rst.restrict(subdomain)
                    return self._restrictions[subdomain]

            # Secondly one tries to get the restriction from one previously
            # defined on a larger domain:
            for dom, ext in self._extensions_graph.items():
                if subdomain in ext._restrictions:
                    res = ext._restrictions_graph[subdomain]
                    self._restrictions[subdomain] = res
                    self._restrictions_graph[subdomain] = res
                    res._extensions_graph.update(self._extensions_graph)
                    for ext in self._extensions_graph.values():
                        ext._restrictions[subdomain] = res
                        ext._restrictions_graph[subdomain] = res
                    return self._restrictions[subdomain]

        # if that fails, create restriction from scratch
        res = self._construct_restriction_(subdomain)
        # update restriction graph
        res._extensions_graph.update(self._extensions_graph)
        for dom, ext in self._extensions_graph.items():
            ext._restrictions[subdomain] = res
            ext._restrictions_graph[subdomain] = res

        for dom, rst in self._restrictions.items():
            if dom.is_subset(subdomain):
                if rst is not res:
                    res._restrictions.update(rst._restrictions)
                res._restrictions_graph.update(rst._restrictions_graph)
                rst._extensions_graph.update(res._extensions_graph)
        self._restrictions[subdomain] = res
        self._restrictions_graph[subdomain] = res
        res._extensions_graph.update(self._extensions_graph)
        return res

    @abstract_method
    def _construct_restriction_(self, subdomain):
        r"""
        Construct the restriction of ``self`` to ``subdomain``.

        """

class SheafSection(PresheafSection):
    r"""

    """
    def __eq__(self, other):
        r"""

        """
        if self is other:
            return True
        if self._domain != other._domain:
            return False
        # compare on open covers until hitting a subdomain on which the
        # sections are actually comparable (e.g. parallelizable subdomains for
        # tensor fields)
        return all(self.restrict(subdom) == other.restrict(subdom)
                   for subdom in self._domain.open_covers(trivial=False))

    def concatenate(self, other):
        r"""
        Concatenate ``self`` with ``other``.

        """
        if self is other:
            return self
        inter = self._domain.intersection(other._domain)
        if self.restrict(inter) != other.restrict(inter):
            raise ValueError(f'{self} and {other} must coincide on intersection')
        return self._construct_concatenation_(other)

    @abstract_method
    def _construct_concatenation_(self):
        r"""
        Construct the concatenation of ``self`` with ``other``.

        """
