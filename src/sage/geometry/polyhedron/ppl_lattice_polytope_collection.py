"""
Collections of fast Lattice Polytopes using PPL.

The :func:`LatticePolytopeCollection` creates a collection of
PPL LatticePolytopes.

AUTHORS:

    - Jan Keitel: initial version, 2013
"""

########################################################################
#       Copyright (C) 2013 Jan Keitel <jan.keitel@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
########################################################################

from sage.structure.all import SageObject
from sage.geometry.polyhedron.ppl_lattice_polytope import (
    LatticePolytope_PPL_class )

def LatticePolytopeCollection(*args):
    r"""
    Return a LatticePolytopeCollection consisting of the lattice polytopes
    passed as arguments.

    INPUT:

    OUTPUT:

    A :class:`LatticePolytopeCollection<LatticePolytopeCollection_class>`.

    EXAMPLES::
    """
    def ppl_instance_check(x):
        return all(isinstance(i, LatticePolytope_PPL_class) for i in x)  
    if len(args) == 1:
        if isinstance(args[0], (list, tuple)) and \
            ppl_instance_check(args[0]):
            return LatticePolytopeCollection_class(*args[0])
    if ppl_instance_check(args):
        return LatticePolytopeCollection_class(*args)
    raise ValueError('Invalid parameters. Only LatticePolytope_PPL or ' +
                     'iterables thereof are allowed.')

class LatticePolytopeCollection_class(SageObject):

    def __init__(self, *args):
        r"""
        Construct an element of class:`LatticePolytopeCollection_class`.

        Do not call this constructor directly, but use 
        func:`LatticePolytopeCollection` instead.

        INPUT:

        See func:`LatticePolytopeCollection` for more information.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.ppl_lattice_polytope_collection \
            ....: import LatticePolytopeCollection_class
            sage: from sage.geometry.polyhedron.ppl_lattice_polytope import \
            ....: LatticePolytope_PPL
            sage: p = LatticePolytope_PPL((1, 0), (0, 1))
            sage: q = LatticePolytope_PPL((1, 0), (0, -1))
            sage: LatticePolytopeCollection_class(p, q)
            A collection of 2 lattice polytopes.
        """
        self._polytopes = tuple(args)
        super(LatticePolytopeCollection_class, self).__init__()

    def _eq_(self, other):
        if type(self) <> type(other):
            return False
        if self.n_polytopes() <> other.n_polytopes():
            return False
        return all(i == j for i, j in zip(self.polytopes(), other.polytopes()))

    def _repr_(self):
        r"""
        Return the string representation of ``self``.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.ppl_lattice_polytope_collection \
            ....: import LatticePolytopeCollection
            sage: from sage.geometry.polyhedron.ppl_lattice_polytope import \
            ....: LatticePolytope_PPL
            sage: p = LatticePolytope_PPL([(1, 0), (0, 1), (-1, -1)])
            sage: q = LatticePolytope_PPL([(1, 0, 0), (0, 1, 0),\
            ....: (0, 2, 3), (-1, -1, 2)])
            sage: LatticePolytopeCollection(p)
            A collection of 1 lattice polytope.
            sage: LatticePolytopeCollection(p, q)
            A collection of 2 lattice polytopes.
        """
        n = self.n_polytopes()
        if n > 1:
            end = 's'
        else:
            end = ''
        out = 'A collection of {0} lattice polytope{1}.'.format(n, end)
        return out

    def _embedding_table(self, other):
        r"""
        Compute the Boolean matrix corresponding to the possible embeddings
        of the polytopes of ``self`` into the polytopes of ``other``.

        INPUT:

        - ``other`` -- a 
            :class:`LatticePolytopeCollection<LatticePolytopeCollection_class>`.

        OUTPUT:

        A tuple of tuples of Booleans.

        EXAMPLES::
        """
        res = []
        for p in self.polytopes():
            tmp = []
            for q in other.polytopes():
                tmp.append(p.is_embeddable_in(q))
            res.append(tuple(tmp))
        return tuple(res)

    def embed_in_strong(self, other):
        r"""
        Embed ``self`` in the LatticePolytopeCollection ``other`` using
        only one linear map.

        INPUT:

        - ``other`` -- a
            :class:`LatticePolytopeCollection<LatticePolytopeCollection_class>`.

        OUTPUT:

        EXAMPLES::
        """
        # The existence of a weak embedding is a necessary condition, so begin looking
        # at embedding tables.
        # Begin by looking for all sub-polytopes of other with the right
        # numbers of vertices and integral points
        from sage.combinat.all import Arrangements
        from sage.geometry.polyhedron.ppl_lattice_polytope import \
            _compute_affine_normal_forms
        from sage.geometry.polyhedron.lattice_euclidean_group_element import \
            LatticePolytopesNotIsomorphicError
        n_self = self.n_polytopes()
        # Order other's polytopes by their number of integral points
        polys_enumerated = list(enumerate(other.polytopes()))
        polys_sorted = sorted(polys_enumerated,
                              key=lambda x: x[1].n_integral_points())
        print 'Sorted polys:', polys_sorted
        sorted_indices = tuple([i[0] for i in polys_sorted])
        maps = Arrangements(range(other.n_polytopes()), n_self)
        if not maps:
            raise ValueError('The embedding should not have fewer polytopes.')
        embedding_table = self._embedding_table(other)
        not_working = True
        for map in maps:
            if not all(embedding_table[i][map[i]] for i in range(n_self)):
                continue
            # We have found a weak embedding.
            print 'Checking weak embedding', map
            # Find the codomain with the least integral points and the index
            # of the polytope to be embedded into it
            index_codomain = [i for i in sorted_indices if i in map][0]
            index_domain = map.index(index_codomain)
            # Find all relevant sub-polytopes in the relevant codomain
            domain = self.polytopes()[index_domain]
            codomain = other.polytopes()[index_codomain]
            n_pts = domain.n_integral_points()
            n_vertices = domain.n_vertices()
            dim = domain.affine_dimension()
            sub_polys = list(codomain.sub_polytope_generator_all(
                            n_pts, n_vertices, dim))
            if not sub_polys:
                print 'No relevant sub-polytopes!'
                continue
            _compute_affine_normal_forms(sub_polys)
            # We must improve the isomorphisms, this is much too slow
            print 'There are', len(sub_polys), 'sub-polytopes to check.'
            for sub in sub_polys:
                try:
                    invhom = domain.find_isomorphism(sub)
                except LatticePolytopesNotIsomorphicError:
                    print 'Not isomorphic.'
                    not_working = True
                    continue
                # Make sure that the other polytopes map correctly, too
                not_working = False
                for i in range(n_self):
                    if i == index_domain:
                        continue
                    do = self.polytopes()[i]
                    co = other.polytopes()[map[i]]
                    if not co.contains(invhom(do)):
                        print 'Isomorphism does not work for', do
                        not_working = True
                        break
                if not_working:
                    continue
            if not_working:
                continue

            '''codomains = [other.polytopes()[map[i]] for i in range(n_self)]
            # Pick the polytope with the least integral points


            embeddings = tuple([p.embed_in(q, **kwds) for p, q in
                          zip(self.polytopes(), codomains)])'''
            yield (tuple(map), invhom)

    def embed_in(self, other, **kwds):
        r"""
        Embed ``self`` in the LatticePolytopeCollection ``other``.

        By an embedding of a
        :class:`LatticePolytopeCollection<LatticePolytopeCollection_class>`
        `A` into another collection `B`, we mean that for every 
        `Q \in A` there exists a polytope `P_Q \in B` with an embedding
        `Q \hookrightarrow P_Q`. Furthermore, we require that
        `P_Q = P_{Q'}` if and only if `Q = Q'`.

        INPUT:

        - ``other`` -- a
            :class:`LatticePolytopeCollection<LatticePolytopeCollection_class>`.

        - for further optional parameters see 
            :meth:`~sage.geometry.lattice_polytope.LatticePolytopePPL_class.embed_in`.

        OUTPUT:

        A tuple consisting of two further tuples, namely

            1. A tuple `(i_1, i_2, \dots, i_n)` indicating that the `j`th
               polytope of ``self`` is to be embedded in the `i_j`th polytope
               of ``other``.

            2. A tuple of the return values from the ``embed_into`` method of
               `Q_j` into ``P_{i_j}``. See
               :meth:`~sage.geometry.lattice_polytope.LatticePolytopePPL_class.embed_in`
               for more information.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.ppl_lattice_polytope_collection \
            ....: import LatticePolytopeCollection
            sage: from sage.geometry.polyhedron.ppl_lattice_polytope import \
            ....: ReflexivePolytope_PPL
            sage: a, b, c = [ReflexivePolytope_PPL(2, i) for i in [0, 3, 15]]
            sage: L1 = LatticePolytopeCollection(a, b)
            sage: L2 = LatticePolytopeCollection(a, c)
            sage: for embedding in L1.embed_in(L2):
            ....:     print embedding
            ((0, 1), (The map A*x+b with A=
            [-1  1]
            [-1  0]
            b = 
            (0, 0), The map A*x+b with A=
            [ 1 -1]
            [ 1  0]
            b = 
            (0, 0)))
            sage: for embedding in L1.embed_in(L2, output='sub_polytope'):
            ....:     print embedding
            ((0, 1), (A 2-dimensional lattice polytope in ZZ^2 with 3 vertices, A
            2-dimensional lattice polytope in ZZ^2 with 4 vertices))
        """
        from sage.combinat.all import Arrangements
        n_self = self.n_polytopes()
        maps = Arrangements(range(other.n_polytopes()), n_self)
        if not maps:
            raise ValueError('The embedding should not have fewer polytopes.')

        embedding_table = self._embedding_table(other)
        # Check whether we can find an embedding
        for map in maps:
            if not all(embedding_table[i][map[i]] for i in range(n_self)):
                continue
            # We have found one. Bring it into desired form:
            codomains = [other.polytopes()[map[i]] for i in range(n_self)]
            embeddings = tuple([p.embed_in(q, **kwds) for p, q in
                          zip(self.polytopes(), codomains)])
            result = tuple([tuple(map), embeddings])
            yield result

    def embed_in_nef_partition(self, other, **kwds):
        r"""
        Embed ``self`` in a LatticePolytopeCollection corresponding to a nef partition.

        TODO: How can we do this effectively?

        INPUT:

        - ``other`` -- `other`` -- a
            :class:`LatticePolytopeCollection<LatticePolytopeCollection_class>`.

        OUTPUT:

        EXAMPLES::
        """
        pass

    def is_embeddable_in(self, other):
        r"""
        Check whether ``self`` can be embedded in ``other``.

        For more information on what is meant by and embedding of one
        :class:`LatticePolytopeCollection<LatticePolytopeCollection_class>`
        into another, see :meth:`embed_in`.

        INPUT:

        - ``other`` -- a
            :class:`LatticePolytopeCollection<LatticePolytopeCollection_class>`.

        OUTPUT:

        A Boolean.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.ppl_lattice_polytope import \
            ....: ReflexivePolytope_PPL
            sage: from sage.geometry.polyhedron.ppl_lattice_polytope_collection \
            ....: import LatticePolytopeCollection
            sage: a, b, c = [ReflexivePolytope_PPL(2, i) for i in [0, 3, 15]]
            sage: L1 = LatticePolytopeCollection(a, b)
            sage: L2 = LatticePolytopeCollection(a, c)
            sage: b.is_embeddable_in(c)
            True
            sage: [c.is_embeddable_in(i) for i in [a, b]]
            [False, False]
            sage: L1.is_embeddable_in(L2)
            True
            sage: L2.is_embeddable_in(L1) 
            False
        """
        embeddings = self.embed_in(other, output='sub_polytope')
        # Check that there is at least one embedding
        try:
            embeddings.next()
            return True
        except StopIteration:
            return False

    def is_isomorphic(self, other):
        r"""
        Check whether ``self`` is isomorphic to ``other``.

        Two 
        :class:`LatticePolytopeCollections<LatticePolytopeCollection_class>`
        are isomorphic if and only if the ordered lists of the affine normal
        forms of their polytopes agree.

        INPUT:

        - ``other`` -- a
            :class:`LatticePolytopeCollection<LatticePolytopeCollection_class>`.

        OUTPUT:

        A Boolean.

        EXAMPLES::
        """
        # A simple (but possibly very slow) way of sorting them is via their affine
        # normal forms
        if self.n_polytopes() <> other.n_polytopes():
            return False
        anfs_self = sorted([i.affine_normal_form() for i in self.polytopes()])
        anfs_other = sorted([i.affine_normal_form() for i in other.polytopes()])
        return anfs_self == anfs_other

    def n_polytopes(self):
        r"""
        Return the number of polytopes ``self`` consists of.

        OUTPUT:

        An integer.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.ppl_lattice_polytope_collection \
            ....: import LatticePolytopeCollection
            sage: from sage.geometry.polyhedron.ppl_lattice_polytope import \
            ....: LatticePolytope_PPL
            sage: p = LatticePolytope_PPL([(1, 0), (0, 1), (-1, -1)])
            sage: q = LatticePolytope_PPL([(1, 0, 0), (0, 1, 0),\
            ....: (0, 2, 3), (-1, -1, 2)])
            sage: L1 = LatticePolytopeCollection(p)
            sage: L1.n_polytopes()
            1
            sage: L2 = LatticePolytopeCollection([p,q]*5)
            sage: L2.n_polytopes()
            10
        """
        return len(self._polytopes)

    def polytopes(self):
        r"""
        Return the polytopes constituting ``self``.

        OUTPUT:

        A tuple of LatticePolytopes.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.ppl_lattice_polytope_collection \
            ....: import LatticePolytopeCollection
            sage: from sage.geometry.polyhedron.ppl_lattice_polytope import \
            ....: LatticePolytope_PPL
            sage: p = LatticePolytope_PPL([(1, 0), (0, 1), (-1, -1)])
            sage: q = LatticePolytope_PPL([(1, 0), (0, 1), (-1, -1), (1, 1)])
            sage: L = LatticePolytopeCollection(p, q)
            sage: L.polytopes()
            (A 2-dimensional lattice polytope in ZZ^2 with 3 vertices, A
            2-dimensional lattice polytope in ZZ^2 with 4 vertices)
        """
        return self._polytopes


#######################################################################################
# Infrastructure for embeddings into nef-partitions of three-dimensional reflexive
# polyhedra.
#######################################################################################
