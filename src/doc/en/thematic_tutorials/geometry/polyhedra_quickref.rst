.. -*- coding: utf-8 -*-
.. linkall

.. _polyhedra_quickref:

=====================================
Quick reference for polyhedra in Sage
=====================================

.. MODULEAUTHOR:: Jean-Philippe Labbé <labbe@math.fu-berlin.de>
                  Vincent Delecroix <vincent.delecroix@u-bordeaux.fr>

List of Polyhedron methods
==========================

**H and V-representation**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.base_ring` | ring on which the polyhedron is defined
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.ambient_space` | ambient vector space or free module
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.Hrepresentation_space` | vector space or free module used for the vectors of the H-representation
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.Vrepresentation_space` | vector space or free module used for the vectors of the V-representation
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.n_Hrepresentation` | number of elements in the H-representation (sum of the number of equations and inequalities)
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.n_Vrepresentation` | number of elements in the V-representation (sum of vertices, rays and lines)
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.n_equations` | number of equations
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.n_inequalities` | number of inequalities
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.n_vertices` | number of vertices
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.n_rays` | number of rays
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.n_lines` | number of lines
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.n_facets` | number of facets

**Polyhedron boolean properties:**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.is_empty` | tests emptyness
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.is_universe` | tests whether a polyhedra is the whole ambient space
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.is_full_dimensional` | tests if the polyhedron has the same dimension as the ambient space
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.is_combinatorially_isomorphic` | tests whether two polyhedra are combinatorially isomorphic
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.is_compact` | tests compactness, or boundedness of a polyhedron
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.is_lattice_polytope` | tests whether a polyhedron is a lattice polytope
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.is_inscribed` | tests whether the polyhedron is inscribed in a sphere
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.is_minkowski_summand` | tests if the polyhedron can be used to produce another given polyhedron using a Minkowski sum.
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.is_neighborly` | tests whether the polyhedron has full skeleton until half of the dimension (or up to a certain dimension)
    :meth:`~sage.geometry.polyhedron.base_ZZ.Polyhedron_ZZ.is_reflexive` | tests if the polar of a lattice polytope is also a lattice polytope (only for :class:`Polyhedron over ZZ <sage.geometry.polyhedron.base_ZZ.Polyhedron_ZZ>`)
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.is_simple` |  checks whether the degree of all vertices is equal to the dimension of the polytope
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.is_simplex` | test whether a polytope is a simplex
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.is_simplicial` |  checks whether all faces of the polyhedron are simplices
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.is_lawrence_polytope` |  tests whether self is a Lawrence polytope
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.is_self_dual` |  tests whether the polytope is self-dual
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.is_pyramid` | test whether the polytope is a pyramid over one of its facets
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.is_bipyramid` | test whether the polytope is combinatorially equivalent to a bipyramid over some polytope
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.is_prism` | test whether the polytope is combinatorially equivalent to a prism of some polytope

**Enumerative properties**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.ambient_dim` |  the dimension of the ambient vector space
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.dim` |  the dimension of the polytope
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.dimension` |  alias of dim
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.f_vector` |  the `f`-vector (number of faces of each dimension)
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.neighborliness` | highest cardinality for which all `k`-subsets of the vertices are faces of the polyhedron

**Implementation properties**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.backend` | gives the backend used
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.base_ring` | gives the base ring used
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.change_ring` | changes the base ring

**Transforming polyhedra**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.minkowski_sum` | Minkowski sum of two polyhedra
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.minkowski_difference` | Minkowski difference of two polyhedra
    :meth:`~sage.geometry.polyhedron.base_ZZ.Polyhedron_ZZ.minkowski_decompositions` | Minkowski decomposition (only for :class:`Polyhedron over ZZ <sage.geometry.polyhedron.base_ZZ.Polyhedron_ZZ>`)
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.product` | cartesian product of two polyhedra
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.intersection` | intersection of two polyhedra
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.join` | join of two polyhedra
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.convex_hull` | convex hull of the union of two polyhedra
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.affine_hull` | constructs an affinely equivalent full dimensional polyhedra
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.barycentric_subdivision` | constructs a geometric realization of the barycentric subdivision
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.dilation` |  scalar dilation
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.face_truncation` | truncates a specific face
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.face_split` | returns the face splitting of a face of self
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.one_point_suspension` | the one-point suspension over a vertex of self (face splitting of a vertex)
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.stack` | stack a face of the polyhedron
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.lattice_polytope` | returns an encompassing lattice polytope.
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.polar` | returns the polar of a polytope (needs to be compact)
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.prism` | prism over a polyhedron (increases both the dimension of the polyhedron and the dimension of the ambient space)
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.pyramid` | pyramid over a polyhedron (increases both the dimension of the polyhedron and the dimension of the ambient space)
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.bipyramid` | bipyramid over a polyhedron (increases both the dimension of the polyhedron and the dimension of the ambient)
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.translation` | translates by a given vector
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.truncation` | truncates all vertices simultaneously
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.lawrence_extension` | returns the Lawrence extension of self on a given point
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.lawrence_polytope` | returns the Lawrence polytope of self
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.wedge` | returns the wedge over a face of self

**Combinatorics**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.face_lattice` | the face lattice
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.combinatorial_automorphism_group` | the automorphism group of the underlying combinatorial polytope
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.graph`, :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.vertex_graph` | underlying graph
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.vertex_digraph` | digraph (orientation of edges determined by a linear form)
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.vertex_facet_graph` | bipartite digraph given vertex-facet adjacency
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.adjacency_matrix` | adjacency matrix
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.incidence_matrix` | incidence matrix
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.facet_adjacency_matrix` | adjacency matrix of the facets
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.vertex_adjacency_matrix` | adjacency matrix of the vertices

**Integral points**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`~sage.geometry.polyhedron.base_ZZ.Polyhedron_ZZ.ehrhart_polynomial` | the Ehrhart polynomial (only for :class:`Polyhedron over ZZ <sage.geometry.polyhedron.base_ZZ.Polyhedron_ZZ>`)
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.integral_points` | list of integral points
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.integral_points_count` | number of integral points
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.get_integral_point` | get the i-th integral point without computing all interior lattice points
    :meth:`~sage.geometry.polyhedron.base_ZZ.Polyhedron_ZZ.has_IP_property` | checks whether the origin is an interior lattice point and compactness (only for :class:`Polyhedron over ZZ <sage.geometry.polyhedron.base_ZZ.Polyhedron_ZZ>`)
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.random_integral_point` | get a random integral point


**Getting related geometric objects**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.center` | returns the average of the vertices of the polyhedron
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.representative_point` | returns the sum of the center and the rays
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.face_fan` | returns the fan spanned by the faces of the polyhedron
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.faces` | the list of faces
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.facets` | the list of facets
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.normal_fan` | returns the fan spanned by the normals of the supporting hyperplanes of the polyhedron
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.gale_transform` | returns the (affine) Gale transform of the vertices of the polyhedron
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.hyperplane_arrangement` | returns the hyperplane arrangement given by the defining facets of the polyhedron
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.to_linear_program` | transform the polyhedra into a Linear Program
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.triangulate` | returns a triangulation of the polyhedron
    :meth:`~sage.geometry.polyhedron.base_ZZ.Polyhedron_ZZ.fibration_generator` | returns an iterator of the fibrations of the lattice polytope (only for :class:`Polyhedron over ZZ <sage.geometry.polyhedron.base_ZZ.Polyhedron_ZZ>`)

**Other**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |


    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.bounded_edges` | generator for bounded edges
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.bounding_box` | returns the vertices of an encompassing cube
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.contains` | tests whether the polyhedron contains a vector
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.interior_contains` | tests whether the polyhedron contains a vector in its interior using the ambient topology
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.relative_interior_contains` | tests whether the polyhedron contains a vector in its relative interior
    :meth:`~sage.geometry.polyhedron.base_ZZ.Polyhedron_ZZ.find_translation` | returns the translation vector between two translation of two polyhedron (only for :class:`Polyhedron over ZZ <sage.geometry.polyhedron.base_ZZ.Polyhedron_ZZ>`)
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.integrate` | computes the integral of a polynomial over the polyhedron
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.radius` | returns the radius of the smallest sphere containing the polyhedron
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.radius_square` | returns the square of the radius of the smallest sphere containing the polyhedron
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.volume` | computes different volumes of the polyhedron
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.restricted_automorphism_group` | returns the restricted automorphism group
    :meth:`~sage.geometry.polyhedron.ppl_lattice_polytope.LatticePolytope_PPL_class.lattice_automorphism_group` | returns the lattice automorphism group. Only for :class:`PPL Lattice Polytope <sage.geometry.polyhedron.ppl_lattice_polytope.LatticePolytope_PPL_class>`

