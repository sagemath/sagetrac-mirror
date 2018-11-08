r"""
Subschemes of affine space
AUTHORS:

- David Kohel (2005): initial version.
- William Stein (2005): initial version.
- Ben Hutz (2013): affine subschemes
"""

#*****************************************************************************
# Copyright (C) 2005 William Stein <wstein@gmail.com>
# Copyright (C) 2013 Ben Hutz <bn4941@gmail.com>
#
# Distributed under the terms of the GNU General Public License (GPL)
# as published by the Free Software Foundation; either version 2 of
# the License, or (at your option) any later version.
# http://www.gnu.org/licenses/
#*****************************************************************************

from sage.categories.fields import Fields
from sage.interfaces.all import singular
from sage.schemes.generic.algebraic_scheme import AlgebraicScheme_subscheme


class AlgebraicScheme_subscheme_affine(AlgebraicScheme_subscheme):
    r"""
    Construct an algebraic subscheme of affine space.

    .. WARNING::

        You should not create objects of this class directly. The
        preferred method to construct such subschemes is to use
        :meth:`~sage.schemes.affine.affine_space.AffineSpace_generic.subscheme`
        method of :class:`affine space
        <sage.schemes.affine.affine_space.AffineSpace_generic>`.

    INPUT:

    - ``A`` -- ambient :class:`affine space
      <sage.schemes.affine.affine_space.AffineSpace_generic>`

    - ``polynomials`` -- single polynomial, ideal or iterable of
      defining polynomials.

    EXAMPLES::

        sage: A3.<x, y, z> = AffineSpace(QQ, 3)
        sage: A3.subscheme([x^2-y*z])
        Closed subscheme of Affine Space of dimension 3 over Rational Field defined by:
          x^2 - y*z

    TESTS::

        sage: from sage.schemes.affine.affine_subscheme import AlgebraicScheme_subscheme_affine
        sage: AlgebraicScheme_subscheme_affine(A3, [x^2-y*z])
        Closed subscheme of Affine Space of dimension 3 over Rational Field defined by:
          x^2 - y*z
    """
    def __init__(self, A, polynomials, embedding_center=None,
                 embedding_codomain=None, embedding_images=None):
        """
        EXAMPLES::

            sage: A.<x,y,z> = AffineSpace(QQ, 3)
            sage: A.subscheme([y^2-x*z-x*y])
            Closed subscheme of Affine Space of dimension 3 over Rational Field defined by:
              -x*y + y^2 - x*z
        """
        AlgebraicScheme_subscheme.__init__(self, A, polynomials)
        if embedding_images is not None:
            self._embedding_morphism = self.hom(embedding_images, embedding_codomain)
        elif A._ambient_projective_space is not None:
            self._embedding_morphism = self.projective_embedding \
                (A._default_embedding_index, A._ambient_projective_space)
        if embedding_center is not None:
            self._embedding_center = self.point(embedding_center)

    def _morphism(self, *args, **kwds):
        r"""
        A morphism between two schemes in your category, usually defined via
        polynomials. Your morphism class should derive from
        :class:`SchemeMorphism_polynomial`. These morphisms will usually be
        elements of the Hom-set
        :class:`~sage.schemes.generic.homset.SchemeHomset_generic`.

        EXAMPLES::

            sage: A3.<x,y,z> = AffineSpace(3, ZZ)
            sage: A3._morphism(A3.Hom(A3), [x,y,z])
            Scheme endomorphism of Affine Space of dimension 3 over Integer Ring
              Defn: Defined on coordinates by sending (x, y, z) to
                    (x, y, z)

        """
        return self.ambient_space()._morphism(*args, **kwds)

    def dimension(self):
        """
        Return the dimension of the affine algebraic subscheme.

        OUTPUT: Integer.

        EXAMPLES::

            sage: A.<x,y> = AffineSpace(2, QQ)
            sage: A.subscheme([]).dimension()
            2
            sage: A.subscheme([x]).dimension()
            1
            sage: A.subscheme([x^5]).dimension()
            1
            sage: A.subscheme([x^2 + y^2 - 1]).dimension()
            1
            sage: A.subscheme([x*(x-1), y*(y-1)]).dimension()
            0

        Something less obvious::

            sage: A.<x,y,z,w> = AffineSpace(4, QQ)
            sage: X = A.subscheme([x^2, x^2*y^2 + z^2, z^2 - w^2, 10*x^2 + w^2 - z^2])
            sage: X
            Closed subscheme of Affine Space of dimension 4 over Rational Field defined by:
              x^2,
              x^2*y^2 + z^2,
              z^2 - w^2,
              10*x^2 - z^2 + w^2
            sage: X.dimension()
            1
        """
        try:
            return self.__dimension
        except AttributeError:
            self.__dimension = self.defining_ideal().dimension()
            return self.__dimension

    def projective_embedding(self, i=None, PP=None):
        """
        Returns a morphism from this affine scheme into an ambient
        projective space of the same dimension.

        The codomain of this morphism is the projective closure of this affine scheme in ``PP``,
        if given, or otherwise in a new projective space that is constructed.

        INPUT:

        -  ``i`` -- integer (default: dimension of self = last
           coordinate) determines which projective embedding to compute. The
           embedding is that which has a 1 in the i-th coordinate, numbered
           from 0.

        -  ``PP`` -- (default: None) ambient projective space, i.e., ambient space
            of codomain of morphism; this is constructed if it is not given.

        EXAMPLES::

            sage: A.<x, y, z> = AffineSpace(3, ZZ)
            sage: S = A.subscheme([x*y-z])
            sage: S.projective_embedding()
            Scheme morphism:
              From: Closed subscheme of Affine Space of dimension 3 over Integer Ring defined by:
              x*y - z
              To:   Closed subscheme of Projective Space of dimension 3 over Integer Ring defined by:
              x0*x1 - x2*x3
              Defn: Defined on coordinates by sending (x, y, z) to
                    (x : y : z : 1)

        ::

            sage: A.<x, y, z> = AffineSpace(3, ZZ)
            sage: P = ProjectiveSpace(3,ZZ,'u')
            sage: S = A.subscheme([x^2-y*z])
            sage: S.projective_embedding(1,P)
            Scheme morphism:
              From: Closed subscheme of Affine Space of dimension 3 over Integer
            Ring defined by:
              x^2 - y*z
              To:   Closed subscheme of Projective Space of dimension 3 over Integer
            Ring defined by:
              u0^2 - u2*u3
              Defn: Defined on coordinates by sending (x, y, z) to
                    (x : 1 : y : z)

        ::

            sage: A.<x,y,z> = AffineSpace(QQ, 3)
            sage: X = A.subscheme([y - x^2, z - x^3])
            sage: X.projective_embedding()
            Scheme morphism:
              From: Closed subscheme of Affine Space of dimension 3 over Rational
            Field defined by:
              -x^2 + y,
              -x^3 + z
              To:   Closed subscheme of Projective Space of dimension 3 over
            Rational Field defined by:
              x0^2 - x1*x3,
              x0*x1 - x2*x3,
              x1^2 - x0*x2
              Defn: Defined on coordinates by sending (x, y, z) to
                    (x : y : z : 1)

        """
        AA = self.ambient_space()
        n = AA.dimension_relative()
        if i is None:
            try:
                return self._embedding_morphism
            except AttributeError:
                i = n
        i = int(i)
        if i < 0 or i > n:
            raise ValueError("Argument i (=%s) must be between 0 and %s, inclusive"%(i, n))
        try:
            phi = self.__projective_embedding[i]
            #assume that if you've passed in a new ambient projective space
            #you want to override the existing embedding
            if PP is None or phi.codomain().ambient_space() == PP:
                return(phi)
        except AttributeError:
            self.__projective_embedding = {}
        except KeyError:
            pass
        if PP is None:
            PP = AA.projective_embedding(i).codomain()
        elif PP.dimension_relative() != n:
            raise ValueError("Projective Space must be of dimension %s"%(n))
        PR = PP.coordinate_ring()
        # Groebner basis w.r.t. a graded monomial order computed here to ensure
        # after homogenization, the basis elements will generate the defining
        # ideal of the projective closure of this affine subscheme
        R = AA.coordinate_ring()
        G = self.defining_ideal().groebner_basis()
        v = list(PP.gens())
        z = v.pop(i)
        phi = R.hom(v,PR)
        v.append(z)
        X = PP.subscheme([phi(f).homogenize(i) for f in G])
        v = list(R.gens())
        v.insert(i, R(1))
        phi = self.hom(v, X)
        self.__projective_embedding[i] = phi
        return phi

    def projective_closure(self, i=None, PP=None):
        r"""
        Return the projective closure of this affine subscheme.

        INPUT:

        - ``i`` -- (default: None) determines the embedding to use to compute the projective
          closure of this affine subscheme. The embedding used is the one which has a 1 in the
          i-th coordinate, numbered from 0.

        -  ``PP`` -- (default: None) ambient projective space, i.e., ambient space
           of codomain of morphism; this is constructed if it is not given.

        OUTPUT:

        - a projective subscheme.

        EXAMPLES::

            sage: A.<x,y,z,w> = AffineSpace(QQ,4)
            sage: X = A.subscheme([x^2 - y, x*y - z, y^2 - w, x*z - w, y*z - x*w, z^2 - y*w])
            sage: X.projective_closure()
            Closed subscheme of Projective Space of dimension 4 over Rational Field
            defined by:
              x0^2 - x1*x4,
              x0*x1 - x2*x4,
              x1^2 - x3*x4,
              x0*x2 - x3*x4,
              x1*x2 - x0*x3,
              x2^2 - x1*x3

        ::

            sage: A.<x,y,z> = AffineSpace(QQ, 3)
            sage: P.<a,b,c,d> = ProjectiveSpace(QQ, 3)
            sage: X = A.subscheme([z - x^2 - y^2])
            sage: X.projective_closure(1, P).ambient_space() == P
            True
        """
        return self.projective_embedding(i, PP).codomain()

    def is_smooth(self, point=None):
        r"""
        Test whether the algebraic subscheme is smooth.

        INPUT:

        - ``point`` -- A point or ``None`` (default). The point to
          test smoothness at.

        OUTPUT:

        Boolean. If no point was specified, returns whether the
        algebraic subscheme is smooth everywhere. Otherwise,
        smoothness at the specified point is tested.

        EXAMPLES::

            sage: A2.<x,y> = AffineSpace(2,QQ)
            sage: cuspidal_curve = A2.subscheme([y^2-x^3])
            sage: cuspidal_curve
            Closed subscheme of Affine Space of dimension 2 over Rational Field defined by:
              -x^3 + y^2
            sage: smooth_point = cuspidal_curve.point([1,1])
            sage: smooth_point in cuspidal_curve
            True
            sage: singular_point = cuspidal_curve.point([0,0])
            sage: singular_point in cuspidal_curve
            True
            sage: cuspidal_curve.is_smooth(smooth_point)
            True
            sage: cuspidal_curve.is_smooth(singular_point)
            False
            sage: cuspidal_curve.is_smooth()
            False
        """
        R = self.ambient_space().coordinate_ring()
        if not point is None:
            self._check_satisfies_equations(point)
            point_subs = dict(zip(R.gens(), point))
            Jac = self.Jacobian().subs(point_subs)
            return not Jac.is_zero()

        # testing smoothness everywhere tends to be expensive
        try:
            return self._smooth
        except AttributeError:
            pass
        sing_dim = self.Jacobian().dimension()
        self._smooth = (sing_dim == -1)
        return self._smooth

    def intersection_multiplicity(self, X, P):
        r"""
        Return the intersection multiplicity of this subscheme and the subscheme ``X`` at the point ``P``.

        The intersection of this subscheme with ``X`` must be proper, that is `\mathrm{codim}(self\cap
        X) = \mathrm{codim}(self) + \mathrm{codim}(X)`, and must also be finite. We use Serre's Tor
        formula to compute the intersection multiplicity. If `I`, `J` are the defining ideals of ``self``, ``X``,
        respectively, then this is `\sum_{i=0}^{\infty}(-1)^i\mathrm{length}(\mathrm{Tor}_{\mathcal{O}_{A,p}}^{i}
        (\mathcal{O}_{A,p}/I,\mathcal{O}_{A,p}/J))` where `A` is the affine ambient space of these subschemes.

        INPUT:

        - ``X`` -- subscheme in the same ambient space as this subscheme.

        - ``P`` -- a point in the intersection of this subscheme with ``X``.

        OUTPUT: An integer.

        EXAMPLES::

            sage: A.<x,y> = AffineSpace(QQ, 2)
            sage: C = Curve([y^2 - x^3 - x^2], A)
            sage: D = Curve([y^2 + x^3], A)
            sage: Q = A([0,0])
            sage: C.intersection_multiplicity(D, Q)
            4

        ::

            sage: R.<a> = QQ[]
            sage: K.<b> = NumberField(a^6 - 3*a^5 + 5*a^4 - 5*a^3 + 5*a^2 - 3*a + 1)
            sage: A.<x,y,z,w> = AffineSpace(K, 4)
            sage: X = A.subscheme([x*y, y*z + 7, w^3 - x^3])
            sage: Y = A.subscheme([x - z^3 + z + 1])
            sage: Q = A([0, -7*b^5 + 21*b^4 - 28*b^3 + 21*b^2 - 21*b + 14, -b^5 + 2*b^4 - 3*b^3 \
            + 2*b^2 - 2*b, 0])
            sage: X.intersection_multiplicity(Y, Q)
            3

        ::

            sage: A.<x,y,z> = AffineSpace(QQ, 3)
            sage: X = A.subscheme([z^2 - 1])
            sage: Y = A.subscheme([z - 1, y - x^2])
            sage: Q = A([1,1,1])
            sage: X.intersection_multiplicity(Y, Q)
            Traceback (most recent call last):
            ...
            TypeError: the intersection of this subscheme and (=Closed subscheme of Affine Space of dimension 3
            over Rational Field defined by: z - 1, -x^2 + y) must be proper and finite

        ::

            sage: A.<x,y,z,w,t> = AffineSpace(QQ, 5)
            sage: X = A.subscheme([x*y, t^2*w, w^3*z])
            sage: Y = A.subscheme([y*w + z])
            sage: Q = A([0,0,0,0,0])
            sage: X.intersection_multiplicity(Y, Q)
            Traceback (most recent call last):
            ...
            TypeError: the intersection of this subscheme and (=Closed subscheme of Affine Space of dimension 5
            over Rational Field defined by: y*w + z) must be proper and finite
        """
        AA = self.ambient_space()
        if AA != X.ambient_space():
            raise TypeError("this subscheme and (=%s) must be defined in the same ambient space"%X)
        W = self.intersection(X)
        try:
            W._check_satisfies_equations(P)
        except TypeError:
            raise TypeError("(=%s) must be a point in the intersection of this subscheme and (=%s)"%(P,X))
        if AA.dimension() != self.dimension() + X.dimension() or W.dimension() != 0:
            raise TypeError("the intersection of this subscheme and (=%s) must be proper and finite"%X)
        I = self.defining_ideal()
        J = X.defining_ideal()
        # move P to the origin and localize
        chng_coords = [AA.gens()[i] + P[i] for i in range(AA.dimension_relative())]
        R = AA.coordinate_ring().change_ring(order="negdegrevlex")
        Iloc = R.ideal([f(chng_coords) for f in I.gens()])
        Jloc = R.ideal([f(chng_coords) for f in J.gens()])
        # compute the intersection multiplicity with Serre's Tor formula using Singular
        singular.lib("homolog.lib")
        i = 0
        s = 0
        t = sum(singular.Tor(i, Iloc, Jloc).std().hilb(2).sage())
        while t != 0:
            s = s + ((-1)**i)*t
            i = i + 1
            t = sum(singular.Tor(i, Iloc, Jloc).std().hilb(2).sage())
        return s

    def multiplicity(self, P):
        r"""
        Return the multiplicity of ``P`` on this subscheme.

        This is computed as the multiplicity of the local ring of this subscheme corresponding to ``P``. This
        subscheme must be defined over a field. An error is raised if ``P`` is not a point on this subscheme.

        INPUT:

        - ``P`` -- a point on this subscheme.

        OUTPUT:

        An integer.

        EXAMPLES::

            sage: A.<x,y,z,w> = AffineSpace(QQ, 4)
            sage: X = A.subscheme([z*y - x^7, w - 2*z])
            sage: Q1 = A([1,1/3,3,6])
            sage: X.multiplicity(Q1)
            1
            sage: Q2 = A([0,0,0,0])
            sage: X.multiplicity(Q2)
            2

        ::

            sage: A.<x,y,z,w,v> = AffineSpace(GF(23), 5)
            sage: C = A.curve([x^8 - y, y^7 - z, z^3 - 1, w^5 - v^3])
            sage: Q = A([22,1,1,0,0])
            sage: C.multiplicity(Q)
            3

        ::

            sage: K.<a> = QuadraticField(-1)
            sage: A.<x,y,z,w,t> = AffineSpace(K, 5)
            sage: X = A.subscheme([y^7 - x^2*z^5 + z^3*t^8 - x^2*y^4*z - t^8])
            sage: Q1 = A([1,1,0,1,-1])
            sage: X.multiplicity(Q1)
            1
            sage: Q2 = A([0,0,0,-a,0])
            sage: X.multiplicity(Q2)
            7
        """
        if not self.base_ring() in Fields():
            raise TypeError("subscheme must be defined over a field")

        # Check whether P is a point on this subscheme
        try:
            P = self(P)
        except TypeError:
            raise TypeError("(=%s) is not a point on (=%s)"%(P,self))

        # Apply a linear change of coordinates to self so that P is sent to the origin
        # and then compute the multiplicity of the local ring of the translated subscheme
        # corresponding to the point (0,...,0)
        AA = self.ambient_space()
        chng_coords = [AA.gens()[i] + P[i] for i in range(AA.dimension_relative())]
        R = AA.coordinate_ring().change_ring(order='negdegrevlex')
        I = R.ideal([f(chng_coords) for f in self.defining_polynomials()])
        return singular.mult(singular.std(I)).sage()
