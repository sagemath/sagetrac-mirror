"""
Divisors

Sage provides extensive computations with divisors on global function fields.

EXAMPLES:

The divisor of an element of the function field is the formal sum of poles and zeros
of the element with multiplicities::

    sage: K.<x> = FunctionField(GF(2)); R.<t> = K[]
    sage: L.<y> = K.extension(t^3 + x^3*t + x)
    sage: f = x/(y+1)
    sage: f.divisor()
    -1*Place (1/x, 1/x^3*y^2 + 1/x)
     + Place (1/x, 1/x^3*y^2 + 1/x^2*y + 1)
     + 3*Place (x, y)
     - Place (x^3 + x + 1, y + 1)

The Riemann-Roch space of a divisor can be computed. We can get a basis of
the space as a vector space over the constant field::

    sage: p = L.places_finite()[0]
    sage: q = L.places_infinite()[0]
    sage: (3*p + 2*q).basis_function_space()
    [1/x*y^2 + x^2, 1, 1/x]

We verify the Riemann-Roch theorem::

    sage: D = 3*p - q
    sage: index_of_speciality = len(D.basis_differential_space())
    sage: D.dimension() == D.degree() - L.genus() + 1 + index_of_speciality
    True

AUTHORS:

- Kwankyu Lee (2017-04-30): initial version

"""
#*****************************************************************************
#       Copyright (C) 2016 Kwankyu Lee <ekwankyu@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import absolute_import

from sage.misc.cachefunc import cached_method
from sage.misc.lazy_import import lazy_import

from sage.arith.all import lcm

from sage.structure.parent import Parent
from sage.structure.element import ModuleElement
from sage.structure.richcmp import richcmp

from sage.categories.homset import Hom
from sage.categories.morphism import SetMorphism
from sage.categories.commutative_additive_groups import CommutativeAdditiveGroups

from sage.matrix.constructor import matrix

from sage.modules.free_module_element import vector

from sage.rings.integer_ring import IntegerRing
from sage.rings.integer import Integer

from .place import FunctionFieldPlace, PlaceSet

def zero_divisor(field):
    """
    Construct a zero divisor.

    INPUT:

    - ``field`` -- function field

    EXAMPLES::

        sage: K.<x> = FunctionField(GF(2)); R.<t> = K[]
        sage: F.<y> = K.extension(t^3 - x^2*(x^2 + x + 1)^2)
        sage: from sage.rings.function_field.divisor import zero_divisor
        sage: z = zero_divisor(F)
        sage: z == z + z
        True
    """
    return FunctionFieldDivisor(field, {})

def prime_divisor(field, place, m=1):
    """
    Construct a prime divisor from the place.

    INPUT:

    - ``field`` -- function field

    - ``place`` -- place of the function field

    - ``m`` -- (default: 1) a positive integer; multiplicity at the place

    EXAMPLES::

        sage: K.<x> = FunctionField(GF(2)); R.<t> = K[]
        sage: F.<y> = K.extension(t^3 - x^2*(x^2 + x + 1)^2)
        sage: p = F.places()[0]
        sage: from sage.rings.function_field.divisor import prime_divisor
        sage: d = prime_divisor(F, p)
        sage: 3 * d == prime_divisor(F, p, 3)
        True
    """
    return FunctionFieldDivisor(field, {place: Integer(m)})

def is_Divisor(x):
    """
    Return ``True`` if ``x`` is a function field divisor.

    EXAMPLES::

        sage: K.<x> = FunctionField(GF(2)); R.<t> = PolynomialRing(K)
        sage: F.<y> = K.extension(t^3-x^2*(x^2+x+1)^2)
        sage: p = F.places()[0]
        sage: from sage.rings.function_field.divisor import prime_divisor
        sage: d = prime_divisor(F, p)
        sage: from sage.rings.function_field.divisor import is_Divisor
        sage: is_Divisor(3 * d)
        True
    """
    return isinstance(x, FunctionFieldDivisor)

class FunctionFieldDivisor(ModuleElement):
    """
    Divisors of function fields.
    """
    def __init__(self, field, data):
        """
        Initialize.

        INPUT:

        - ``field`` -- functon field

        - ``data`` -- dict of place and multiplicity pairs

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); R.<t> = K[]
            sage: F.<y> = K.extension(t^3 - x^2*(x^2 + x + 1)^2)
            sage: f = x/(y+1)
            sage: f.divisor()
            Place (1/x, 1/x^4*y^2 + 1/x^2*y + 1)
             + Place (1/x, 1/x^2*y + 1)
             + 3*Place (x, (1/(x^3 + x^2 + x))*y^2)
             - 6*Place (x + 1, y + 1)
        """
        ModuleElement.__init__(self, field.divisor_group())
        self._field = field
        self._data = data

    def _repr_(self, split=True):
        """
        Return a string representation of the divisor.

        INPUT:

        - ``split`` -- boolean; if ``True``, split at the end of each place

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); R.<t> = PolynomialRing(K)
            sage: F.<y> = K.extension(t^3-x^2*(x^2+x+1)^2)
            sage: f = x/(y+1)
            sage: d = f.divisor()
            sage: d._repr_(split=False)
            'Place (1/x, 1/x^4*y^2 + 1/x^2*y + 1) + Place (1/x, 1/x^2*y + 1)
            + 3*Place (x, (1/(x^3 + x^2 + x))*y^2) - 6*Place (x + 1, y + 1)'
        """
        mul = '*'
        plus = ' + '
        minus = ' - '

        if split:
            cr = '\n'
        else:
            cr = ''

        places = sorted(self._data.keys())

        if len(places) == 0:
            return '0'

        p = places.pop(0)
        m = self._data[p]
        if m == 1:
            r = repr(p)
        elif m == -1:
            r = '-1*' + repr(p) # seems more readable then `-`
        else: # nonzero
            r = repr(m) + mul + repr(p)
        for p in places:
            m = self._data[p]
            if m == 1:
                r += cr + plus + repr(p)
            elif m == -1:
                r += cr + minus + repr(p)
            elif m > 0:
                r += cr + plus + repr(m) + mul + repr(p)
            elif m < 0:
                r += cr + minus + repr(-m) + mul + repr(p)
        return r

    def _richcmp_(self, other, op):
        """
        Compare the divisor and the other divisor with respect to the operator.

        INPUT:

        - ``other`` -- divisor

        - ``op`` -- comparison operator

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(4)); _.<Y>=K[]
            sage: L.<y> = K.extension(Y^3+x+x^3*Y)
            sage: pls1 = L.places()
            sage: D1 = pls1[0] + pls1[1]
            sage: D2 = pls1[1] + 2*pls1[2]
            sage: (D1 < D2) == (not D2 < D1)
            True
            sage: D1 + D2 == D2 + D1
            True
        """
        return richcmp(self._data, other._data, op)

    def _neg_(self):
        """
        Return the additive inverse of the divisor.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(4)); _.<Y>=K[]
            sage: L.<y> = K.extension(Y^3+x+x^3*Y)
            sage: f = x/(y+1)
            sage: D = f.divisor()
            sage: D
            -1*Place (1/x, 1/x^3*y^2 + 1/x)
             + Place (1/x, 1/x^3*y^2 + 1/x^2*y + 1)
             + 3*Place (x, y)
             - Place (x^3 + x + 1, y + 1)
            sage: -D
            Place (1/x, 1/x^3*y^2 + 1/x)
             - Place (1/x, 1/x^3*y^2 + 1/x^2*y + 1)
             - 3*Place (x, y)
             + Place (x^3 + x + 1, y + 1)
        """
        data = {}
        for place in self._data:
            data[place] = -self._data[place]
        return FunctionFieldDivisor(self._field, data)

    def _add_(self, other):
        """
        Add the divisor to the other divisor.

        INPUT:

        - ``other`` -- divisor

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(4)); _.<Y>=K[]
            sage: L.<y> = K.extension(Y^3+x+x^3*Y)
            sage: f = x/(y+1)
            sage: D = f.divisor()
            sage: D + 2*D == 3*D
            True
            sage: 3*D - D == 2*D
            True
        """
        places = set(self.support()).union( set(other.support()))
        data = {}
        for place in places:
            m = self.multiplicity(place) + other.multiplicity(place)
            if m != 0:
                data[place] = m
        return FunctionFieldDivisor(self._field, data)

    def _rmul_(self, i):
        """
        Multiply integer `i` to the divisor.

        INPUT:

        - ``i`` -- integer

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(4)); _.<Y>=K[]
            sage: L.<y> = K.extension(Y^3+x+x^3*Y)
            sage: f = x/(y+1)
            sage: D = f.divisor()
            sage: (-3)*(2*D) == -6*D
            True
        """
        data = {}
        for place in self._data:
            m = i * self._data[place]
            if m != 0:
                data[place] = m
        return FunctionFieldDivisor(self._field, data)

    def dict(self):
        """
        Return the dictionary representing the divisor.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(4)); _.<Y>=K[]
            sage: L.<y> = K.extension(Y^3+x+x^3*Y)
            sage: f = x/(y+1)
            sage: D = f.divisor()
            sage: D.dict()
            {Place (1/x, 1/x^3*y^2 + 1/x): -1,
             Place (1/x, 1/x^3*y^2 + 1/x^2*y + 1): 1,
             Place (x, y): 3,
             Place (x^3 + x + 1, y + 1): -1}
        """
        return self._data

    def list(self):
        """
        Return the list of place and multiplicity pairs of the divisor.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(4)); _.<Y>=K[]
            sage: L.<y> = K.extension(Y^3+x+x^3*Y)
            sage: f = x/(y+1)
            sage: D = f.divisor()
            sage: D.list()
            [(Place (1/x, 1/x^3*y^2 + 1/x), -1),
             (Place (1/x, 1/x^3*y^2 + 1/x^2*y + 1), 1),
             (Place (x, y), 3),
             (Place (x^3 + x + 1, y + 1), -1)]
        """
        return sorted(self._data.items())

    def support(self):
        """
        Return the support of the divisor.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(4)); _.<Y>=K[]
            sage: L.<y> = K.extension(Y^3+x+x^3*Y)
            sage: f = x/(y+1)
            sage: D = f.divisor()
            sage: D.support()
            [Place (1/x, 1/x^3*y^2 + 1/x),
             Place (1/x, 1/x^3*y^2 + 1/x^2*y + 1),
             Place (x, y),
             Place (x^3 + x + 1, y + 1)]
        """
        return sorted(self._data.keys())

    def multiplicity(self, place):
        """
        Return the multiplicity of the divisor at the place.

        INPUT:

        - ``place`` -- place of a function field

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(4)); _.<Y>=K[]
            sage: L.<y> = K.extension(Y^3+x+x^3*Y)
            sage: p1,p2 = L.places()[:2]
            sage: D = 2*p1 - 3*p2
            sage: D.multiplicity(p1)
            2
            sage: D.multiplicity(p2)
            -3
        """
        if not self._data.has_key(place):
            return 0
        return self._data[place]

    def degree(self):
        """
        Return the degree of the divisor.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(4)); _.<Y>=K[]
            sage: L.<y> = K.extension(Y^3+x+x^3*Y)
            sage: p1,p2 = L.places()[:2]
            sage: D = 2*p1 - 3*p2
            sage: D.degree()
            -1
        """
        return sum([p.degree() * m for p, m in self.list()])

    def dimension(self):
        """
        Return the dimension of the Riemann-Roch space of the divisor.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(5)); R.<t> = PolynomialRing(K)
            sage: F.<y> = K.extension(t^2-x^3-1)
            sage: O = F.maximal_order()
            sage: I = O.ideal(x-2)
            sage: P1 = I.divisor().support()[0]
            sage: Pinf = F.places_infinite()[0]
            sage: D = 3*Pinf+2*P1
            sage: D.dimension()
            5
        """
        return len(self.basis_function_space())

    def basis_function_space(self):
        """
        Return a basis of the Riemann-Roch space of the divisor.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(5)); R.<t> = K[]
            sage: F.<y> = K.extension(t^2 - x^3 - 1)
            sage: O = F.maximal_order()
            sage: I = O.ideal(x-2)
            sage: D = I.divisor()
            sage: D.basis_function_space()
            [x/(x + 3), 1/(x + 3)]
        """
        basis,_ = self._function_space()
        return basis

    @cached_method
    def function_space(self):
        """
        Return the vector space of the Riemann-Roch space of the divisor.

        OUTPUT:

        - a vector space, an isomorphism from the vector space
          to the Riemann-Roch space, and its inverse.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(5)); R.<t> = PolynomialRing(K)
            sage: F.<y> = K.extension(t^2-x^3-1)
            sage: O = F.maximal_order()
            sage: I = O.ideal(x-2)
            sage: D = I.divisor()
            sage: V, from_V, to_V = D.function_space()
            sage: all(to_V(from_V(e)) == e for e in V)
            True
        """
        F = self._field
        k = F.constant_base_field()

        basis, coordinates = self._function_space()

        n = len(basis)
        V = k ** n

        def from_V(v):
            return sum(v[i] * basis[i] for i in range(n))

        def to_V(f):
            return vector(coordinates(f))

        from sage.rings.function_field.maps import (
            FunctionFieldLinearMap, FunctionFieldPartiallyDefinedLinearMap)

        mor_from_V = FunctionFieldLinearMap(Hom(V,F), from_V)
        mor_to_V = FunctionFieldPartiallyDefinedLinearMap(Hom(F,V), to_V)

        return V, mor_from_V, mor_to_V

    @cached_method
    def _function_space(self):
        """
        Return an (echelon) basis and coordinates function for the Riemann-Roch
        space of the divisor.

        The return values are cached so that :meth:`basis_function_space` and
        :meth:`function_space` methods give consistent outputs.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(5)); R.<t> = PolynomialRing(K)
            sage: F.<y> = K.extension(t^2-x^3-1)
            sage: O = F.maximal_order()
            sage: I = O.ideal(x-2)
            sage: D = I.divisor()
            sage: basis, coordinates = D._function_space()
            sage: basis
            [x/(x + 3), 1/(x + 3)]
            sage: coordinates(basis[0])
            [1, 0]
            sage: coordinates(basis[1])
            [0, 1]
            sage: coordinates(basis[0] + basis[1])
            [1, 1]
            sage: coordinates((x + 4)/(x + 3))
            [1, 4]
        """
        basis, coordinates_func = self._echelon_basis(self._basis())

        return basis, coordinates_func

    def _basis(self):
        """
        Return a basis of the Riemann-Roch space of the divisor.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(5)); R.<t> = PolynomialRing(K)
            sage: F.<y> = K.extension(t^2-x^3-1)
            sage: O = F.maximal_order()
            sage: I = O.ideal(x-2)
            sage: D = I.divisor()
            sage: D._basis()
            [1/(x + 3), x/(x + 3)]

        This implements Hess' algorithm 6.1 in [Hes2002]_
        """
        F = self._field
        n = F.degree()
        O = F.maximal_order()
        Oinf = F.maximal_order_infinite()

        # Step 1: The ideal I is the inverse of the product of prime ideals attached
        # to the finite places in the divisor while the ideal J corresponds
        # to the infinite places in the divisor. The later steps are basically
        # to compute the intersection of the ideals I and J.
        I = O.ideal(1)
        J = Oinf.ideal(1)
        for p in self._data:
            m = self._data[p]
            if p.is_infinite_place():
                J *= p.prime_ideal() ** (-m)
            else:
                I *= p.prime_ideal() ** (-m)

        # Step 2: construct matrix M of rational functions in x such that
        # M * B == C where B = [b1,b1,...,bn], C =[v1,v2,...,vn]
        V,fr,to = F.vector_space()
        B = matrix([to(b) for b in J.gens_over_base()])
        C = matrix([to(v) for v in I.gens_over_base()])
        M = C * B.inverse()

        # Step 2.5: get the denonimator d of M and set mat = d * M
        den = lcm([e.denominator() for e in M.list()])
        R = den.parent() # polynomial ring
        one = R.one()
        mat = matrix(R, n, [e.numerator() for e in (den*M).list()])
        gens = list(I.gens_over_base())

        # Step 3: transform mat to a weak Popov form, together with gens

        # initialise pivot_row and conflicts list
        pivot_row = [[] for i in range(n)]
        conflicts = []
        for i in range(n):
            bestp = -1
            best = -1
            for c in range(n):
                d = mat[i,c].degree()
                if d >= best:
                    bestp = c
                    best = d

            if best >= 0:
                pivot_row[bestp].append((i,best))
                if len(pivot_row[bestp]) > 1:
                    conflicts.append(bestp)

        # while there is a conflict, do a simple transformation
        while conflicts:
            c = conflicts.pop()
            row = pivot_row[c]
            i,ideg = row.pop()
            j,jdeg = row.pop()

            if jdeg > ideg:
                i,j = j,i
                ideg,jdeg = jdeg,ideg

            coeff = - mat[i,c].lc() / mat[j,c].lc()
            s = coeff * one.shift(ideg - jdeg)

            mat.add_multiple_of_row(i, j, s)
            gens[i] += s * gens[j]

            row.append((j,jdeg))

            bestp = -1
            best = -1
            for c in range(n):
                d = mat[i,c].degree()
                if d >= best:
                    bestp = c
                    best = d

            if best >= 0:
                pivot_row[bestp].append((i,best))
                if len(pivot_row[bestp]) > 1:
                    conflicts.append(bestp)

        # Step 4: build a Riemann-Roch basis from the data in mat and gens.
        # Note that the values mat[i,j].degree() - den.degree() are known as
        # invariants of M.
        basis = []
        for j in range(n):
            i,ideg = pivot_row[j][0]
            for k in range( den.degree() - ideg + 1 ):
                basis.append(one.shift(k) * gens[i])
        # Done!
        return basis

    def _echelon_basis(self, basis):
        """
        that is a helper method to compute an echelonized basis of the subspace
        generated by ``basis`` over `k`.

        The entries of ``basis`` vectors are function field elements, viewed
        as vectors of rational functions over `k`.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(5)); _.<t> = PolynomialRing(K)
            sage: F.<y> = K.extension(t^2-x^3-1)
            sage: D = F.divisor_group()(0)
            sage: echelon_basis, coordinates = D._echelon_basis([x/y,(x+1)/y])
            sage: echelon_basis
            [(x/(x^3 + 1))*y, (1/(x^3 + 1))*y]
            sage: f1, f2 = echelon_basis
            sage: coordinates(f1)
            [1, 0]
            sage: coordinates(f2)
            [0, 1]
            sage: coordinates(x/y)
            [1, 0]
            sage: coordinates((x+1)/y)
            [1, 1]
            sage: 1 * f1 + 1 * f2 == (x + 1)/y
            True
        """
        F = self._field
        k = F.constant_base_field()
        V, fr_V, to_V = F.vector_space()
        n = V.degree()
        m = len(basis)

        vbasis = [to_V(f) for f in basis]

        # function to compute the pivot position
        def pivot(v): # v - nonzero vector over rational functions
            for i in range(n):
                e = v[i]
                if e != 0:
                    return (i,e.numerator().degree() - e.denominator().degree())

        def greater(v,w): # v and w are not equal
            return v[0] < w[0] or v[0] == w[0] and v[1] > w[1]

        # collate rows by their pivot position
        pivot_rows = {}
        for i in range(m):
            p = pivot(vbasis[i])
            if p in pivot_rows:
                pivot_rows[p].append(i)
            else:
                pivot_rows[p] = [i]

        # compute "echelon" basis in which pivot positions
        # decrease strictly
        nbasis = []
        npivots = []
        while pivot_rows:
            pivots = pivot_rows.keys()

            head = pivots[0]
            for p in pivots[1:]:
                if not greater(head,p):
                    head = p

            rows = pivot_rows[head]
            if len(rows) > 1:
                r = rows[0]
                vr = vbasis[r][head[0]]
                cr = vr.numerator().lc() / vr.denominator().lc()
                for i in rows[1:]:
                    vi = vbasis[i][head[0]]
                    ci = vi.numerator().lc() / vi.denominator().lc()
                    vbasis[i] -= ci/cr * vbasis[r]
                    p = pivot(vbasis[i])
                    if p in pivot_rows:
                        pivot_rows[p].append(i)
                    else:
                        pivot_rows[p] = [i]
            nbasis.append(vbasis[rows[0]])
            npivots.append(head)
            del pivot_rows[head]

        def coordinates(f):
            v = to_V(f)
            coords = [k(0) for i in range(m)]
            while v != 0:
                p = pivot(v)
                ind = npivots.index(p) # exception implies x is not in the domain
                w = nbasis[ind]
                cv = v[p[0]].numerator().lc() / v[p[0]].denominator().lc()
                cw = w[p[0]].numerator().lc() / w[p[0]].denominator().lc()
                c = cv/cw
                v -= c * w
                coords[ind] = c
            return coords

        newbasis = [fr_V(f) for f in nbasis]

        return newbasis, coordinates

class DivisorGroup(Parent):
    """
    Groups of divisors of function fields.
    """
    def __init__(self, field):
        """
        Initialize.

        INPUT:

        - ``field`` -- function field

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(5)); _.<t> = PolynomialRing(K)
            sage: F.<y> = K.extension(t^2-x^3-1)
            sage: F.divisor_group()
            Divisor group of Function field in y defined by y^2 + 4*x^3 + 4
        """
        Parent.__init__(self, base=IntegerRing(), category=CommutativeAdditiveGroups())

        self._field = field # function field

    def _repr_(self):
        """
        Return the string representation of the group of divisors.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(5)); _.<t> = PolynomialRing(K)
            sage: F.<y> = K.extension(t^2-x^3-1)
            sage: F.divisor_group()
            Divisor group of Function field in y defined by y^2 + 4*x^3 + 4
        """
        return "Divisor group of %s"%(self._field,)

    def _element_constructor_(self, x):
        """
        Construct a divisor from ``x``.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(5)); _.<t> = PolynomialRing(K)
            sage: F.<y> = K.extension(t^2-x^3-1)
            sage: G = F.divisor_group()
            sage: G(0)
            0
        """
        if x == 0:
            return zero_divisor(self._field)
        raise NotImplementedError

    def _coerce_map_from_(self, S):
        """
        Define coercions.

        EXAMPLES:

        A place is converted to a prime divisor::

            sage: K.<x> = FunctionField(GF(5)); R.<t> = PolynomialRing(K)
            sage: F.<y> = K.extension(t^2-x^3-1)
            sage: O = F.maximal_order()
            sage: I = O.ideal(x+1,y)
            sage: P = I.place()
            sage: y.divisor() + P
            -3*Place (1/x, 1/x^2*y)
             + 2*Place (x + 1, y)
             + Place (x^2 + 4*x + 1, y)
        """
        if isinstance(S, PlaceSet):
            func =  lambda place: prime_divisor(self._field, place)
            return SetMorphism(Hom(S,self), func)

    def function_field(self):
        """
        Return the function field to which the divisor group is attached.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(5)); _.<t> = PolynomialRing(K)
            sage: F.<y> = K.extension(t^2-x^3-1)
            sage: G = F.divisor_group()
            sage: G.function_field()
            Function field in y defined by y^2 + 4*x^3 + 4
        """
        return self._field
