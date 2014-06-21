r"""
Commutative Differential Graded Algebras.


AUTHORS:

- Miguel Marco (2014-06-20): initial version

Commutative graded differential algebras are graded commutative algebras endowed
with a differential.


EXAMPLES:

The CDGA is constructed with the command CDGAlgebra::

    sage: A = CDGAlgebra(QQ, 'x,y,z')
    sage: A
    Commutative Graded Differential Algebra over Rational Field with generators ('x', 'y', 'z')

If the degrees are not given, they are assumed to be all one::

    sage: A.<x,y,z> = CDGAlgebra(QQ)
    sage: x.degree()
    1
    sage: B.<a,b> = CDGAlgebra(QQ, degrees = (2,3))
    sage: a.degree()
    2
    sage: b.degree()
    3

The differential is defined as a separated object::

    sage: A.<x,y> = CDGAlgebra(QQ)
    sage: D = A.differential({x: x*y, y: x*y})
    sage: D
    Differential map in Commutative Graded Differential Algebra over Rational Field with generators ('x', 'y')
    sending:
        x --> x*y
        y --> x*y

This way the same algebra structire can be used to compute cohomologies
in different cases::

    sage: A.<x,y,z,t> = CDGAlgebra(QQ)
    sage: D = A.differential({x: x*y, y: x*y, z: z*t, t: t*z})
    sage: [A.cohomology(i, D).dimension() for i in range(6)]
    [1, 2, 1, 0, 0, 0]
    sage: D2 = A.differential({x:0, y: 0, z: 0, t: x*y})
    sage: [A.cohomology(i, D2).dimension() for i in range(6)]
    [1, 3, 4, 3, 1, 0]




"""

#*****************************************************************************
#       Copyright (C) 2014 Miguel Marco <mmarco@unizar.es>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.rings.quotient_ring import QuotientRing_nc
from sage.rings.quotient_ring_element import QuotientRingElement
from sage.misc.functional import is_odd
from sage.algebras.free_algebra import FreeAlgebra
from sage.combinat.integer_vector_weighted import WeightedIntegerVectors
from sage.interfaces.singular import is_SingularElement
from sage.rings.polynomial.term_order import TermOrder
from sage.misc.misc_c import prod
from sage.matrix.constructor import matrix
from sage.structure.unique_representation import UniqueRepresentation
from sage.categories.map import Map
from sage.categories.homset import Hom
from sage.categories.sets_cat import Sets
from sage.misc.functional import coerce



class CDGA_Differential(UniqueRepresentation, Map):
    r"""
    The class of differentials over CDGA's

    EXAMPLES::

        sage: A.<x,y>=CDGAlgebra(QQ, degrees=(1,1))
        sage: D = A.differential({x: x*y, y: x*y})
        sage: D
        Differential map in Commutative Graded Differential Algebra over Rational Field with generators ('x', 'y')
        sending:
            x --> x*y
            y --> x*y
        sage: D(x)
        x*y

    """
    def __init__(self, A, im_gens):
        r"""
        Python constructor

        INPUT:

        - ``A`` -- The algebra where the differential is defined.

        - ``im_gens`` -- A tuple containing the image of each generator.

        EXAMPLES::

            sage: A.<x,y,z,t> = CDGAlgebra(QQ)
            sage: D = A.differential({x: x*y, y: x*y, z: z*t, t: t*z})
            sage: [A.cohomology(i, D).dimension() for i in range(6)]
            [1, 2, 1, 0, 0, 0]
            sage: D2 = A.differential({x:0, y: 0, z: 0, t: x*y})
            sage: [A.cohomology(i, D2).dimension() for i in range(6)]
            [1, 3, 4, 3, 1, 0]
            sage: D
            Differential map in Commutative Graded Differential Algebra over Rational Field with generators ('x', 'y', 'z', 't')
            sending:
                x --> x*y
                y --> x*y
                z --> z*t
                t --> -z*t
            sage: D2
            Differential map in Commutative Graded Differential Algebra over Rational Field with generators ('x', 'y', 'z', 't')
            sending:
                x --> 0
                y --> 0
                z --> 0
                t --> x*y

        """
        Map.__init__(self, Hom(A, A, category=Sets()))
        dic = {A.gen(i): im_gens[i] for i in range(A.ngens())}
        self._dic_ = dic
        for i in dic.keys():
            if not dic[i].is_zero():
                if (not dic[i].is_homogeneous() or
                        dic[i].degree() != i.degree()+1):
                    raise ValueError("The given dictionary does not determine a degree 1 map")
            for i in A.gens():
                if not self(self(i)).is_zero():
                    raise ValueError("The given dictionary does not determine a valid differential")

    def _call_(self, x):
        r"""
        Apply the differential to ``x``.

        TESTS::

            sage: A.<x,y,z,t> = CDGAlgebra(QQ)
            sage: D = A.differential({x: x*y, y: x*y, z: z*t, t: t*z})
            sage: D(x*t+1/2*t*x*y) # indirect doctesting
            -1/2*x*y*z*t + x*y*t + x*z*t

        """
        if x.is_zero():
            return self.codomain().zero()
        res = self.codomain().zero()
        dic = x.dict()
        for key in dic:
            keyl = list(key)
            mon = prod([x.parent().gen(i)**key[i] for i in range(len(key))])
            coef = dic[key]
            while sum(keyl) > 0:
                nonz = [a for a in keyl if a != 0]
                first = keyl.index(nonz[0])
                keyl[first] -= 1
                v1 = self._dic_[x.parent().gen(first)]
                v2 = prod(x.parent().gen(i)**keyl[i] for i in
                          range(len(keyl)))
                res += coef*v1*v2
                coef *= (-1) ** x.parent()._degrees[first] * x.parent().gen(first)
        return res

    def _repr_(self):
        r"""
        TESTS::

            sage: A.<x,y,z,t> = CDGAlgebra(QQ)
            sage: D = A.differential({x: x*y, y: x*y, z: z*t, t: t*z})
            sage: D
            Differential map in Commutative Graded Differential Algebra over Rational Field with generators ('x', 'y', 'z', 't')
            sending:
                x --> x*y
                y --> x*y
                z --> z*t
                t --> -z*t
            sage: D2 = A.differential({t: x*y})
            sage: D2
            Differential map in Commutative Graded Differential Algebra over Rational Field with generators ('x', 'y', 'z', 't')
            sending:
                x --> 0
                y --> 0
                z --> 0
                t --> x*y
        """
        string = "Differential map in {}".format(self.domain()) + "\n"
        string += "sending:"
        for i in self.domain().gens():
            string += "\n    " + str(i)+" --> "+str(self(i))
        return string



class CDGAElement(QuotientRingElement):
    r"""
    An element of the CDGA.
    """
    def __init__(self, A, rep):
        r"""
        Python constructor

        INPUT:

        - ``parent`` - the CDGA

        - ``rep`` - a representative of the element in `R`; this is used
        as the internal representation of the element

        EXAMPLES::


            sage: B.<x,y>=CDGAlgebra(QQ, degrees=(2, 2))
            sage: a = B({(1,1): -3, (2,5): 1/2})
            sage: a # indirect doctest
            1/2*x^2*y^5 - 3*x*y
            sage: b = x^2*y^3+2
            sage: b
            x^2*y^3 + 2

        """
        QuotientRingElement.__init__(self, A, rep)

    def degree(self):
        r"""
        Return the degree of self.

        It is the maximum of the degrees of its monomials.

        EXAMPLES::

            sage: A.<x,y,z,t>=CDGAlgebra(QQ, degrees = (2,3,3,1))
            sage: el = y*z+2*x*t-x^2*y
            sage: el.degree()
            7
            sage: el.monomials()
            [x^2*y, y*z, x*t]
            sage: [i.degree() for i in el.monomials()]
            [7, 6, 3]

        """
        exps = self.lift().dict().keys()
        degrees = self.parent()._degrees
        n = self.parent().ngens()
        l = [sum(e[i] * degrees[i] for i in range(n)) for e in exps]
        return max(l)

    def  is_homogeneous(self):
        r"""
        Determine if the element is homogenous or not.

        EXAMPLES::

            sage: A.<x,y,z,t>=CDGAlgebra(QQ, degrees = (2,3,3,1))
            sage: el = y*z + 2*x*t - x^2*y
            sage: el.degree()
            7
            sage: el.monomials()
            [x^2*y, y*z, x*t]
            sage: [i.degree() for i in el.monomials()]
            [7, 6, 3]
            sage: el.is_homogeneous()
            False
            sage: em = x^3 - 5*y*z + 3/2*x*z*t
            sage: em.is_homogeneous()
            True
            sage: em.monomials()
            [x^3, y*z, x*z*t]
            sage: [i.degree() for i in em.monomials()]
            [6, 6, 6]

        """
        if self.is_zero():
            return True
        if len(set([a.degree() for a in self.monomials()]))==1:
            return True
        return False

    def differential(self, D):
        r"""
        Return the differential of the element.

        EXAMPLES::

            sage: A.<x,y,z> = CDGAlgebra(QQ)
            sage: D = A.differential({x: y*z, y: x*z, z: x*y})
            sage: x.differential(D)
            y*z
            sage: y.differential(D)
            x*z
            sage: (x + y).differential(D)
            x*z + y*z
            sage: (x*y).differential(D)
            0

        """
        return D(self)

    def dict(self):
        r"""
        A dictionary that determines the element.

        The keys of this dictionary are the tuples of exponents of each monomial,
        and the keys are the corresponding coefficients.

        EXAMPLES::

            sage: A.<x,y,z,t> = CDGAlgebra(QQ, degrees=(1, 2, 2, 3))
            sage: dic = (x*y - 5*y*z + 7*x*y^2*z^3*t).dict()
            sage: sorted(dic.items())
            [((0, 1, 1, 0), -5), ((1, 1, 0, 0), 1), ((1, 2, 3, 1), 7)]

        """
        return self.lift().dict()



class CDGAlgebra(UniqueRepresentation, QuotientRing_nc):
    r"""
    The class of CDGA's

    """
    Element = CDGAElement

    @staticmethod
    def __classcall__(cls, base, names, degrees=None, R = None, I = None):
        r"""
        Normalize the input for the __init__ method.

        INPUT:

        - ``base```-- The base ring of the algebra.

        - ``names`` -- The names of the variables.

        - ``degrees`` -- The degrees of the generators. By default they are set to 1.

        - ``R`` -- An underlying g-algebra. Only meant to be used by the
        .quotient method.

        - ``I`` -- A twosided ideal in R, with the desired relations. Only
        meant to be used by the .quotient method.

        TESTS::

        sage: CDGAlgebra(GF(2), 'x,y', (3, 6)) # indirect doctesting
        Commutative Graded Differential Algebra over Finite Field of size 2 with generators ('x', 'y')

        """
        if isinstance(names, basestring):
            n = len(names.split(','))
        else:
            n = len(names)
            names = tuple(names)
        if not degrees:
            degrees = tuple([1 for i in range(n)])
        else:
            degrees = tuple(degrees)
        if not R or not I:
            F = FreeAlgebra(base, n, names)
            gens = F.gens()
            rels = {}
            for i in range(len(gens)-1):
                for j in range(i+1, len(gens)):
                    rels[gens[j]*gens[i]] = (-1) ** (degrees[i] * degrees[j]) * gens[i] * gens[j]
            R = F.g_algebra(rels, order = TermOrder('wdegrevlex', degrees))
            I = R.ideal([R.gen(i)**2 for i in range(n) if is_odd(degrees[i])], side='twosided')
        return super(CDGAlgebra, cls).__classcall__(cls, base, R, I, names, degrees)


    def __init__(self, base, R, I, names, degrees=None):
        """
        Python constructor

        INPUT:

        - ``R`` -- The ring over which the algebra is defined.

        - ``I`` -- An ideal over the corresponding g-algebra. It is meant to
        include, among other relations, the squares of the generators of
        odd degree.

        - ``names`` -- The names of the generators

        -- ``differential`` -- A tuple of dictionaries. They determine the
        differential of each generator.

        EXAMPLES::

            sage: A.<x,y,z,t> = CDGAlgebra(QQ)
            sage: A # indirect doctesting
            Commutative Graded Differential Algebra over Rational Field with generators ('x', 'y', 'z', 't')
            sage: CDGAlgebra(QQ, ('x','y','z'), [3,4,2])
            Commutative Graded Differential Algebra over Rational Field with generators ('x', 'y', 'z')
            sage: CDGAlgebra(QQ, ('x','y','z', 't'), [3, 4, 2, 1])
            Commutative Graded Differential Algebra over Rational Field with generators ('x', 'y', 'z', 't')

        """
        self._degrees = tuple(degrees)
        QuotientRing_nc.__init__(self, R, I, names)




    def _repr_(self):
        """
        TESTS::

            sage: CDGAlgebra(QQ, ('x','y','z', 't'), [3, 4, 2, 1]) # indirect doctesting
            Commutative Graded Differential Algebra over Rational Field with generators ('x', 'y', 'z', 't')

        """
        return "Commutative Graded Differential Algebra over {} with generators {}".format(self.base_ring(), self._names)


    def homogeneous_part(self, n):
        """
        Return a basis of the n'th homogeneous part of the algebra.

        EXAMPLES::

            sage: A.<x,y,z,t> = CDGAlgebra(QQ, degrees=(1, 2, 2, 3))
            sage: A.homogeneous_part(2)
            [z, y]
            sage: A.homogeneous_part(3)
            [t, x*z, x*y]
            sage: A.homogeneous_part(4)
            [x*t, z^2, y*z, y^2]
            sage: A.homogeneous_part(5)
            [z*t, y*t, x*z^2, x*y*z, x*y^2]
            sage: A.homogeneous_part(6)
            [x*z*t, x*y*t, z^3, y*z^2, y^2*z, y^3]

        """
        weivec = WeightedIntegerVectors(n, self._degrees)
        basis = []
        for v in weivec:
            el = prod([self.gen(i)**v[i] for i in range(len(v))])
            di = el.dict()
            if len(di) == 1:
                if list(di.keys()[0]) == v:
                    basis.append(el)
        return basis



    def differential_matrix(self, n, D):
        """
        Return the matrix that gives the differential on the n'th degree.

        INPUT:

        - ``n`` -- Integer

        - ``D`` -- The differential.

        EXAMPLES::

            sage: A.<x,y,z,t> = CDGAlgebra(GF(5), degrees=(2, 3, 2, 4))
            sage: D = A.differential({t: x*y, x: y, z: y})
            sage: A.differential_matrix(4, D)
            [0 1]
            [2 0]
            [1 1]
            [0 2]
            sage: A.homogeneous_part(4)
            [t, z^2, x*z, x^2]
            sage: A.homogeneous_part(5)
            [y*z, x*y]
            sage: D(t)
            x*y
            sage: D(z^2)
            2*y*z
            sage: D(x*z)
            x*y + y*z
            sage: D(x^2)
            2*x*y


        """
        dom = self.homogeneous_part(n)
        cod = self.homogeneous_part(n+1)
        cokeys = [a.lift().dict().keys()[0] for a in cod]
        m = matrix(self.base_ring(), len(dom), len(cod))
        for i in range(len(dom)):
            im = dom[i].differential(D)
            dic = im.lift().dict()
            for j in dic.keys():
                k = cokeys.index(j)
                m[i,k] = dic[j]
        return m

    def quotient(self, I, check=True):
        """
        Create the quotient of this algebra by a twosided ideal ``I``


        EXAMPLES::

            sage: A.<x,y,z,t> = CDGAlgebra(GF(5), degrees=(2, 3, 2, 4))
            sage: I = A.ideal([x*t+y^2, x*z - t])
            sage: B = A.quotient(I)
            sage: B
            Commutative Graded Differential Algebra over Finite Field of size 5 with generators ('x', 'y', 'z', 't')
            sage: B(x*t)
            0
            sage: B(x*z)
            t
            sage: A.homogeneous_part(7)
            [y*t, y*z^2, x*y*z, x^2*y]
            sage: B.homogeneous_part(7)
            [y*t, y*z^2, x^2*y]


        """
        NCR = self._QuotientRing_nc__R
        gens1 = list(self._QuotientRing_nc__I.gens())
        gens2 = [i.lift() for i in I.gens()]
        J = NCR.ideal(gens1+gens2, side='twosided')
        if check:
            for i in I.gens():
                if not i.is_homogeneous():
                    raise ValueError("The ideal must be homogeneous")
        CA = CDGAlgebra(self.base_ring(), self._names, self._degrees,  NCR, J)
        return CA


    def _element_constructor_(self, x):
        r"""
        TESTS::

            sage: A.<x,y,z,t> = CDGAlgebra(QQ, degrees = (2, 3, 2, 4))
            sage: A({(1,0,3,1): 2, (2,1,2,2): 3}) # indirect doctesting
            3*x^2*y*z^2*t^2 + 2*x*z^3*t
            sage: A.<x,y,z,t> = CDGAlgebra(GF(5))
            sage: A({(1,0,3,1): 2, (2,1,2,2): 3})
            0
        """
        if isinstance(x, QuotientRingElement):
            if x.parent() is self:
                return x
            x = x.lift()
        if is_SingularElement(x):
            #self._singular_().set_ring()
            x = self.element_class(self, x.sage_poly(self.cover_ring()))
            return x
        if isinstance(x, dict):
            res = self.zero()
            for i in x.keys():
                mon = prod(self.gen(j)**i[j] for j in range(len(i)))
                res += x[i]*mon
            return res
        if coerce:
            R = self.cover_ring()
            x = R(x)
        return self.element_class(self, x)

    def cohomology(self, n, differential):
        """
        Return the n'th cohomology group of the algebra. This is in fact a
        vector space over the base ring.

        INPUT:

        - ``n`` -- The degree in which the cohomology is computed.

        - ``differential`` -- The differential to be used.

        EXAMPLES::

            sage: A.<x,y,z,t> = CDGAlgebra(QQ, degrees = (2, 3, 2, 4))
            sage: D = A.differential({t: x*y, x: y, z: y})
            sage: A.cohomology(4, D)
            Vector space quotient V/W of dimension 2 over Rational Field where
            V: Vector space of degree 4 and dimension 2 over Rational Field
            Basis matrix:
            [   1    0    0 -1/2]
            [   0    1   -2    1]
            W: Vector space of degree 4 and dimension 0 over Rational Field
            Basis matrix:
            []



        """
        N = self.differential_matrix(n, differential)
        if n == 1:
            M = matrix(self.base_ring(), 0, N.nrows())
        else:
            M = self.differential_matrix(n-1, differential)
        V0 = self.base_ring()**M.nrows()
        V1 = self.base_ring()**M.ncols()
        V2 = self.base_ring()**N.ncols()
        h1 = V0.Hom(V1)(M)
        h2 = V1.Hom(V2)(N)
        return h2.kernel().quotient(h1.image())

    def differential(self, dic):
        r"""
        Return the differential of self determined by the dictionary ``dic``

        INPUT:

        - ``dic`` -- A dictionary with the images of the generators.

        The generators that don't appear in the dictionary are sent to zero.

        EXAMPLES::

            sage: A.<x,y,z,t> = CDGAlgebra(QQ, degrees = (2, 3, 2, 4))
            sage: D = A.differential({t: x*y, x: y, z: y})
            sage: D
            Differential map in Commutative Graded Differential Algebra over Rational Field with generators ('x', 'y', 'z', 't')
            sending:
                x --> y
                y --> 0
                z --> y
                t --> x*y

        """
        aux = [self.zero() for i in range(self.ngens())]
        for k in dic:
            i = self.gens().index(k)
            aux[i] = coerce(self, dic[k])
        return CDGA_Differential(self, tuple(aux))




