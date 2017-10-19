## TODO
# 1) Most classes and methods need examples, tests, I/O documentation

import warnings

from sage.structure.sage_object import SageObject
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
# from sage.rings.polynomial.polynomial_element import *
# from sage.rings.all import *
from sage.rings.ring import Ring
from sage.symbolic.ring import SymbolicRing
# from sage.rings.complex_field import *
from sage.rings.real_mpfr import RR
from sage.misc.flatten import flatten
from sage.rings.polynomial.laurent_polynomial_ring import LaurentPolynomialRing

## CLASSES

class PolynomialSystem(SageObject):
    """
    A class for systems of polynomials, primarily for numerical purposes.
    """
    def __init__(self, polys, var_order=None):
        """
        This is a constructor that takes a list of polynomials and
        returns an object of class PolynomialSystem.

        EXAMPLES::

            sage: from sage.numerical.phcpy_interface import PolynomialSystem
            sage: R.<x,y>=PolynomialRing(QQ,2)
            sage: F=[x^2-y,x+y]
            sage: A=PolynomialSystem(F)
            sage: B=PolynomialSystem(F,var_order=[y,x])
            sage: A.ring
            Multivariate Laurent Polynomial Ring in x, y over Rational Field
            sage: B.ring
            Multivariate Laurent Polynomial Ring in y, x over Rational Field

            sage: S.<a,b,c>=PolynomialRing(RealField(100),3)
            sage: G=[(a+b+c)^2-a-b-c,(a-b)^2+(b-c)^2]
            sage: C=PolynomialSystem(G,var_order=[b,a,c])
            sage: C.ring
            Multivariate Laurent Polynomial Ring in b, a, c over Real Field \
            with 100 bits of precision
            sage: var('w,z')
            (w, z)
            sage: D=PolynomialSystem([w^7-z])
            doctest:...: RuntimeWarning: SymbolicRing expressions not checked \
            for consistency. Precision may be lost due to conversion \
            of rationals.
            sage: D.ring
            Multivariate Laurent Polynomial Ring in z, w over Complex \
            Field with 64 bits of precision
            sage: D.polys
            [w^7 - z]

        """
        if isinstance(polys, PolynomialSystem):
            self = polys
        else:
            if not isinstance(polys, list):
                raise TypeError("incorrect input")
            ring_list = list(set([p.parent() if initially_symbolic == False \
                else list(p.variables()) for p in polys]))
            if len(ring_list) != 1 or not isinstance(ring_list[0], Ring):
                raise TypeError("polynomials don't have same parent ring")
            initially_symbolic = isinstance(ring_list[0].base_ring(), SymbolicRing)
            if initially_symbolic:
                good_base_ring = ComplexField(64)
                warnings.warn("SymbolicRing expressions not checked for consistency. \
                Precision may be lost due to conversion of rationals.", RuntimeWarning)
            elif ring_list[0].base_ring() in set([RR, CC, QQ, ZZ]) \
                 or ring_list[0].base_ring().parent() == RR.parent() or \
                 ring_list[0].base_ring().parent() == CC.parent():
                good_base_ring = ring_list[0].base_ring()
            else:
                raise TypeError("coefficient ring")
            myvars = list(reversed(list(set(flatten([list(p.parent().gens()) for p in polys])))))
            if var_order == None:
                var_order = myvars
            if var_order != None:
                if set(var_order) == set(myvars):
                    myvars = var_order
                    self.ring = LaurentPolynomialRing(good_base_ring, len(myvars), myvars)
                    self.polys = [(self.ring)(p)  for p in polys]
                else:
                    raise TypeError("Variable order is not the exact list of variables involved")
            if self.ring.base_ring() != QQ and self.ring.base_ring() != ZZ:
                self.prec = self.ring.base_ring().precision()
    def evaluate(self, npoint):
        if isinstance(npoint, list):
            npoint = NumericalPoint(npoint, ring=self.ring)
        if not isinstance(npoint, NumericalPoint):
            raise TypeError("point provided must be of type NumericalPoint")
        return([f.subs(npoint.to_dict(temp_ring=self.ring)) for f in self.polys])
    def __str__(self):
        return "%s over %s. " %(self.polys, self.ring)
    def __repr__(self):
        return "%s over %s. " %(self.polys, self.ring)
class NumericalPoint(SageObject):
    """
    A class for representing points numerically
    """
    def __init__(self, coords, ring=None, multiplicity=None, rco=None, err=None, res=None):
        """
        Construct from list of coordinates.

        EXAMPLES::

            sage: from sage.numerical.phcpy_interface import NumericalPoint
            sage: R.<x,y,z> = PolynomialRing(CC,3)
            sage: p = NumericalPoint([2,3,4+I],ring=R)
            sage: p.dict
            {z: I + 4, y: 3, x: 2}
            sage: p.coordinates
            [2, 3, I + 4]
        """
        if isinstance(coords, list):
            self.coordinates = coords
            self.ring = ring
            self.multiplicity = multiplicity
            self.rco = rco
            self.err = err
            self.res = res
            if ring != None:
                self.dict = self.to_dict()
            else:
                self.dict = None
        if isinstance(coords, dict):
            ring_list = list(set([p.parent() for p in coords.keys()]))
            if len(ring_list) != 1 or not isinstance(ring_list[0], Ring):
                raise TypeError("Keys are not all in the same ring")
            if len(coords.keys()) != coords.keys()[0].parent().ngens():
                raise TypeError("Number of keys disagrees with number of variables of ring")
            self.ring = coords.keys()[0].parent()
            self.coordinates = [coords[i] for i in self.ring.gens()]
            self.multiplicity = multiplicity
            self.rco = rco
            self.err = err
            self.res = res
            self.dict = coords
        # and so on as more args are added
    def to_dict(self, temp_ring=None):
        if self.ring != None and temp_ring is None:
            temp_ring = self.ring
        if temp_ring != None:
            new_dictionary = dict([(temp_ring.gens()[i], \
            self.coordinates[i]) for i in range(0, \
            len(self.coordinates))])
        if not self.ring is None:
            self.dict = new_dictionary
            return new_dictionary
        else:
            raise AttributeError("please set a ring")
    def __str__(self):
        return "A numerical point in CC^%s." %(len(self.coordinates))
    def __repr__(self):
        return "A numerical point in CC^%s." %(len(self.coordinates))

class WitnessSet(SageObject):
    def __init__(self, poly_sys, forms, points):
        """
        INPUT:
            *) poly_sys, an object of type PolynomialSystem
            *) forms, an object of type PolynomialSystem consisting of linear \
            forms w/ the same ring as poly_sys
            *) points --- a list of objects of type NumericalPoint

        EXAMPLES:

            sage: from sage.numerical.phcpy_interface import PolynomialSystem
            sage: from sage.numerical.phcpy_interface import NumericalPoint
            sage: from sage.numerical.phcpy_interface import WitnessSet
            sage: R.<x,y>=PolynomialRing(QQ,2)
            sage: F=PolynomialSystem([y-x^2])
            sage: L=PolynomialSystem([y-25])
            sage: pts=[NumericalPoint([5,25]),NumericalPoint([-5,25])]
            sage: W=WitnessSet(F,L,pts)
            sage: W.check_validity()
            True
        """
        if not isinstance(poly_sys, PolynomialSystem):
            raise TypeError("first argument should be a PolynomialSystem")
        if not isinstance(forms, PolynomialSystem):
            raise TypeError("second argument should be a polynomial system")
        points = [NumericalPoint(p, ring=poly_sys.ring) \
                  if isinstance(p, list) else p for p in points]
        if not (isinstance(points, list) and len(set([p.parent() \
            for p in points])) == 1 and isinstance(points[0], NumericalPoint)):
            raise TypeError("third argument should be a list of NumericalPoints")
        if not poly_sys.ring == forms.ring:
            raise TypeError("make sure first two arguments share a common ring")
        self.polynomials = poly_sys
        self.slices = forms
        self.points = points
        self.dimension = len(forms.polys)
    def check_validity(self, tolerance=0.000001):
        evaluations = flatten([[(f.evaluate(p)) for f in \
            ([self.polynomials, self.slices])] for p in self.points])
        if False in set([(abs(q) < tolerance) for q in evaluations]):
            return False
        else:
            return True
    def __str__(self):
        return "A witness set for a dimension-%s component with %s points." \
            %(self.dimension, len(self.points))
    def __repr__(self):
        return "A witness set for a dimension-%s component with %s points." \
            %(self.dimension, len(self.points))


class NumericalIrreducibleDecomposition(SageObject):
    """
    a class which organizes the witness sets appearing in the NID of a variety
    """
    def __init__(self):
        self.components = dict()
    def append_witness_set(self, wset):
        if not isinstance(wset, WitnessSet):
            raise TypeError("must append with a witness set")
        if wset.dimension in self.components.keys():
            self.components[wset.dimension].append(wset)
        else:
            self.components[wset.dimension] = [wset]
    def __str__(self):
        return_string = ""
        for i in self.components.keys():
            return_string += "Dimension "+str(i) + ":" + "\n"
            for j in self.components[i]:
                return_string += "    Component of degree "+str(len(j.witness_points)) + "\n "
        return return_string
    def __repr__(self):
        return_string = ""
        for i in self.components.keys():
            return_string += "Dimension "+str(i) + ":" + "\n"
            for j in self.components[i]:
                return_string += "    Component of degree "+str(len(j.witness_points)) + "\n "
        return return_string


class ParametrizedPolynomialSystem(PolynomialSystem):
    def __init__(self, system, params):
        super(ParametrizedPolynomialSystem, self).__init__(system)
        if (not isinstance(params, list)) or (False in \
            set([g in self.ring.gens() for g in params])):
            raise TypeError("Parameters must be a list of variables in the ring.")
        self.params = params
        self.variables = (set(self.ring.gens())).difference(self.params)
    def specialize(self, sub_dict, specialize_ring=True):
        if False in set([g in sub_dict.keys() for g in self.params]):
            raise TypeError("Specialization keys should be parameters.")
        special_self = self
        special_self.polys = [f.subs(sub_dict) for f in self.polys]
        if specialize_ring:
            new_vars = list(set(self.ring.gens()).difference(sub_dict.keys()))
            specializedring = PolynomialRing(self.polys[0].base_ring(), len(new_vars), new_vars)
            special_self.polys = [specializedring(str(f)) for f in special_self.polys]
            special_self.ring = specializedring
        return special_self
