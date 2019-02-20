'''

IndRing

Construct sequential diagrams of finitely generated commutative rings and maps between such diagrams.

'''

from sage.rings.integer_ring import *
from sage.rings.rational_field import *
from sage.rings.finite_rings.integer_mod_ring import *
from sage.categories.homset import *
from sage.rings.power_series_ring import *
from sage.rings.finite_rings.finite_field_constructor import *
from sage.arith.misc import *
#Python/Sage can only compare the addresses of function pointers. It cannot compare if two function pointer point to functions that do the same thing.
#This function is defined now so that its address may be compared to later.
def universal_identity(integer):
    return integer

#User has three options to define an IndRing
#Option 1) Pass a string name of an IndRing we have already constructed. That IndRing is returned by the constructor
#Option 2) Pass a ring without a homs function.
#In this case the constructor builds an IndRing whose self.rings function sends every natural to number user passed ring object
#self.homs function sends every natural number to the identity on the user passed ring
#Option 3) User passes function pointers to a rings and a homs function
class IndRing:
    def __init__(self, rings, homs = None):
        # self.rings is a function from the natural numbers to rings
        # self.homs is a function from the natural numbers to homomorphisms of rings where homs(n): rings(n) -> rings(n+1)
        if homs == None:
            if str(type(rings)) == "<type 'str'>":
                # built-in cases
                if rings == "lazard":
                    self.rings = lambda n: PolynomialRing(ZZ, n, 'U')
                    self.homs = lambda n: (self.rings(n)).hom((self.rings(n+1).gens())[:-1], self.rings(n+1))
                    self.name = "Lazard Ring"
                elif rings == "rational_lazard":
                    self.rings = lambda n: PolynomialRing(QQ, n, 'U')
                    self.homs = lambda n: (self.rings(n)).hom((self.rings(n+1).gens())[:-1], self.rings(n+1))
                    self.name = "Rational Lazard Ring"
            else: # assume the argument given is a ring
                self.storageVar = rings
                self.storageVar2 = End(self.storageVar)
                self.rings = lambda n: self.storageVar
                self.homs = lambda n: self.storageVar2.identity()
        else: #rings and homs are pointers to functions as described above
            self.rings = rings
            self.homs = homs

    def __repr__(self):
        if self.id != None:
            return self.name
        else:
            return "User-defined IndRing" #TODO: figure out if there is a __name__ method
#This data type is a map between IndRings
#Domain and codomain are both IndRings
#Maps is a function pointer to a function from the natural numbers to homomorphisms from rings that make up the domain to rings that make up the codomain
#The default implementation is for maps(i) to be a homomorphism from domain.rings(i) to codomain.rings(i)
#rstep and sstep are optional parameters that allow the user to reparameterize the domain and codomain for the homomorphisms from the map functions
#maps(i) is a homomorphism from domain.rings(rstep(i)) to codomain.rings(sstep(i))
#It is assumed that the user is careful enough that the appropriate squares commute
class IndRingHom:
    def __init__(self, domain, codomain, maps, rstep = universal_identity, sstep = universal_identity):
        self.domain = domain
        self.codomain = codomain
        self.maps = maps
        self.rstep = rstep
        self.sstep = sstep

    def __call__(series):
        return series.change_base(self)

Lazard = IndRing('lazard')
Rational_Lazard = IndRing('rational_lazard')
Integer_Ring = IndRing(ZZ)
Rational_Ring = IndRing(QQ)
Integer_7_Ring = IndRing(Integers(7))

