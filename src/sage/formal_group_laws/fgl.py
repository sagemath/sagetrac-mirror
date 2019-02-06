from sage.rings.integer_ring import *
from sage.rings.rational_field import *
from sage.rings.finite_rings.integer_mod_ring import *
from sage.categories.homset import *
from sage.rings.power_series_ring import *
from sage.rings.finite_rings.finite_field_constructor import *
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

#Used to multiply power series faster.
def fast_mult_pow(poly1, poly2, prec):
    val = poly1.valuation()
    terms = poly1.coefficients()
    ans = 0
    for mon in terms.keys():
        ans += (terms[mon] * mon * poly2.add_bigoh(prec + 1 - mon.degree()))
    return ans


class FPS:
    def __init__(self, ring, n, coeffs, b_func = universal_identity, var = 'x', const = None, children = None, operation = None):
        self.ring = ring              # IndRing in which the formal power series resides
        self.n = n                    # Number of variables of the formal power series
        self.b_func = b_func          # Takes total variable precision as argument, and returns how far up into the ring chain one must go; default is identity
        self.coeffs = coeffs          # Coefficient function; takes the exponents of variables in list form as the argument and returns the coefficient of the
                                           # corresponding monomial
        self.var = var                # String which represents the formal power series variables when the formal power series is viewed; default is 'x'
                                           # Note: if len(self.var) == self.n and self.var has no duplicate characters, then the variables will appear as
                                           # self.var[0], self.var[1], ..., self.var[n - 1] instead of the standard self.var + '0', ..., self.var + str(n - 1)
        self.const = const            # Optional argument which keeps track of the formal power series's constant term
        self.children = children      # (Will only have a value other than None on internal calls of the constructor) Points to other FPS objects, which together
                                           # with self.operation define the formal power series under consideration
        self.operation = operation    # (Will only have a value other than None on internal calls of the constructor) self.operation is a list, the first member
                                           # of which is a string denoting the operation and which any additional members of are extra necessary data, i.e. the
                                           # index of the variable with respect to which one is differentiating

    #Overide the repr function
    def __repr__(self):
        return "Formal Power Series in " + str(self.n) + " variables over " + repr(self.ring)

    def __add__(self, ps): # addition operator
        if self.ring != ps.ring:
            raise ValueError("Different Rings")
        elif self.n != ps.n:
            raise ValueError("Incompatible Dimensions")
        elif self.var != ps.var:
            raise ValueError("Different variables")
        else:
            return FPS(self.ring, self.n, None, None, self.var, None, [self,ps], ["+"])

    def __sub__(self, ps): # subtraction operator
        if self.ring != ps.ring:
            raise ValueError("Different Rings")
        elif self.n != ps.n:
            raise ValueError("Incompatible Dimensions")
        elif self.var != ps.var:
            raise ValueError("Different variables")
        else:
            return FPS(self.ring, self.n, None, None, self.var, None, [self,ps], ["-"])

    def __mul__(self, ps): # multiplication operator
        if self.ring != ps.ring:
            raise ValueError("Different Rings")
        elif self.n != ps.n:
            raise ValueError("Incompatible Dimensions")
        elif self.var != ps.var:
            raise ValueError("Different variables")
        else:
            return FPS(self.ring, self.n, None, None, self.var, None, [self,ps], ["fast*"])

    def __div__(self, ps): # division operator
        return (self * ps.inverse())

    #Returns self raised to the power e
    def __pow__(self, e):
        
        if e < 0:
            pow = self.inverse()
            if e < -1:
                pow = pow**(-e)
            return pow
        return FPS(self.ring, self.n, None, None, self.var, None, [self], ["pow", e])

    def inverse(self): # multiplicative inverse
        try:
            unit = self.getConst().is_unit()
        except NotImplementedError:
            unit = True
        if unit:
            return FPS(self.ring, self.n, None, None, self.var, None, [self], ["inv"])
        else:
            raise ValueError("Constant term is not a unit")

    def reversion(self): # compositional inverse
        if self.n != 1:
            raise ValueError("More than one variable")
        poly = self.view(1)
        if poly[0] != 0:
            raise ValueError("Constant term is nonzero")
        # should also check that poly[1] is a unit
        return FPS(self.ring, self.n, None, None, self.var, None, [self], ["rev"])

    def derive(self, i): # differentiate with respect the the ith variable
        # i is the index of the variable with respect to which self is being differentiated, 0 <= i < n because this is how sage will enumerate the variables
        return FPS(self.ring, self.n, None, None, self.var, None, [self], ["d/dx", i])

    def include(self, g, m, newVar = None):
        # g is a list representing a function which takes self.children[0].n arguments to range((new self.n))
        # m is the new self.n
        # the optional newVar argument renames the power series variables
        if newVar == None:
            newVar = self.var
        return FPS(self.ring, m, None, None, newVar, None, [self], ["inc", g])

    def change_base(self, hom):
        return FPS(hom.codomain, self.n, None, None, self.var, None, [self], ["cb", hom])

    # composition of power series
    def __call__(self, seriesList, positions = None):
        # self is the outer (composing) function; seriesList is a list containing the power series that are being used as arguments of self
        # positions is a list of integers which must have an entry for each variable of each member of seriesList. For example, if g and h are
        # power series in one variable, and we are computing f([g, h]), then positions = [0, 1] yields an answer semantically equivalent to
        # f(g(x), h(y)) and if positions = [0, 0] the result is semantically equivalent to f(g(x), h(x)). The default is that positions is an
        # increasing list [0, 1, 2, ... n-1], where n is the total number of variables among every member of seriesList; this corresponds to
        # every variable being independent.
        # For a more complicated example, if g1 and g2 are in two variables and g3 is in three, then changing the position argument to
        # [0,1,1,0,0,1,2] corresponds to f(g1(x0,x1), g2(x1,x0), g3(x0,x1,x2)) whereas
        # [0,2,0,1,0,2,0] corresponds to f(g1(x0,x2), g2(x0,x1), g3(x0,x2,x0))
        # The result will be considered as living under the smallest number of dimensions required by the position argument;
        # returning to the first example, f(g(x), h(x)) is considered to have one variable, if we change the position to [1, 1] then we have
        # an answer corresponding to f(g(y), h(y)); this lives in two variables.
        # Variables can be added and moved around in the include method.
        sumdim = sum([seriesList[i].n for i in range(len(seriesList))])
        if positions == None: # default case is that each of the arguments' variables are independent
            positions = list(range(sumdim))
        if len(seriesList) != self.n:
            raise ValueError("Improper number of arguments")
        if len(positions) != sumdim:
            raise ValueError("Improper number of positions given")
        for ps in seriesList:
            if self.ring != ps.ring:
                raise ValueError("Different Rings")
            if ps.getConst() != 0:
                raise ValueError("Argument has non-zero constant coefficient")
        newdim = max(positions) + 1
        newTrees = []
        divider = 0
        # everything is moved into the appropriate space
        for i in range(len(seriesList)):
            divider2 = divider + (seriesList[i]).n
            # we need to extract the right part of the list positions for each of the power series in seriesList
            newTrees.append(seriesList[i].include(positions[divider:divider2], newdim))
            divider = divider2
        return FPS(self.ring, newdim, None, None, self.var, None, [self] + newTrees, ["o"])

    #Slow multiplication of formal power series
    #The original implementation of multiplication for this package used the multiplication function for the already
    #existent multivariable power series package for Sage.
    #However, this implementation performed several unecessary computations.
    #This function has been kept in the case the user wants to use the older version of multiplication
    def slow_mult(self, ps):
        if self.ring != ps.ring:
            raise ValueError("Different Rings")
        elif self.n != ps.n:
            raise ValueError("Incompatible Dimensions")
        elif self.var != ps.var:
            raise ValueError("Different variables")
        else:
            return FPS(self.ring, self.n, None, None, self.var, None, [self,ps], ["slow*"])

    #Faster composition with outside power series of degree 1 using fast multiplication
    #Will hopefully be updated in a future iteration for a more general case
    def fast_comp(self, ps):
        if self.ring != ps.ring:
            raise ValueError("Different Rings")
        elif self.n != 1:
            raise ValueError("Incompatible Dimensions")
        elif self.var != ps.var:
            raise ValueError("Different variables")
        else:
            return FPS(self.ring, self.n, None, None, self.var, None, [self,ps], ["fast_o", 0])


    #Return constant of the formal power series
    def getConst(self):
        if self.const == None:
            if self.operation == None:
                self.const = self.coeffs([0]*self.n)
            elif self.operation[0] == '+':
                self.const = self.children[0].getConst() + self.children[1].getConst()
            elif self.operation[0] == '-':
                self.const = self.children[0].getConst() - self.children[1].getConst()
            elif self.operation[0] == '*':
                self.const = self.children[0].getConst() * self.children[1].getConst()
            elif self.operation[0] == 'inv':
                self.const = self.children[0].getConst() **(-1)
            elif self.operation[0] == 'rev':
                self.const = 0
            elif self.operation[0] == 'd/dx':
                # TODO
                self.const = self.view(0) # might need a typecast
            elif self.operation[0] == 'inc':
                self.const = self.children[0].getConst()
            elif self.operation[0] == 'cb':
                pass # TODO
            elif self.operation[0] == 'o':
                self.const = self.children[0].getConst() # all the other children should have const == 0
        return self.const

    #Depth-first Search to find the maximal index needed for the base IndRing for a given precision based on data on nodes, and the locations of derivatives and base changes
    #Derivative is the only implemented operation which requires a different level of precision for the children to get the correct precision for the parent
    #Change base has some weirdness due to the rstep and sstep functions if they are not the identity.
    def ring_size(self, prec, derivative_counter):
        #Is a leaf in this case
        if self.children == None:
            return self.b_func(prec + derivative_counter)
        #Is an operation
        else:
            #derivative_counter keeps track of how many derivatives been taken
            new_counter = derivative_counter
            if self.operation[0] == "d/dx":
                new_counter += 1
            #Check children and take max of all of them
            t = 0
            for child in self.children:
                child_size = child.ring_size(prec, new_counter)
                t = max(child_size, t)
            #Ring size might need to be adjusted if we pass any change base operations
            #Need to go through the domain IndRing until we find an index where there is a homomorphism whose range is passed the index of t
            #This will guarantee that we are at least as far out as we need in the IndRing calculations after we apply the change base operation
            if self.operation[0] == "cb":
                if self.operation[1].rstep == universal_identity:
                    t = self.operation[1].sstep(t)
                else:
                    i = 0
                    while self.operation[1].rstep(i) < t:
                        i += 1
                    t = self.operation[1].sstep(i)
            return t
    def view_helper(self, prec, target_ring):
        if self.children == None: # we are at a leaf node with a coefficient function
            H = PowerSeriesRing(self.ring.rings(target_ring), self.n, self.var)
            x = H.gens()
            #Would like to iterate over every non-negative list of length self.n whose components sum to at most prec
            #There are prec + self.n - l choose self.n such lists
            degrees = [0] * self.n
            #Constant term
            poly = self.coeffs(degrees)
            #Initial index of the IndRing for the constant term
            start_ring = self.b_func(0)
            #Push term through to target index in the IndRing
            for i in range(start_ring, target_ring):
                poly = self.ring.homs(i)(poly)
            #Iterate over all degree lists whose sum is less than or equal to the prec
            #Terminates when degree list is all 0's except the last entry which is prec
            while degrees[-1] != prec:
                #Add 1 to the first element in the list until the sum is equal to the prec
                #When the sum is equal to the prec, find the first nonzero term. Change it to 0. 
                #Then increment the following position by 1. Repeat.
                if sum(degrees) == prec:
                    #Find first nonzero coefficient
                    index = 0
                    while degrees[index] == 0:
                        index += 1
                    degrees[index] = 0
                    degrees[index + 1] += 1
                else:
                    degrees[0] += 1
                #Create the monomial term corresponding to the degree list
                product = 1
                #Multiply the powers of x together
                for j in range(self.n):
                    product *= (x[j]**degrees[j])
                #Evaluate the coefficients function to get the coefficient for this monomial
                term = self.coeffs(degrees)
                #Push coefficient to target index in the IndRing
                start_ring = self.b_func(sum(degrees))
                for j in range(start_ring, target_ring):
                    term = (self.ring.homs(j))(term)
                #Add monomial to what has already been constructed
                poly = poly + term*product
            return poly.add_bigoh(prec + 1)
        else: # self is an internal node with an operation
            #Uses existing code for adding multivariable power series in SageMath
            if self.operation[0] == "+":
                return self.children[0].view_helper(prec, target_ring) + self.children[1].view_helper(prec, target_ring)
            #Uses existing code for subtracting multivariable power series in SageMath
            elif self.operation[0] == "-":
                return self.children[0].view_helper(prec, target_ring) - self.children[1].view_helper(prec, target_ring)
            #Uses existing code for multiplying multivariable power series in Sagemath
            #There is also an updated version below
            elif self.operation[0] == "slow*":
                multiplicand = self.children[0].view_helper(prec, target_ring)
                val = multiplicand.valuation()
                return(multiplicand *(self.children[1].view_helper(prec-val, target_ring))).add_bigoh(prec + 1)
            #Faster version of multiplication. 
            #Only computes products of coefficients that are needed to determine the power series up to a certain degree.
            #If possible the algorithm favors the left side of the multiplication having a higher degree.
            elif self.operation[0] == "fast*":
                #Return left side of the multiplication up to degree prec
                multiplicand = self.children[0].view_helper(prec, target_ring)
                #Determine the degree of the first nonzero coefficient
                val = multiplicand.valuation()
                #Dictionary of monomials and their coefficients
                terms = multiplicand.coefficients()
                #Return right side of the multiplication up to degree prec minus the valuation of the left side of the multiplication
                right = self.children[1].view_helper(prec - val, target_ring)
                #Construct product up to degree prec by only multiplying terms whose degree is less than or equal to prec
                #Note right.add_bigoh(prec + 1 - mon.degree()) consists only of terms whose degree when multiplied by the monomial mon 
                #will have degree less than or equal to prec
                ans = 0
                for mon in terms.keys():
                    ans += (terms[mon] * mon * right.add_bigoh(prec + 1 - mon.degree()))
                return ans
            #Works only if the outer function is of degree 1
            #Algorithm is under a similar premise as fast mult
            #Only computes products of coefficients that are needed to determine the power series up to a certain degree.
            elif self.operation[0] == "fast_o":
                outside = self.children[0].view_helper(prec, target_ring)
                base_ring = self.ring.rings(target_ring)
                A = PowerSeriesRing(base_ring, self.var)
                outside_poly = A(outside.polynomial())
                inside_poly = self.children[1].view_helper(prec, target_ring)
                current_power = inside_poly
                total_poly = outside_poly[0]
                if prec == 0:
                    return total_poly
                for i in range(1, prec + 1):
                    total_poly += (outside_poly[i] * current_power)
                    if i < prec:
                        current_power = fast_mult_pow(inside_poly, current_power, prec)
                return total_poly


            elif self.operation[0] == "pow":
                return ((self.children[0].view_helper(prec, target_ring))**self.operation[1]).add_bigoh(prec + 1)
            elif self.operation[0] == "inv":
                return (self.children[0].view_helper(prec, target_ring).inverse().add_bigoh(prec + 1))
            elif self.operation[0] == "rev":
                poly = self.children[0].view_helper(prec, target_ring)
                A = poly.parent()
                B = PowerSeriesRing(self.ring.rings(target_ring), self.var)
                poly2 = B(poly.polynomial())
                new_poly = A((poly2.reverse(prec + 1)).polynomial())
                return new_poly.add_bigoh(prec + 1)
            elif self.operation[0] == "d/dx":
                x = PowerSeriesRing(self.ring.rings(target_ring), self.n, self.var).gens()
                return (self.children[0].view_helper(prec + 1, target_ring).derivative(x[self.operation[1]])).add_bigoh(prec + 1)
            elif self.operation[0] == "cb":
                if(self.operation[1].rstep == universal_identity and self.operation[1].sstep == universal_identity):
                    poly = self.children[0].view_helper(prec, target_ring)
                    mon_dict = poly.coefficients()
                    t = 0
                    for term in mon_dict:
                        t = t + term*(self.operation[1].maps(target_ring)(mon_dict[term]))
                    return t.add_bigoh(prec+1)
                else:
                    #Make polynomial below cb node and get coefficients
                    new_target = self.children[0].ring_size(prec, 0)
                    poly = self.children[0].view_helper(prec, new_target)
                    mon_dict = poly.coefficients()
                    #Obtain information about which rings everything is in so the correct homomorphisms can be applied
                    j = 0
                    while(self.operation[1].rstep(j) < new_target):
                        j += 1
                    final_base = self.operation[1].rstep(j)
                    returned_ring = self.operation[1].sstep(j)
                    #Loop everything together and apply the correct homomorphisms
                    t = 0
                    for term in mon_dict:
                        element = mon_dict[term]
                        #Apply correct number of homomorphisms in domain
                        for i in range(new_target, final_base):
                            element = self.operation[1].domain.homs(i)(element)
                        #Switch Rings
                        element = self.operation[1].maps(j)(element)
                        #Apply correct number of homomorphisms in codomain
                        for i in range(returned_ring, target_ring):
                            element = self.operation[1].codomain.homs(i)(element)
                        t = t + term*element
                    return t.add_bigoh(prec +1)
            elif self.operation[0] == "inc":
                if(self.children[0].n == 1):
                    poly = self.children[0].view_helper(prec, target_ring)
                    H = poly.parent()
                    x = H.gens()
                    y = PowerSeriesRing(self.ring.rings(target_ring), self.n, self.var).gens()
                    poly2 = 0
                    poly_coeffs = poly.coefficients()
                    for i in range(0, prec + 1):
                        if (x[0]**i in poly_coeffs.keys()):
                            poly2 += poly_coeffs[x[0]**i] * y[(self.operation[1])[0]]**i
                    return poly2.add_bigoh(prec + 1)
                K = PowerSeriesRing(self.ring.rings(self.ring_size(prec, 0)), self.n, 'zzz')
                y = K.gens()
                H = PowerSeriesRing(self.ring.rings(self.ring_size(prec, 0)), self.children[0].n, self.var)
                x = H.gens()
                image = []
                for i in range(len(self.operation[1])):
                    image += [y[(self.operation[1])[i]]]
                phi = Hom(H, K)(image)
                poly = phi((self.children[0]).view_helper(prec, target_ring))
                # can discard H now
                H = PowerSeriesRing(self.ring.rings(self.ring_size(prec, 0)), self.n, self.var)
                x = H.gens()
                id = Hom(K, H)(x)
                return id(poly)
            elif self.operation[0] == "o":
                if self.children[0].n == 1:
                    outside_poly = (self.children[0].view_helper(prec, target_ring)).polynomial()
                    old_series = self.children[1].view_helper(prec, target_ring)
                    H = old_series.parent()
                    inside_poly = (old_series.polynomial())
                    comp_poly = outside_poly(inside_poly)
                    new_series = H(comp_poly).add_bigoh(prec + 1)
                    return new_series
                outside_poly = self.children[0].view_helper(prec, target_ring)
                #Can reduce precision needed for the inside power series based on the valuation of the outside power series disregarding the constant term
                val = (outside_poly).valuation() # TODO: make sure cases when val != 1 work
                viewList = []
                for i in range(1, len(self.children)):
                    viewList.append(self.children[i].view_helper(prec - val + 1, target_ring))
                return (outside_poly(*viewList)).add_bigoh(prec + 1)

    def view(self, prec):
        target = self.ring_size(prec, 0)
        return self.view_helper(prec, target)

class UFGL(FPS):
    def __init__(self):
        def v(n):
            fctr = list(factor(n))
            if len(fctr) == 1: # len(fctr) is 0 for n = 1; this checks if n is prime
                return fctr[0][0]
            else:
                return 1
        def c(p, d):
            q = v(d)
            if q == 1 or q == p:
                return 1
            else:
                R = Integers(p)
                ans = (int(R(q) ^-1) * q)
                return ans
        def mu(n, d):
            ans = 1
            fctr = list(factor(n))
            for num in fctr:
                ans *= c(num[0], d)
            return ans
        def log_coeffs(l):
            def b(n, d = None):
                if n == 0:
                    return 0
                if n == 1:
                    return 1
                if d == None:
                    d = n
                R = PolynomialRing(QQ, d, 'U')
                U = R.gens()
                ans = U[n - 1]
                div = (divisors(n))[1:-1]
                for num in div:
                    ans += ((mu(n, num)*v(n) / v(num)) * b(n / num, d) * U[num - 1]**(n/num))
                ans /= v(n)
                return ans
            n = l[0]
            ans = b(n)
            return ans
        ring = Rational_Lazard
        logarithm = FPS(ring, 1, log_coeffs)
        log1 = logarithm.include([0], 2)
        log2 = logarithm.include([1], 2)
        add = log1 + log2
        rev = logarithm.reversion()
        law = rev([add])
        #print("law")
        #print(law.view(4))
        self.ring = ring
        self.n = 2
        self.b_func = None
        self.coeffs = None
        self.children = [rev, add]
        self.operation = ["fast_o", 0]
        self.var = 'x'
        self.const = 0


def Poly_To_FPS(poly, big_ring, size, var_name = 'x'):
    #returns list of generators
    x = (poly.parent()).gens()
    #Figure out how many generators
    n = len(x)
    #Make line of code to find a coefficient given a list of powers of variables l
    code = "poly.coefficient({"
    for i in range(0, n):
        code = code + "x[" +str(i) + "]:l[" + str(i) + "]"
        if i < n-1:
            code = code + ","
    code = code + "})"
    return FPS(big_ring, n, lambda l: eval(code) if max(l) < poly.degree() + 1 else 0, lambda k: size, var_name)




#Example Coefficient function
def f_coeffs(powers):
    n = powers[0] + powers[1]
    if n == 0:
        return 0
    R = PolynomialRing(ZZ, n, 'U')
    U = R.gens()
    return (powers[0]*U[0] + powers[1]*U[-1])

f = FPS(Lazard, 2, f_coeffs, universal_identity,)
#repr(f); f
print("f view 3:")
print(f.view(3))
print("\n")

##Another example
def g_b(size):
    return size**2

def g_coeffs(powers):
    n = (powers[0] + powers[1]) ** 2
    if n == 0:
        return 0
    R = PolynomialRing(ZZ, n, 'U')
    U = R.gens()
    return (powers[0]*powers[1]*U[-1])
g = FPS(Lazard, 2, g_coeffs, g_b)


print("g view 3:")
print(g.view(10))
print("\n")

#Addition
h = g+f
print("g+f view 3:")
print(h.view(3))
print("g+f view 4:")
print(h.view(4))
print("\n")

#Derivatives
h = (g+f).derive(0).derive(0)
print("Two derivatives of g+f with respect to x0:")
print(h.view(3))
print("\n")

#Multipication
h = f*g
print("f*g view 5:")
print(h.view(5))
print("\n")

#Slow Mult
h = f.slow_mult(g)
print("f slow mult g:")
print(h.view(5))
print("\n")
#Example of user defined IndRing
def perfect_closure_ring_7(n):
    prime = 3
    if n == 0:
        return PolynomialRing(GF(prime), 1, 'c0')
    else:
        R = PolynomialRing(GF(prime), 2, "xy")
        x, y = R.gens()
        S = R.quotient(y**(prime**n) - x, names=('c0', 'c' + str(n),))
        return S
def perfect_closure_homs_7(n):
    prime = 3
    if n == 0:
        return perfect_closure_ring_7(0).hom((perfect_closure_ring_7(1).gens())[:-1], perfect_closure_ring_7(1))
    else:
        generators = perfect_closure_ring_7(n+1).gens()
        return perfect_closure_ring_7(n).hom([generators[0],generators[1]**prime], perfect_closure_ring_7(n+1))
perfect_closure_7 = IndRing(perfect_closure_ring_7, perfect_closure_homs_7)

#Example of a power series of this IndRing plus finding its inverse
def k_coeffs(natural):
    if natural[0] == 0:
        R = perfect_closure_ring_7(0)
        c = R.gens()
        return(1 + 0*c[0])
    else:
        R = perfect_closure_ring_7(natural[0] - 1)
        c = R.gens()
        return(1*c[-1])

def k_b (natural):
    if natural == 0:
        return 0
    else:
        return natural - 1

def l_coeffs(natural):
    prime = 3
    R = perfect_closure_ring_7(natural[0])
    c = R.gens()
    total = 0
    for i in range(0, natural[0] + 1):
        total += (c[-1]**(prime**(i)))
    return total
def l_b(natural):
    return natural

print("Inverses and perfect closure ring:")
k = FPS(perfect_closure_7, 1, k_coeffs, k_b)
l = FPS(perfect_closure_7, 1, l_coeffs, l_b)
print("k view 5:")
print(k.view(5))
print("")
print("k inverse view 5:")
print((k.inverse()).view(5))
print("l view 5:")
print(l.view(5))
print("k * l view 5:")
print((k*l).view(5))
print("\n")


#Reversion
print("Reversion")
def p_func(integer):
    return integer[0]
p = FPS(Integer_Ring, 1,  p_func, universal_identity)
print("p view 3:")
print(p.view(3))
print("reverse of p view 3:")
print(p.reversion().view(3))
print("\n")
#Change Base Examples
def lazard_to_integers_f(n):
    ring = PolynomialRing(ZZ, n, 'U')
    eval_list = [1]*(n)
    hom = ring.hom(eval_list, ZZ)
    return hom

print("H is a homomorphism that sends all of the ci's in the Lazard ring to 1")
H = IndRingHom(Lazard, Integer_Ring, lazard_to_integers_f)
h = f.change_base(H)
print(" ")
print("H(f): ")
print(h.view(3))

j = g.change_base(H)
print("H(g):")
print(j.view(3))

def lazard_to_integers_g(n):
    ring = PolynomialRing(ZZ, n**2, 'c')
    eval_list = []
    for i in range(0, n**2):
        eval_list.append(i)
    hom = ring.hom(eval_list, ZZ)
    return hom
K = IndRingHom(Lazard, Integer_Ring, lazard_to_integers_g, lambda n: n**2, lambda n: 2*n)
j = g.change_base(K)
print("K(g):")
print(j.view(10))
print("")


print("Inclusion: variable swap of f; third variable added")
j = f.include([1, 0], 3)
print(j.view(4))


def f_coeffs(powers):
    n = powers[0]
    if n == 0:
        return 0
    R = PolynomialRing(ZZ, n, 'U')
    U = R.gens()
    return (powers[0]*U[0] + U[-1])

f = FPS(Lazard, 1, f_coeffs, universal_identity)
print("")
print("one-dimensional f")
print(f.view(3))
k = f.include([1], 2)
print("f in two dimensions")
print(k.view(3))
print("f([g]) with f in two dimensions")
print(k([g, g]).view(4))
print("One variable composition")
print(f([g]).view(4))


#Turn Multivariable Power Series into FPS data type
#Make example multivariable power series
Lazard_3 = Lazard.rings(3)
U = Lazard_3.gens()
L = PolynomialRing(Lazard_3, 2, 'xy')
x = L.gens()
poly = 0

for i in range(0, 3):
    for j in range(0, 3):
        poly += (i+j+1)*U[0]*U[1]*U[i]^j*x[0]^i*x[1]^j

print("poly:")
print(poly)

#Turn this example into an FPS data type
q = Poly_To_FPS(poly, Lazard, 3, 'xy')
print("q:")
print(q.view(4))
print("\n")


#Test how sage treats rationals to zmodp
print("This is a test.")
R = PolynomialRing(QQ, 2, 'xy')
x = R.gens()
poly2 = 7*x[0]*x[1] + 8*x[0]^3 + 9
print("poly2:")
print(poly2)
q = Poly_To_FPS(poly2, Rational_Ring, 3, 'xy')

def integers_to_modp(natural):
    S = ZZ.quo(7*ZZ)
    pi = S.cover()
    return pi

K = IndRingHom(Rationals, Integer_7_Ring, integers_to_modp)
print("q:")
print(q.view(5))
print("t:")
t = q.change_base(K)
print(t.view(5))
print("\n")
      
#Universal Formal Group Law
ufgl = UFGL()
print("UFGL")
print(ufgl.view(4))
print("\n\n")
#print(ufgl.view(3).parent())
