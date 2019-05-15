from formal_power_series import *
import indring


# F(X, Y) = F(Y, X) ## so we only need half of the coefficients
# F(X, 0) = X
# F(X, F(Y, Z)) = F(F(X, Y), Z)
# implies XY divides (F-X-Y)

# want: base change
# homomorphisms: f(F(X, Y)) = G(f(X),f(Y)) where F,G are formal group laws
class FormalGroupLaw(FPS):
    # check that there are two variables, constant coefficient is zero, and degree 1 coefficients are 1
    # user can pass either a coefficient and IndRing depth function, or an FPS object
    def __init__(self, ring, coeffs, b_func = lambda n: n, var='XY', children = None, operation = None):
    # Python does not allow overloading, so the constructor has to be prepared for any arguments FPS() might take
    # except for n, which is always equal to 2, and const, which is always equal to 0
        if coeffs != None:
            if not (coeffs([0,1]) == 1 and coeffs([1,0]) == 1):
                raise ValueError("Improper coefficient function")
            elif coeffs([0, 0]) != 0:
                raise ValueError("Nonzero constant term")
        self.ring = ring
        self.n = 2
        self.b_func = b_func
        self.coeffs = coeffs
        self.children = children
        self.operation = operation
        self.var = var
	self.memoized_degree = None
	self.memoized_poly = None

    def __repr__(self):
        return "some Formal Group Law"

    # i-series, aka n-series
    # returns a power series in one variable
    def iseries(self, i): # having X as an optional argument avoids it being declared as a variable in every call
        # returns f composed with itself to the power of n
        # f(x, f(x, f(x, f(x, ... f(x, x)))))
        # f.series(1) = f(X, X)
        # f(x, f.iseries(-1)) == 1
	zero = FPS(self.ring, 1, lambda n: 0, self.b_func)
        def twovar(i, poly):
            if i < 0:
                # maybe define this someday, inverse of x in a FGL
                raise ValueError("i must be positive")
            elif i == 0:
                return zero
            elif i == 1:
                return poly
            else:
                return self([poly, twovar(i-1, poly)]).include([0,0], 1, "x")
        Xcoeffs = lambda l: 1 if l == [1] else 0
        x_lin = FPS(self.ring, 1, Xcoeffs, lambda l: 1)
        ps = twovar(i, x_lin)
        return ps
    # conjugate takes some f in R[[z]] and a FGL G(x, y) in R[[x,y]] and conjugates G wrt f:
    # G.conjugate(f) = f^(-1) (G(f(x), f(y)))
    # where f^(-1) is the compositional inverse (reversion).
    # the value returned is a FGL
    def conjugate(self, f): # f should be in one variable
        frev = f.reversion()
        # intermediate: self(f(x), f(y))
        intermediate = self.__call__(seriesList = [f, f])
        #test = grouplaw(intermediate)
        ps = frev([intermediate])
        ps.view(4)
        return grouplaw(ps)

    def forget(self):
        return FPS(self.ring, self.n, self.b_func, self.coeffs, self.var, self.children, self.operation)

    def log(self):
        # logarithm: find  (d/dy F(x, y)) evaluated at y = 0 and return the integral of this dx.

        # differentiate with respect to y
        derivative = self.derive(1)
        # evaluate by setting y = 0
        
        evaluated = derivative([linear(self.ring), zero(self.ring, 2, "xy")], [0, 0, 1, 1])
        # multiplicative inverse
        inverted = evaluated.inverse()
        # integrate with respect to x
        integrated = inverted.integrate(0)
        return integrated

def zero (ring, n, var):
    return FPS(ring, n, lambda l: 0, lambda n: 1, var)

def linear(ring):
    return FPS(ring, 2,lambda l: 1 if l == [1,0] else 0, lambda n: 1, "xy")

def ptypical(p, gen = "hazelwinkel"):
    if gen == "hazelwinkel":
        def L(n):
            if n == 0:
                return 1
            else:
                nu = PolynomialRing(QQ, n+1, 'V').gens()
                prev = [1]
                for j in range(1, n+1):
                    nx = nu[j+1]
                    for k in range(j):
                        nx += prev[k] * v[n+1-k]**(p**k)
                    nx /= p
                    prev.append(nx)
                return prev[-1]
        def logcoeffs(l):
            n = l[0]
            k = 0
            while n % p == 0:
                n //= p
                k += 1
            if n == 1:
                return L(k)
            else:
                return 0
        def discretelog(n):
            x = 0
            while n > 1:
                x += 1
                n //= p
            return x
        lazard = IndRing('lazard')
        logr = FPS(ring = lazard, n = 1, coeffs = logcoeffs, b_func = discretelog, var = 'XY')
        addition = FormalGroupLaw(ring = lazard, coeffs = lambda l: int(l == [0,1] or l == [1,0]), b_func = lambda n: 0)
        return addition.conjugate(logr)
    elif gen == 'araki':
        pass


def grouplaw(ps):
    if ps.n != 2:
        raise ValueError("Wrong number of variables for FormalGroupLaw: " + str(ps.n))
    
    # check that ps.view(1) == x+y
    poly = ps.view(1).polynomial()
    A = poly.parent()
    x = A.gens()
    if (poly.coefficient({x[0]:1, x[1]:0}) != 1) or (poly.coefficient({x[0]:1, x[1]:0}) != 1):
        raise ValueError("Not a formal group law")
    return FormalGroupLaw(ps.ring, ps.coeffs, ps.b_func, ps.var, ps.children, ps.operation)

#def forget(fgl):
#    return FPS(fgl.ring, fgl.n, fgl.b_func, fgl.coeffs, fgl.var, fgl.const, fgl.children, fgl.operation)

# print("debug1")
# hz = ptypical(p=2)
# hz.view(3)

rationals = IndRing(QQ)
#basic = FormalGroupLaw(rationals, lambda l: int(l in [[0,1], [1,0], [1,1]]))
#basic.view(4)
#
#ps = FPS(rationals, 1, lambda l: int(l[0] in [1, 2]), lambda n: n)
#ps.view(4)
#
#rev = ps.reversion()
#print("rev", rev.view(4))
#intermediate = basic([ps, ps])
#print("intermediate", intermediate.view(4))
#after_conj = rev([intermediate])
#print("after_conj", after_conj.view(4))
#conj = basic.conjugate(ps)
#conj.view(4)
#
#
#print('it has run.')

#Example Coefficient function
def universal_identity(natural):
    return natural
Lazard = IndRing('lazard')
Integer_Ring = IndRing(ZZ)
#def f_coeffs(powers):
#    n = powers[0] + powers[1]
#    R = PolynomialRing(ZZ, n, 'U')
#    U = R.gens()
#    if n == 0:
#        return 2
#    return (powers[0]*U[0] + powers[1]*U[-1])
#
#f = FPS(Lazard, 2, f_coeffs, universal_identity)
#print("f:",f.view(3))
#def lazard_to_integers_f(n):
#    ring = PolynomialRing(ZZ, n, 'U')
#    eval_list = [2]*(n)
#    hom = ring.hom(eval_list, ZZ)
#    return hom
#
#print("H is a homomorphism that sends all of the ci's in the Lazard ring to 1")
#H = IndRingHom(Lazard, Integer_Ring, lazard_to_integers_f)
#h = f.change_base(H)
#print("h", h.view(3))
