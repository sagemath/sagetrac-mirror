'''

Examples

You can import this from python to see some examples run.

'''


from sage.rings.integer_ring import *
from sage.rings.rational_field import *
from sage.rings.finite_rings.integer_mod_ring import *
from sage.categories.homset import *
from sage.rings.power_series_ring import *
from sage.rings.finite_rings.finite_field_constructor import *
from sage.arith.misc import *

from formal_power_series import *

print('1/(1-x) as a formal power series:')

def spacing():
	print('\n'+'#'*80 + '\n')

geometric = FPS(Rational_Ring, 1, lambda n: Integer(1))
print(geometric.view(5))

print('Inverse of the previous formal power series, up to x^60:')

print(geometric.inverse().view(60))

spacing()

print('Consider the example formal power series g and f.')

#Example Coefficient function
def f_coeffs(powers):
    n = powers[0] + powers[1]
    if n == 0:
        return 0
    R = PolynomialRing(ZZ, n, 'U')
    U = R.gens()
    return (powers[0]*U[0] + powers[1]*U[-1])

# constructs an element of Lazard[[x, y]] with coefficients given by f_coeffs
f = FPS(Lazard, 2, f_coeffs, universal_identity,)
#repr(f); f
print("first 3 coefficients of f:")
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

print("first 3 coefficients of g:")
print(g.view(3))
print("\n")

#Addition
h = g+f
print("g+f up to 3rd coefficient:")
print(h.view(3))
print("\n")

#Derivatives
h = (g+f).derive(0).derive(0)
print("Second derivatives of g+f with respect to x0:")
print(h.view(3))
print("\n")

#Multipication
h = f*g
print("f*g up to 5:")
print(h.view(5))
print("\n")

##Slow Mult
#h = f.slow_mult(g)
#print("f slow mult g:")
#print(h.view(5))
#print("\n")

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

spacing()

print('''
Inverse example: perfect closure of a ring.

k is GF(7). R is the colimit along k[c_0, c_i]/(c_i**(p**i) - c_0) --> k[c_0, c_{i+1}]/(c_{i+1}**(p**i) - c_0)
sending c_i to c_{i+1}^p and c_0 to c_0.

This demonstrates you can take inverses over IndRings.
''')

f_perfect = FPS(perfect_closure_7, 1, k_coeffs, k_b)
l = FPS(perfect_closure_7, 1, l_coeffs, l_b)
print("First 5 coefficients of f:")
print(f_perfect.view(5))
print('')
print("f inverse view 5:")
print((f_perfect.inverse()).view(5))
print("l up to 5:")
print(l.view(5))
print("f * l up to 5:")
print((f_perfect*l).view(5))
print("\n")


#Reversion
spacing()
print("Reversion example")
def p_func(integer):
    return integer[0]
p = FPS(Integer_Ring, 1,  p_func, universal_identity)
print("first 3 coefficients of p:")
print(p.view(3))
print("reversopn of p up to coefficient 3:")
print(p.reversion().view(3))
print("\n")
#Change Base Examples
def lazard_to_integers_f(n):
    ring = PolynomialRing(ZZ, n, 'U')
    eval_list = [1]*(n)
    hom = ring.hom(eval_list, ZZ)
    return hom

spacing()

print('You can also create maps of IndRings and push coefficients along them.')
print("H is a homomorphism that sends all of the c_i's in the Lazard ring to 1")
print('Example of applying H to f and g:')
H = IndRingHom(Lazard, Integer_Ring, lazard_to_integers_f)
h = f.change_base(H)
print(" ")
print("H(f): ")
print(h.view(3))

j = g.change_base(H)
print("H(g):")
print(j.view(3))

print('\nK maps each U_i to i where U_i is the polynomial generator of the Lazard ring.')

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

spacing()

print("Inclusion example:")
print("swapping variables in generators of polynomial ring and add a new one (x2):")
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
        poly += (i+j+1)*U[0]*U[1]*U[i]**j*x[0]**i*x[1]**j

print("poly:")
print(poly)

#Turn this example into an FPS data type
q = Poly_To_FPS(poly, Lazard, 3, 'xy')
print("q:")
print(q.view(4))
print("\n")


##Test how sage treats rationals to zmodp
#print("This is a test.")
#R = PolynomialRing(QQ, 2, 'xy')
#x = R.gens()
#poly2 = 7*x[0]*x[1] + 8*x[0]**3 + 9
#print("poly2:")
#print(poly2)
#q = Poly_To_FPS(poly2, Rational_Ring, 3, 'xy')
#
#def integers_to_modp(natural):
#    S = ZZ.quo(7*ZZ)
#    pi = S.cover()
#    return pi
#
#K = IndRingHom(QQ, Integer_7_Ring, integers_to_modp)
#print("q:")
#print(q.view(5))
#print("t:")
#t = q.change_base(K)
#print(t.view(5))
#print("\n")

spacing()
      
#Universal Formal Group Law
ufgl = UFGL()
print("Universal formal group law (first 7 coefficients):")
print(ufgl.view(7))
