# This file is tested with SageMath version 9.4
# open a sage session and type
# attach(".../K3SurfaceAutGrp.py")
# with the 3 dots replaced by the path to this file
# to load the database


from sage.rings.all import PolynomialRing,ZZ
from sage.misc.functional import cyclotomic_polynomial
from sage.functions.other import sqrt
from sage.misc.sage_eval import sage_eval
from sage.modules.all import IntegralLattice
from sage.matrix.all import matrix
from sage.arith.misc import euler_phi

script_path = ""
for f in sage.repl.attach.attached_files():
    if f[-18:] == "K3SurfaceAutGrp.py":
      script_path = f[:-18]
      break
else:
    print("use attach to load K3SurfaceAutGrp.py")

script_path = "k3s/new/"
def kernel_sublattice(L, f):
    r"""
    Return the the kernel of `f` as a sublattice of `L`.

    Input:

    - ``L`` -- an integral lattice
    - ``f`` -- a matrix of an element of the orthogonal group `O(L)`

    Examples::

        sage: f = matrix(ZZ, 2, 2, [0,1,1,0])
        sage: L = IntegralLattice("A2")
        sage: kernel_sublattice(L, f - 1)
        Lattice of degree 2 and rank 1 over Integer Ring
        Basis matrix:
        [1 1]
        Inner product matrix:
        [ 2 -1]
        [-1  2]
    """
    try:
        f = f.matrix()
    except AttributeError:
        pass
    K = L.ambient_vector_space().span(f.kernel().gens())
    return L.sublattice((L & K).gens())


def invariant_lattice(L, G):
    r"""
    Return the invariant lattice `L^G`.

    Input:

    - ``L`` -- a lattice
    - ``G`` -- either a list of matrices or
               a subgroup of the orthogonal group `O(L)` of `L`

        sage: L = IntegralLattice("A2")
        sage: G = L.orthogonal_group()
        sage: S = G.subgroup(G.gens()[:2])
        sage: invariant_lattice(L, S)
        Lattice of degree 2 and rank 0 over Integer Ring
        Basis matrix:
        []
        Inner product matrix:
        [ 2 -1]
        [-1  2]
    """
    F = L
    try:
        gens = [g.matrix() for g in G.gens()]
    except AttributeError:
        gens = G
    for g in gens:
        F = kernel_sublattice(F, g - 1)
    return F


def maximum(self):
    r"""
    Return the maximum of this lattice.

    EXAMPLES::

        sage: L = IntegralLattice('A2')
        sage: maximum(L)
        +Infinity
        sage: maximum(L.twist(-1))
        -2
    """
    p, n = self.signature_pair()
    if p != 0:
        from sage.rings.infinity import PlusInfinity
        return PlusInfinity()
    mpari = (-self.gram_matrix()).__pari__().qfminim(None,0)[1]
    return -mpari

def small_local_type(L, f, n=None):
    r"""
    Return the small local type of the isometry `f` of `L`.

    Input:

    - ``L`` -- an integral lattice
    - ``f`` -- an isometry of `L`

    Output:

    A tuple `((1,G(L_1)),...,(n,G(L_n))` with
    `L_i = L^{f^i}` and `G(N)` denoting the genus of `N`.
    """
    if n is None:
        n = f.change_ring(ZZ).multiplicative_order()
    D = n.divisors()
    return tuple((i, kernel_sublattice(L,f**i-1)) for i in D)



class K3SurfaceAutGrp(object):
    r"""
    This class represents a finite group of automorphisms of
    a complex K3 surface with a non-sympletic action.

    Input:

    - ``L`` -- an even unimodular lattice of signature (3,19)
    - ``H0`` -- a list of matrices in the orthogonal group `O(L)`
      inducing symplectic automorphisms, that is,
      the coinvariant lattice L_H0 is negative definite and root free
    - ``g`` -- a matrix representing an isometry of `L` such that `g` respects the
      orthogonal decomposition `L = `L_H0` + `L^H0` and `g` is non-symplectic
    - ``n`` -- the transcendental value, i.e. the order of `<H0,g> / H0`
    - ``name`` -- (default: ``None``) optional name used for printing
    - ``check`` -- (default: ``True``) check if the input actually defines a
      group of automorphisms of some complex K3 surface.
    - ``components`` -- (default: ``None``) if ``n==2``, the number of connected
      components of the moduli space of `H`-marked K3 surfaces
    """

    def __init__(self, L, H0, g, n, name=None, check=True, components=None, properties=[]):
        self._id = name
        self._L = L
        self._H0 = H0
        g.set_immutable()
        self._g = g
        self._H = H0 + [g]
        self._n = n
        self._prop = dict(properties)
        if self._n != 2:
          components = ZZ(1)
        self._components = ZZ(components)
        if check:
          if not (L.determinant()==-1 and L.is_even() and L.signature_pair() == (3,19)):
              raise ValueError('L must be an even unimodular lattice '
                                + 'of signature (3,19)')
          if maximum(self.NS_coinvariant()) == -2:
              raise ValueError('Does not perserve an ample class')
          if self.H().order()/self.H0().order() != n:
              raise ValueError('n does not match the order of G/H0')
          if self.invariant_lattice().signature_pair()[0] == 0:
              raise ValueError('Does not preserve an ample class')
          if self.transcendental_lattice().signature_pair()[0] !=2:
              raise ValueError('g is not non-sympletic')

    def __repr__(self):
      if self._id is not None:
        return f"K3Group {self._id}"
      s = f"K3Group of order {self.H().order()} = {self.H0().order()} * {self.n()}"
      return s

    def distinguished_generator(self):
        return self._g

    def id(self):
        return self._id

    def properties(self):
      return self._prop
    end

    def connected_components(self):
        r"""
        Return the number of connected components of the moduli
        space of `(H,\chi)`-marked K3 surfaces.

        If the transcendental value is `>2` this is automatically one.
        """
        if self.n() > 2:
            return 1
        return self._components

    def L(self):
        r"""
        Return the K3 lattice of `self`.
        """
        return self._L

    k3_lattice = L

    def transcendental_value(self):
        r"""
        Return the order of `H/H0`.
        """
        return self._n

    n = transcendental_value

    @cached_method
    def H0(self):
        """
        Return the symplectic subgroup of `H`.
        """
        return self.L().orthogonal_group(self._H0)

    symplectic_subgroup = H0

    @cached_method
    def head(self):
        F = self.symplectic_invariant_lattice()
        C = self.heart()
        f = self.distinguished_generator()
        # restrict f to F
        # then continue it as the identity on F^\perp
        assert F.degree()==22
        B = F.basis_matrix().stack(C.basis_matrix())
        r = F.rank()
        I = matrix.identity(22-r)
        fB = B *f * B.inverse()
        fB[r:,r:] = I
        f = B.inverse() * fB * B
        assert f*F.inner_product_matrix()*f.T == F.inner_product_matrix()
        assert all(b*f in F for b in F.basis())
        f.set_immutable()
        return (F,f)

    @cached_method
    def H(self):
        r"""
        Return `H`.

        Examples::

            sage: k3 = K3SurfaceAutGrp_from_database(30,False)[0]
            sage: k3.H()
            Group of isometries with 2 generators (
            [ 1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0]
            [ 0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0]
            [ 0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0]
            [ 0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0]
            [ 0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0]
            [ 0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0]
            [ 0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0]
            [ 0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0]
            [ 0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0]
            [ 0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0]
            [ 0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0]
            [ 0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0]
            [ 0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0]
            [ 0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0]
            [ 0  0  0  0  0  0  0  0  0  0  0  0  0  0 -1  0  0  0  0  0  0  0]
            [ 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 -1  0  0  0  0  0  0]
            [ 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 -1  0  0  0  0  0]
            [ 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 -1  0  0  0  0]
            [ 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 -1  0  0  0]
            [ 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 -1  0  0]
            [ 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 -1  0]
            [ 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 -1],
            <BLANKLINE>
            [ 0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0]
            [ 0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0]
            [ 0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0]
            [ 0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0]
            [ 0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0]
            [ 0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0]
            [ 0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0]
            [-1 -1  0  1  1  1  0 -1  0  0  0  0  0  0  0  0  0  0  0  0  0  0]
            [ 0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0]
            [ 0  0  0  0  0  0  0  0 -1  1  0  0  0  0  0  0  0  0  0  0  0  0]
            [ 0  0  0  0  0  0  0  0  0  0 -1  0  0  0  0  0  0  0  0  0  0  0]
            [ 0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0]
            [ 0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0]
            [ 0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0]
            [ 0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  1  0  1  0]
            [ 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 -1  0]
            [ 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 -1  0  0 -1  0 -1]
            [ 0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  1  0  0 -1]
            [ 0  0  0  0  0  0  0  0  0  0  0  0  0  0 -1 -1  1  1  0  1  1  1]
            [ 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0 -1 -1 -1 -1  0]
            [ 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 -1  0 -1  0 -1  0]
            [ 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 -1  1  0  1  1  1  1]
            )
            """
        return self.L().orthogonal_group(self._H)

    group = H

    @cached_method
    def transcendental_lattice(self):
        r"""
        Return the transcendental lattice of a very general K3 surface admitting
        the given group action.

        Examples::

            sage: k3 = K3SurfaceAutGrp_from_database(30,False)[0]
            sage: k3.transcendental_lattice().genus()
            Genus of
            [-4  4 -4  0  2 -2  2 -4]
            [ 4 -4  4 -4  0  2 -2  2]
            [-4  4 -4  4 -4  0  2 -2]
            [ 0 -4  4 -4  4 -4  0  2]
            [ 2  0 -4  4 -4  4 -4  0]
            [-2  2  0 -4  4 -4  4 -4]
            [ 2 -2  2  0 -4  4 -4  4]
            [-4  2 -2  2  0 -4  4 -4]
            Signature:  (2, 6)
            Genus symbol at 2:    2^8
            Genus symbol at 5:     1^-6 5^-2
        """
        R = PolynomialRing(QQ,"x")
        x = R.gen()
        cn = cyclotomic_polynomial(self._n,'x')
        g = self._g
        S = self.symplectic_invariant_lattice()
        K = kernel_sublattice(self._L, cn(g))
        return self._L.sublattice((K & S).gens())

    T = transcendental_lattice

    @cached_method
    def NS(self):
        r"""
        Return the generic Neron-Severi lattice of a K3 surface admitting
        the given group action.

        Examples::

            sage: k3 = K3SurfaceAutGrp_from_database(30,False)[0]
            sage: k3.NS().genus()
            Genus of
            [ -4  -5  -1   0  -1   0   0   0   0   0   0   0   0   0]
            [ -5 -10   0   0   0   0   0   0   0   0   0   0   0   0]
            [ -1   0  -2   0  -1   0   0   0   0   0   0   0   0   0]
            [  0   0   0   0   0   1   0   0   0   0   0   0   0   0]
            [ -1   0  -1   0  -2   0   0   0   0   0   0   0   0   0]
            [  0   0   0   1   0   0   0   0   0   0   0   0   0   0]
            [  0   0   0   0   0   0  -4   2  -2  -2   2   2   2  -2]
            [  0   0   0   0   0   0   2  -4   2   0   0  -2  -2   0]
            [  0   0   0   0   0   0  -2   2  -4  -2   2   2   2   0]
            [  0   0   0   0   0   0  -2   0  -2  -4   2   2   2  -2]
            [  0   0   0   0   0   0   2   0   2   2  -4   0   0   0]
            [  0   0   0   0   0   0   2  -2   2   2   0  -4  -2   2]
            [  0   0   0   0   0   0   2  -2   2   2   0  -2  -4   2]
            [  0   0   0   0   0   0  -2   0   0  -2   0   2   2  -4]
            Signature:  (1, 13)
            Genus symbol at 2:    1^6:2^8
            Genus symbol at 5:     1^-12 5^-2
        """
        T = self.transcendental_lattice()
        return self._L.orthogonal_complement(T)

    @cached_method
    def symplectic_invariant_lattice(self):
        r"""
        Return the invariant lattice of the subgroup of symplectic automorphisms.

        Examples::

            sage: k3 = K3SurfaceAutGrp_from_database(30,False)[0]
            sage: k3.symplectic_invariant_lattice().genus()
            Genus of
            [ -6  -7  -2   2  -2  -2   2  -2  -4  -6  -1   0  -1   0]
            [ -7 -12  -2   2  -4   2  -2   0  -5  -9  -1   0   0   0]
            [ -2  -2  -4   4  -4   0   2  -2   0   0   0   0   0   0]
            [  2   2   4  -4   4  -4   0   2   0   0   0   0   0   0]
            [ -2  -4  -4   4  -4   4  -4   0   0   0   0   0   0   0]
            [ -2   2   0  -4   4  -4   4  -4   0   0   0   0   0   0]
            [  2  -2   2   0  -4   4  -4   4   0   0   0   0   0   0]
            [ -2   0  -2   2   0  -4   4  -4   0   0   0   0   0   0]
            [ -4  -5   0   0   0   0   0   0  -4  -5  -1   0  -1   0]
            [ -6  -9   0   0   0   0   0   0  -5 -10   0   0   0   0]
            [ -1  -1   0   0   0   0   0   0  -1   0  -2   0  -1   0]
            [  0   0   0   0   0   0   0   0   0   0   0   0   0   1]
            [ -1   0   0   0   0   0   0   0  -1   0  -1   0  -2   0]
            [  0   0   0   0   0   0   0   0   0   0   0   1   0   0]
            Signature:  (3, 11)
            Genus symbol at 2:    1^6:2^8
        """
        V = self._L
        for g in self._H0:
            V = V & kernel_sublattice(V, g - 1)
        return self._L.sublattice(V.gens())

    @cached_method
    def invariant_lattice(self):
        r"""
        Return the invariant lattice `L^H`.

        Examples::

            sage: k3 = K3SurfaceAutGrp_from_database(30,False)[0]
            sage: k3.invariant_lattice()
            Lattice of degree 22 and rank 3 over Integer Ring
            Basis matrix:
            [0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0]
            Inner product matrix:
            [ -4   4  -4   0   2  -2   2  -4   0   0   0   0   0   0   0   0   0   0   0   0   0   0]
            [  4  -4   4  -4   0   2  -2   2   0   0   0   0   0   0   0   0   0   0   0   0   0   0]
            [ -4   4  -4   4  -4   0   2  -2   0   0   0   0   0   0   0   0   0   0   0   0   0   0]
            [  0  -4   4  -4   4  -4   0   2   0   0   0   0   0   0   0   0   0   0   0   0   0   0]
            [  2   0  -4   4  -4   4  -4   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0]
            [ -2   2   0  -4   4  -4   4  -4   0   0   0   0   0   0   0   0   0   0   0   0   0   0]
            [  2  -2   2   0  -4   4  -4   4   0   0   0   0   0   0   0   0   0   0   0   0   0   0]
            [ -4   2  -2   2   0  -4   4  -4   0   0   0   0   0   0   0   0   0   0   0   0   0   0]
            [  0   0   0   0   0   0   0   0 -10  -5   0   0   0   0   0   0   0   0   0   0   0   0]
            [  0   0   0   0   0   0   0   0  -5 -10   0   0   0   0   0   0   0   0   0   0   0   0]
            [  0   0   0   0   0   0   0   0   0   0  -6   0   0   0   0   0   0   0   0   0   0   0]
            [  0   0   0   0   0   0   0   0   0   0   0   0   0   1   0   0   0   0   0   0   0   0]
            [  0   0   0   0   0   0   0   0   0   0   0   0  -2   0   0   0   0   0   0   0   0   0]
            [  0   0   0   0   0   0   0   0   0   0   0   1   0   0   0   0   0   0   0   0   0   0]
            [  0   0   0   0   0   0   0   0   0   0   0   0   0   0  -4   2  -2  -2   2   2   2  -2]
            [  0   0   0   0   0   0   0   0   0   0   0   0   0   0   2  -4   2   0   0  -2  -2   0]
            [  0   0   0   0   0   0   0   0   0   0   0   0   0   0  -2   2  -4  -2   2   2   2   0]
            [  0   0   0   0   0   0   0   0   0   0   0   0   0   0  -2   0  -2  -4   2   2   2  -2]
            [  0   0   0   0   0   0   0   0   0   0   0   0   0   0   2   0   2   2  -4   0   0   0]
            [  0   0   0   0   0   0   0   0   0   0   0   0   0   0   2  -2   2   2   0  -4  -2   2]
            [  0   0   0   0   0   0   0   0   0   0   0   0   0   0   2  -2   2   2   0  -2  -4   2]
            [  0   0   0   0   0   0   0   0   0   0   0   0   0   0  -2   0   0  -2   0   2   2  -4]

        """
        HG = self.symplectic_invariant_lattice()
        return kernel_sublattice(HG, self._g - 1)

    @cached_method
    def NS_coinvariant(self):
        """
        Return the intersection of the coinvariant lattice and the Neron-Severi
        lattice of a very general member.

        Examples::

            sage: k3 = K3SurfaceAutGrp_from_database(30,False)[0]
            sage: k3.NS_coinvariant().genus()
            Genus of
            [ -6  -5  -4   0   0   0   0   0   0   0   0]
            [ -5 -10   0   0   0   0   0   0   0   0   0]
            [ -4   0  -6   0   0   0   0   0   0   0   0]
            [  0   0   0  -4   2  -2  -2   2   2   2  -2]
            [  0   0   0   2  -4   2   0   0  -2  -2   0]
            [  0   0   0  -2   2  -4  -2   2   2   2   0]
            [  0   0   0  -2   0  -2  -4   2   2   2  -2]
            [  0   0   0   2   0   2   2  -4   0   0   0]
            [  0   0   0   2  -2   2   2   0  -4  -2   2]
            [  0   0   0   2  -2   2   2   0  -2  -4   2]
            [  0   0   0  -2   0   0  -2   0   2   2  -4]
            Signature:  (0, 11)
            Genus symbol at 2:    1^2 [2^9]_1
            Genus symbol at 5:     1^9 5^-2
        """
        L = self._L
        NS = self.NS()
        inv = self.invariant_lattice()
        return L.sublattice((
            L.orthogonal_complement(inv) & NS).gens())

    @cached_method
    def symplectic_co_invariant_lattice(self):
        r"""
        Return the coinvariant lattice `L_{H0}`. Where ``H0`` is the group of
        symplectic automorphisms.

        It is defined as the orthogonal complement of
        the symplectic invariant lattice.

        Examples::

            sage: k3 = K3SurfaceAutGrp_from_database(30,False)[0]
            sage: k3.symplectic_co_invariant_lattice().genus()
            Genus of
            [-4  2 -2 -2  2  2  2 -2]
            [ 2 -4  2  0  0 -2 -2  0]
            [-2  2 -4 -2  2  2  2  0]
            [-2  0 -2 -4  2  2  2 -2]
            [ 2  0  2  2 -4  0  0  0]
            [ 2 -2  2  2  0 -4 -2  2]
            [ 2 -2  2  2  0 -2 -4  2]
            [-2  0  0 -2  0  2  2 -4]
            Signature:  (0, 8)
            Genus symbol at 2:    2^8
        """
        return self._L.orthogonal_complement(self.symplectic_invariant_lattice())

    heart = symplectic_co_invariant_lattice

    def str_old(self):
        r"""
        A string representation of this K3 group which is used for storage.
        """
        s =  "['"+ str(self._id)+"', " + str(self._n) + ", " + str(self.connected_components()) +", ["
        s += str(self._L.basis_matrix().list()) + ","
        s += str(self._L.inner_product_matrix().list()) + "],"
        s += str([g.list() for g in self._H0]) + ","
        s += str(self._g.list()) + "]"
        return s

    def str(self):
        r"""
        A string representation of this K3 group which is used for storage.
        """
        s =  "["
        s += f"['id', '{self.id()}'],"
        s += f"['transcendental_value', {self.transcendental_value()}],"
        s += f"['connected_components', {self.connected_components()}],"
        s += f"['basis', {self.L().basis_matrix().list()}],"
        s += f"['inner_product_matrix', {self.L().inner_product_matrix().list()}],"
        s += f"['symplectic', {[g.matrix().list() for g in self.symplectic_subgroup().gens()]}],"
        s += f"['distinguished_generator', {self._g.list()}],"
        s += f"['properties',{str([(k,self._prop[k]) for k in self._prop.keys()])}]"
        s += "]"
        return s

def K3SurfaceAutGrp_from_str(s, check=False):
    r"""
    Return the K3SurfaceAutGrp defined by the given string.

    Used to read in the output of :meth:`str`.

    Examples::

        sage: k3 = K3SurfaceAutGrp_from_database(30,False)[0]
        sage: K3SurfaceAutGrp_from_str(k3.str())
        K3Group 1.30.0.1
    """
    components = None
    D0 = sage_eval(s)
    D = dict(D0)
    i = D['id']
    cc = D['connected_components']
    basis = D['basis']
    ipm = D['inner_product_matrix']
    Hs = D['symplectic']
    g = D['distinguished_generator']
    n = D['transcendental_value']
    properties = D['properties']
    k = ZZ(sqrt(ZZ(len(ipm))))
    ipm = matrix(QQ, k, k, ipm)
    l = len(basis)//k
    basis = matrix(QQ, l, k, basis)
    L = IntegralLattice(ipm, basis)
    Hs = [matrix(QQ, k, k, h) for h in Hs]
    g = matrix(QQ, k, k, g)
    return K3SurfaceAutGrp(L, Hs, g, n, check=check, name=i,
                           components=components, properties=properties)

def K3SurfaceAutGrp_from_str_old(s, check=False):
    r"""
    Return the K3SurfaceAutGrp defined by the given string.

    Used to read in the output of :meth:`str`.

    Examples::

        sage: k3 = K3SurfaceAutGrp_from_database(30,False)[0]
        sage: K3SurfaceAutGrp_from_str(k3.str())
        K3Group 1.30.0.1
    """
    components = None
    D = sage_eval(s)
    if type(D[0])==str:
        name = D[0]
        D = D[1:]
    else:
        name=None
    if len(D) == 5:
      components = D[1]
      D = [D[0]]+D[2:]
    n = D[0]
    k = ZZ(sqrt(ZZ(len(D[1][1]))))
    gram = matrix(QQ, k, k,D[1][1])
    l = len(D[1][0])//k
    basis = matrix(QQ, l, k, D[1][0])
    L = IntegralLattice(gram, basis)
    H0 = [matrix(QQ, k, k, h) for h in D[2]]
    g = matrix(QQ, k, k, D[3])
    return K3SurfaceAutGrp(L, H0, g, n, check=check, name=name,components=components)

@cached_function
def K3SurfaceAutGrp_from_database(n, pure, from1=0, to1=10000, path=None):
    r"""
    Return all effective characters of the given
    transcendental value `n`.

    Input:

      - ``n`` -- the transcendental value
      - ``pure`` -- boolean - if ``True`` returns only purely non-symplectic groups
                            - if ``False`` returns only impure groups
      - ``K3id`` -- (default:``None``) return only the K3Group of the given Id.
      - ``from1`` - start reading the database at this line
      - ``to1`` - stop reading the data base at this line

    Examples::

        sage: K3SurfaceAutGrp_from_database(7, True)
        [K3Group 0.7.0.1,
         K3Group 0.7.1.1,
         K3Group 0.7.1.2,
         K3Group 0.7.2.1,
         K3Group 0.7.2.2]


        sage: K3SurfaceAutGrp_from_database(7, False)
        [K3Group 1.7.0.1,
         K3Group 1.7.0.2,
         K3Group 1.7.1.1,
         K3Group 1.7.1.2,
         K3Group 2.7.0.1,
         K3Group 21.7.0.1,
         K3Group 21.7.0.2]
    """
    if n==2 and not pure:
        print("this may take a while")
    if euler_phi(n)>20 or n==60:
      raise ValueError("not a valid transcendental value")
    if euler_phi(n)>12 and not pure:
      return []
    if path is None:
      path = script_path
    s = os.path.join(path, "results")
    if pure:
        s = os.path.join(s, "purely_ns")
    filename = "order"
    if n < 10:
        filename+= "0"
    filename+= str(n)
    filename+= ".txt"
    s = os.path.join(s,filename)
    fi = open(s)
    output = []
    n = 0
    for f in fi:
        n+=1
        if n > to1:
          break
        if from1 <= n:
          if f[:8]!="complete":
            output.append(K3SurfaceAutGrp_from_str(f,check=False))
    fi.close()
    return output

def K3SurfaceAutGrp_from_id(K3id):
    r"""
    Return the K3Group of the given id.

    Examples::

        sage: K3SurfaceAutGrp_from_id("1.20.0.1")
        K3Group 1.20.0.1
    """
    idsplit = K3id.split('.')
    n = ZZ(idsplit[1])
    pure = (ZZ(idsplit[0])==0)

    s = os.path.join(script_path, "results")
    if pure:
        s = os.path.join(s, "purely_ns")
    filename = "order"
    if n < 10:
        filename+= "0"
    filename+= str(n)
    filename+= ".txt"
    s = os.path.join(s,filename)
    fi = open(s)
    for f in fi:
        if -1 == f[:20].find(K3id):
          continue
        fi.close()
        return K3SurfaceAutGrp_from_str(f, check=False)
    fi.close()
    raise ValueError("id not found")

def small_local_type(k3):
    assert k3.symplectic_subgroup().order()==1
    n = k3.transcendental_value()
    L = k3.L()
    g = k3._g
    return _slt(n, L, g)

def _slt(n, L, g):
    t = []
    for k in divisors(n):
      if k == 1:
        continue
      t.append((k,L.kernel_sublattice(g**(n/k)-1).genus()))
    return t

from sage.misc.cachefunc import cached_function
@cached_function
def _slt_id(k3id):
  k = K3SurfaceAutGrp_from_id(k3id)
  return small_local_type(k)
  for d in data_slt:
    if d[0] == k3id:
      s = d[1]
      return load_slt(s)

def save_slt(t):
    return [(i[0],i[1].representative().rows()) for i in t]

def load_slt(t):
    return [(i[0],Genus(matrix(i[1]))) for i in t]

def lefschetz_equations(n, euler_top, l, g):
    s = (n-1)//2
    # v_1, ... v_s, k

    K = CyclotomicField(n,"a")
    a = K.gen()
    v = zero_matrix(QQ, K.degree(), s+1)
    # setup the equation in ZZ[a]
    for i in range(s):
      t = 1/((1-a**(i+2))*(1-a**-(i+1)))
      v[:,i] = t.vector()
    v[:, s] = ((1+a)/(1-a)**2).vector()
    c = (1+a**-1)
    c-= l*(1-g)*(1+a)/(1-a)**2
    euler_top -= 2*(1-g)*l
    # set up the linear equations over ZZ
    c = c.vector()
    d = lcm(denominator(v),denominator(c))
    v = (v*d).change_ring(ZZ)
    c = vector(ZZ,(c*d).list()+[euler_top])
    a = matrix(ZZ, s*[1] + [2])
    v = v.stack(a)

    vc = v.augment(c)
    # solutions via smith form!!
    D,U,V = v.smith_form()
    cc = U*c
    y0 = zero_vector(QQ,V.ncols())
    for i in range(D.nrows()):
      if cc[i] != 0:
        y0[i] = cc[i]/D[i,i]
    x0 = V*y0
    assert v*x0 == c

    K = v.right_kernel().basis_matrix()
    ieqs=matrix(x0).stack(K).change_ring(ZZ)
    P = Polyhedron(ieqs=ieqs.columns(), backend="ppl")
    if euler_top==11 and n==42 and l==0 and g==1:
      L = [vector(p) for p in [[3,2,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1],
          [3,3,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,1]]]
    else:
      points = P.integral_points()
      L = [p*K+x0 for p in points]
    for v in L:
      assert satisfies_holomorphic_lefschetz(n=n, v=v[:-1], k=v[-1], l=l, g=g)
    return L

def satisfies_holomorphic_lefschetz(n, v, k, l, g):
    # assumes the action on the 2-form
    # is by a primitive n-th root of unity
    # order n
    # [v_0], ... v_s, s=floor(n/2) points of type i
    # diag(a^{i+1}, a^-i)
    # a curve of genus g
    # k rational curves
    K = CyclotomicField(n,"a")
    a = K.gen()
    s = (n-1)//2
    isolated_contribution = sum(v[i]/((1-a**(i+2))*(1-a**-(i+1))) for i in range(s))
    curve_contribution = (k + l*(1-g))*(1+a)/(1-a)**2
    return 1+a**-1 == curve_contribution + isolated_contribution

def satisfies_lefschetz(v, k, l, g):
    euler = 2 + g.trace()
    return euler == sum(v) + k*2 + l*(2-2*g)

def possible_type(k3):
    n = k3.n()
    L = k3.L()
    g = k3._g
    return possible_types(n,L,g)

def possible_fixed_points_curve(p, genus):
    # 2*genus - 2 = p*(2g_X-2) +(p-1) fixed_points
    # how about non-prime order? might be useful
    res = []
    for n in range((2*genus/(p-1)+2).floor()+1):
      gy = (ZZ(2*genus-2 - (p-1)*n)/p+2)/2
      if gy.denominator()==1:
        res.append(n)
    return res

def fixed_points(k3):
    n = k3.n()
    L = k3.L()
    g = k3._g
    return _fixed_points(n, L, g)

def _fixed_points(n, L, g):
    if is_prime(n):
      return _type_prime(n, L, g)
    T = possible_types(n, L, g)
    if  len(T) == 1:
      return T[0]
    t = _slt(n, L , g)  # any problems if slt is not good enough??
    # some fixed points require geometric arguments
    # these are hardcoded
    if n == 4:
      ambig = [
        # g^2-fixed elliptic fibration
        ['0.4.1.3', ((4,), 0, 0, 1)],
        ['0.4.1.4', ((4,), 0, 0, 1)],
        ['0.4.1.5', ((4,), 0, 0, 1)],
        ['0.4.1.6', ((4,), 0, 1, 1)],

        ['0.4.2.3', ((4,), 0, 1, 1)],
        ['0.4.2.4', ((4,), 0, 0, 1)],
        ['0.4.2.5', ((4,), 0, 0, 1)],

        ['0.4.3.6', ((4,), 0, 0, 1)],
        ['0.4.3.7', ((4,), 0, 0, 1)],
        ['0.4.3.8', ((4,), 0, 1, 1)],

        ['0.4.4.7', ((4,), 0, 1, 1)],
        ['0.4.4.8', ((4,), 0, 0, 1)],
        ['0.4.4.9', ((4,), 0, 0, 1)],

        ['0.4.3.9',  ((6,), 1, 0, 1)],
        ['0.4.3.10', ((6,), 1, 0, 1)],
        ['0.4.3.11', ((6,), 1, 1 ,1)],
        # genus 2
        ['0.4.1.1', ((4,), 0, 0, 1)],
        ['0.4.2.1', ((4,), 0, 0, 1)],
        ['0.4.3.3', ((4,), 0, 0, 1)],
        ['0.4.4.4', ((4,), 0, 0, 1)],
        # two g^2-fixed elliptic fibrations
        ['0.4.5.12',((4,), 0, 0, 1)],
        ['0.4.5.14',((4,), 0, 1, 1)]]
    if n == 6:
      ambig = [
        ["0.6.2.29", ((2,2) ,0, 1, 1)],
        ["0.6.3.36", ((3, 0), 0, 1, 1)]]
    if n == 8:
      t = _slt(n, L, g)
      ambig=[
        ['0.8.1.7', ((1, 1, 2), 0, 0, 1)],
        ['0.8.1.8', ((1, 1, 2), 0, 1, 1)],
        ['0.8.2.8', ((2, 0, 0), 0, 1, 1)],
        ['0.8.2.10', ((2, 0, 0), 0, 0, 1)]]
    if n == 9:
      ambig = [
        ["0.9.1.3", ((1,1,1,0),0,0,1)],
        ["0.9.1.4", ((1,1,1,0),0,1,1)]]
    if n == 10:
      ambig = [
      ["0.10.2.1",((0,0,1,6), 0, 0, 1)],
      ["0.10.1.11",((1,1,0,0),0, 1, 1)]]
    if n == 12:
      ambig = [['0.12.1.21', ((1, 0, 0, 0, 0), 0, 1, 1)]]
    for g in ambig:
      if _slt_id(g[0]) == t:
        return g[1]
    assert False

def _type_prime(p, L, f):
    euler = 2 + f.trace()
    s = (p-1)//2
    pfix = dict()
    F = L.kernel_sublattice(f-1)
    r = F.rank()
    a = F.det().valuation(p)
    g = (22 - r - (p-1)*a)/(2*(p-1))
    o = F.genus().local_symbol(2).symbol(1)[3]
    l = 1
    if p == 2:
      nisolated = 0
      k = (r-a)/2
      if r == 10 and a==10 and o == 0:
        return ((), 0, 0, 1)
      elif r == 10 and a==8 and o == 0:
        return ((), 0, 2, 1)
    if p in [3,5,7]:
      nisolated = (-2 +(p-2)*r)/(p-1)
      k = (2 + r - (p-1)*a)/(2*(p-1))
    if p == 11:
      nisolated = (2 +(p-2)*r)/(p-1)
      k = (-2 + r - (p-1)*a)/(2*(p-1))
    if p == 13:
      nisolated = 9
      k = 0 # there seems to be an error in [AST] Thm 0.1
    if p == 17:
      nisolated = 7
      k = -1
    if p == 19:
      nisolated = 5
      k = -1
    if k == -1:
      k = 0
      l = 0
      g = 1
    if g == 0:
      k += 1
      l = 0
      g = 1
    leq = lefschetz_equations(p, euler,l=l, g=g)
    assert len(leq)==1
    assert leq[0][-1] == k
    assert sum(leq[0][:-1])==nisolated
    return leq[0][:-1], k, l, g


def _possible_types_n(n, L, g):
    # assume for now g purely non-symplectic of order n
    assert not is_prime(n)

    euler = 2 + g.trace()
    # bound on number of rational curves fixed
    ratsbound = 100
    genera = set()
    l = 2
    for k in divisors(n):
      if k in [1, n]:
        continue
      vk, nratk, lk, genusk = _fixed_points(k, L, g**(n/k))
      ratsbound = min(nratk, ratsbound)
      l = min(l,lk)
      genera.add(genusk)

    if len(genera) > 1:
      genus = 1
      l = 0
    else:
      genus = genera.pop()

    data = []
    for nrational in range(ratsbound+1):
      nisolated = euler - 2*nrational
      n1 = nisolated
      # no fixed curve
      if n1 >=0:
        data.append([nisolated, nrational, 0, 1])
      # l fixed curves of genus g
      n2 = nisolated - l*(2-2*genus)
      if n2 >= 0:
        for ll in range(1,l+1):
          data.append([n2, nrational, ll, genus])
        # no fixed curve

    res = []
    for nisolated, nrational, l, genus in data:
      for v in lefschetz_equations(n=n, euler_top=euler, l=l,g=genus):
        if nrational == v[-1]: #doppelt gemoppelt
          res.append((v[:-1],v[-1],l, genus))
    return res

def possible_types(n, L, g):
    if is_prime(n):
      return _type_prime(n, L, g)
    types = _possible_types_n(n, L, g)
    R = PolynomialRing(QQ,"x"); x = R.gen()
    rF = L.kernel_sublattice(g-1).rank()
    # compatibility checks
    res = []
    for v,rat,l,genus in types:
      iscompatible = True
      for p in prime_divisors(n):
        vp,ratp, lp, genusp = _fixed_points(n//p, L, g**p)
        ftp = fixed_types_to_nn(n, v, n//p)  # fixed points of f_p = f^p
        if not all(ftp[i]<=vp[i] for i in range(len(ftp))):
          iscompatible = False
          break
        if not all(mod(ftp[i]-vp[i],p)==0 for i in range(len(ftp))):
          iscompatible = False
          break

        # number of isolated fixed points of f
        # lying on curves fixed by f_p = g^p
        fixed_points_on_curve = sum(v) -  sum(ftp)
        # f acts on the fixed locus of f^p
        # with order dividing p
        a = g.charpoly(x).valuation(x-1)
        b = g.charpoly(x).valuation(cyclotomic_polynomial(p,x))

        # upper bound on the number of cycles
        c_ub = min( a - rat - (0 if genus is None else 1), b)
        #(ratp - rat) = f + c <= a + cc
        finv_lb = ratp - rat - c_ub*p
        finv_lb = max(finv_lb, (ratp-rat) % p) #could be made more precise
        # we know finv \equiv ratp-rat mod p
        # 2 fixed points per not-fixed invariant rational curve
        finv_ub = ratp - rat
        nfix_ub = 2*(finv_ub)
        nfix_lb = 2*(finv_lb)
        # additional fixed points on the genus g >=1 curve (if it exists)
        if lp > l:
          # a genus g curve is not fixed entirely by f
          assert genusp > 0
          FixC = possible_fixed_points_curve(p, genusp)
          if len(FixC) == 0:
            iscompatible = False
          nfix_ub += (lp-l)*max(FixC)
          nfix_lb += (lp-l)*min(FixC)

        # check the bounds
        if fixed_points_on_curve > nfix_ub:
          iscompatible = False
          break
        if nfix_lb > fixed_points_on_curve:
          iscompatible = False
          break
      if iscompatible:
        res.append((v,rat,l,genus))
    return res

# the new fixed points are upper bounds
def fixed_types_to_nn(n, v, nn):
  assert ZZ(nn).divides(n)
  ss = (nn-1)//2
  vv = ss*[0]
  for i in range(1,len(v)+1):
    ai = v[i-1]
    ii1 = (i+1) % nn
    ii2 = -i % nn
    if ii1 > ii2:
      ii1,ii2 = ii2,ii1
    if ii1 != 0:
      vv[ii1-2] += ai
  return vv

def connected_comps(L):
  B = list(L.gens())
  C = []
  while len(B)>0:
    b = B.pop()
    C.append(b)
    for c in B:
      if b.inner_product(c)!=0:
        C.append(c)
        B.remove(c)
  assert len(C)==L.rank()
  return L.sublattice(C)

def connected_comps(L):
  L = L.lll()
  B = list(L.gens())
  C = []
  ADE = []
  while len(B)>0:
    CC = []
    b = B.pop()
    CC.append(b)
    while True:
      for c in B:
        if any([a.inner_product(c)!=0 for a in CC]):
          CC.append(c)
          B.remove(c)
          break
      else:
        break
    S = L.sublattice(CC)
    ADE.append(ADE_type(S.gram_matrix()))
    C += CC
  assert len(C)==L.rank()
  return L.sublattice(C),ADE

def ADE_type(G):
  r = rank(G)
  d = G.det().abs()
  if r == 8 and d==1:
    return "E8"
  if r == 7 and d == 2:
    return "E7"
  if r == 6 and d ==3:
    return "E6"
  if d == r + 1:
    return "A%s"%r
  if d == 4:
    return "D%s"%r
  raise ValueError("Not an ADE lattice")

def fiber_types(L, f):
  assert f in L
  assert f.inner_product(f)==0
  K = L.orthogonal_submodule([f])/L.span([f])
  K = L.sublattice([g.lift() for g in K.gens()])
  R = K.root_sublattice()
  return connected_comps(R)[1], (K/R).invariants()

def random_isotropic(F):
  assert prod(F.signature_pair())!=0
  while True:
    e = F.random_element()
    if e!=0 and e.inner_product(e) == 0:
      return F(e/gcd(F.coordinates(e)))

def root_type_f(N, f):
  f = N(f)
  R = N.orthogonal_submodule([f])
  R = R/R.submodule([f])
  R = N.sublattice([g.lift() for g in R.gens()])
  R = R.root_sublattice().lll()
  return connected_comps(R)[1],R.rank()


def has_invariant_section(k3,k):
  F = k3.invariant_lattice()
  NS = k3.NS()
  N = NS.kernel_sublattice(k3._g**k-1)
  if not F.discriminant_group().is_genus((0,F.rank()-2)):
    return False
  while True:
    f = random_isotropic(F)
    rt,rank = root_type_f(N,f)
    print(rt,N.rank())
    if rank < N.rank()-2 or (rank ==0):
      break
  s = gcd([f.inner_product(b) for b in F.basis()])
  return s==1

def fixed_points_power(k3,k):
  assert k.divides(k3.transcendental_value())
  return _fixed_points(k3.transcendental_value()//k, k3.L(), k3._g**k)

def restrict(F,f):
  return matrix([F.coordinates(b*f) for b in F.basis()])
