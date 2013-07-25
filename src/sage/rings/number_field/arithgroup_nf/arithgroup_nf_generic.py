r"""
Arithmetic (finite index) subgroups of `{\rm SL}_2(\mathcal{O})`
and `{\rm GL}_2(\mathcal{O})`) where `\mathcal{O}` is an (maximal) order in a number field `K`.

AUTHORS:

- Fredrik Stromberg (2013): initial version based on arithgroup_generic.py

"""

################################################################################
#
#       Copyright (C) 2009, The Sage Group -- http://www.sagemath.org/
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#
#
################################################################################


from sage.rings.all import ZZ,Integer
from sage.misc.cachefunc import cached_method
from copy import copy # for making copies of lists of cusps
from sage.modular.cusps_nf import NFCusp
from sage.misc.latex import latex
from sage.groups.matrix_gps.linear import LinearMatrixGroup_generic


class ArithmeticSubgroup_NF(LinearMatrixGroup_generic):
    r"""
    Base class for arithmetic subgroups of `{\rm SL}_2(K)`. Not
    intended to be used directly, but still includes quite a few
    general-purpose routines which compute data about an arithmetic subgroup.

    """

    def __init__(self,ring,special=True,name='',ltx=''):
        r"""
        Standard init routine.

        EXAMPLE:


        """
        degree = 2
        if name=='':
            name = 'Arithmetic Subgroup of the Special Linear Group of degree {0} over {1}'.format(degree, ring)
        if ltx=='':
            ltx  = 'GL({0}, {1})'.format(degree, latex(ring))
  
        super(LinearMatrixGroup_generic,self).__init__(Integer(degree),ring,special,name,ltx)
        self._base_ring = ring
        if ring == ZZ:
            self._base_field = QQ 
        elif hasattr(ring,"number_field"):
            self._base_field = ring.number_field()
        else:
            raise NotImplementedError


    def __reduce__(self):
        r"""
        Used for pickling self.

        NOTE: This function should be overridden by all subclasses.

        EXAMPLES::

            sage: sage.modular.arithgroup.arithgroup_generic.ArithmeticSubgroup().__reduce__()
            Traceback (most recent call last):
            ...
            NotImplementedError: all subclasses must define a __reduce__ method
        """
        raise NotImplementedError, "all subclasses must define a __reduce__ method"

    def __hash__(self):
        r"""
        Return a hash of self.
        
        EXAMPLES::

          
        """
        return hash(str(self))

    def coset_reps(self, G=None):
        r"""
        Return right coset representatives for self \\ G, where G is another
        arithmetic subgroup that contains self.  If G = None, default to G =
        SL2Z.
        
        For generic arithmetic subgroups G this is carried out by Todd-Coxeter
        enumeration; here G is treated as a black box, implementing nothing but
        membership testing.

        EXAMPLES::

        """
        raise NotImplementedError

    def nu(self,order=2):
        r"""
        Return the number of orbits of elliptic points of given order for this
        arithmetic subgroup.        

        """
        raise NotImplementedError
            
    def nu2(self):
        r"""
        Return the number of orbits of elliptic points of order 3 for this
        arithmetic subgroup.
        """
        return self.nu(2)    
    def nu3(self):
        r"""
        Return the number of orbits of elliptic points of order 3 for this
        arithmetic subgroup.
        """
        return self.nu(3)    

    def orders_of_elliptic_elements(self):
        r"""
        Returns the possible orders of elliptic elements.
        """
        raise NotImplementedError
    
    def __cmp__(self, other):
        r"""
        Compare self to other.

        NOTE: This function must be overridden by all subclasses.
        
        EXAMPLES::

            sage: sage.modular.arithgroup.arithgroup_generic.ArithmeticSubgroup().__cmp__(ZZ)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def is_abelian(self):
        r"""
        Return True if this arithmetic subgroup is abelian.

        Since arithmetic subgroups are always nonabelian, this always
        returns False.

        EXAMPLES::


        """
        return False

    def is_finite(self):
        r"""
        Return True if this arithmetic subgroup is finite.

        Since arithmetic subgroups are always infinite, this always
        returns False.

        EXAMPLES::

        """
        return False

    def is_subgroup(self, right):
        r"""
        Return True if self is a subgroup of right, and False otherwise. For
        generic arithmetic subgroups this is done by the absurdly slow
        algorithm of checking all of the generators of self to see if they are
        in right.

        EXAMPLES::

            sage: sage.modular.arithgroup.arithgroup_generic.ArithmeticSubgroup().is_subgroup(SL2Z)
            Traceback (most recent call last):
            ...
            NotImplementedError
            sage: sage.modular.arithgroup.arithgroup_generic.ArithmeticSubgroup.is_subgroup(Gamma1(18), Gamma0(6))
            True
        """
        # ridiculously slow generic algorithm
        w = self.gens()
        for g in w:
            if not (g in right):
                return False
        return True

    def is_normal(self):
        r"""
        Return True precisely if this subgroup is a normal subgroup of SL_2(R).

        EXAMPLES::
            

        """
        G = ArithmeticSubgroup_NF(self.base_ring(),self.is_special())
        for x in self.gens():
            for y in G.gens():
                if y*G(x)*(~y) not in self:
                    return False
        return True

    def is_odd(self):
        r"""
        Return True precisely if this subgroup does not contain the
        matrix -1.

        EXAMPLES::


        """
        return not self.is_even()

    def is_even(self):
        r"""
        Return True precisely if this subgroup contains the matrix -1.

        EXAMPLES::

            sage: SL2Z.is_even()
            True
            sage: Gamma0(20).is_even()
            True
            sage: Gamma1(5).is_even()
            False
            sage: GammaH(11, [3]).is_even()
            False
        """
        G = ArithmeticSubgroup_NF(self.base_ring(),self.is_special())
        minus_one = G([-1,0,0,-1])
        return not minus_one in self

    def to_even_subgroup(self):
        r"""
        Return the smallest even subgroup of `SL(2, \ZZ)` containing self.

        EXAMPLE::

            sage: sage.modular.arithgroup.arithgroup_generic.ArithmeticSubgroup().to_even_subgroup()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        if self.is_even():
            return self
        else:
            raise NotImplementedError

    def order(self):
        r"""
        Return the number of elements in this arithmetic subgroup.

        Since arithmetic subgroups are always infinite, this always returns
        infinity.

        EXAMPLES::

            sage: SL2Z.order()
            +Infinity
            sage: Gamma0(5).order()
            +Infinity
            sage: Gamma1(2).order()
            +Infinity
            sage: GammaH(12, [5]).order()
            +Infinity
        """
        from sage.rings.infinity import infinity
        return infinity

    def reduce_cusp(self, c):
        r"""
        Given a cusp `c \in \mathbb{P}^1(\QQ)`, return the unique reduced cusp
        equivalent to c under the action of self, where a reduced cusp is an
        element `\tfrac{r}{s}` with r,s coprime non-negative integers, s as
        small as possible, and r as small as possible for that s.
        
        NOTE: This function should be overridden by all subclasses.

        EXAMPLES::

            sage: sage.modular.arithgroup.arithgroup_generic.ArithmeticSubgroup().reduce_cusp(1/4)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def cusps(self, algorithm='default'):
        r"""
        Return a sorted list of inequivalent cusps for self, i.e. a set of
        representatives for the orbits of self on `\mathbb{P}^1(\QQ)`.
        These should be returned in a reduced form where this makes sense.

        INPUTS:
        
        - ``algorithm`` -- which algorithm to use to compute the cusps of self.
          ``'default'`` finds representatives for a known complete set of
          cusps. ``'modsym'`` computes the boundary map on the space of weight
          two modular symbols associated to self, which finds the cusps for
          self in the process.

        EXAMPLES::


        """
        try:
            return copy(self._cusp_list[algorithm])
        except (AttributeError,KeyError):
            self._cusp_list = {}
        if self.ncusps()==1:
            s = [NFCusp(1,0)]

        if algorithm == 'default':
            s = self._find_cusps()
        else:
            raise ValueError, "unknown algorithm: %s"%algorithm

        self._cusp_list[algorithm] = s
        return copy(s)

    def _find_cusps(self):
        r"""
        Calculate a list of inequivalent cusps. 

        EXAMPLES::



        NOTE: This should be implemented in subclasses.
              (or if a generic algorithm 

        """
        raise NotImplementedError

            
    def are_equivalent(self, x, y, trans = False):
        r""" 
        Test whether or not cusps x and y are equivalent modulo self.
        
        If self
        has a reduce_cusp() method, use that; otherwise do a slow explicit
        test. 

        If trans = False, returns True or False. If trans = True, then return
        either False or an element of self mapping x onto y.

        EXAMPLE::


        """
        if hasattr(self,'reduce_cusp'):
            if self.reduce_cusp(x)==self.reduce_cusp(y):
                return True
            return False
        else:
            return NotImplementedError

        
    def cusp_data(self, c):
        r"""
        Return a triple (g, w, t) where g is an element of self generating the
        stabiliser of the given cusp, w is the width of the cusp, and t is 1 if
        the cusp is regular and -1 if not.

        EXAMPLES::

            sage: Gamma1(4).cusp_data(Cusps(1/2))
            (
            [ 1 -1]
            [ 4 -3], 1, -1
            )
        """
        raise NotImplementedError

            
    def index(self):
        r"""
        Return the index of self in the full modular group.

        EXAMPLES::

            sage: Gamma0(17).index()
            18
            sage: sage.modular.arithgroup.congroup_generic.CongruenceSubgroup(5).index()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        return len(list(self.coset_reps()))

    # Question: Does any of these algorithms apply for number fields.
    #
    # def is_regular_cusp(self, c):
    #     r"""
    #     Return True if the orbit of the given cusp is a regular cusp for self,
    #     otherwise False. This is automatically true if -1 is in self.

    #     EXAMPLES::

    #         sage: Gamma1(4).is_regular_cusp(Cusps(1/2))
    #         False
    #         sage: Gamma1(4).is_regular_cusp(Cusps(oo))
    #         True
    #     """
    #     if self.is_even(): return True
    #     return (self.cusp_data(c)[2] == 1)    
    # def cusp_width(self, c):
    #     r"""
    #     Return the width of the orbit of cusps represented by c.

    #     EXAMPLES::

    #         sage: Gamma0(11).cusp_width(Cusps(oo))
    #         1
    #         sage: Gamma0(11).cusp_width(0)
    #         11
    #         sage: [Gamma0(100).cusp_width(c) for c in Gamma0(100).cusps()] 
    #         [100, 1, 4, 1, 1, 1, 4, 25, 1, 1, 4, 1, 25, 4, 1, 4, 1, 1]
    #     """
    #     return self.cusp_data(c)[1]

    
    # def generalised_level(self):
    #     r"""
    #     Return the generalised level of self, i.e. the least common multiple of
    #     the widths of all cusps. 
        
    #     If self is *even*, Wohlfart's theorem tells us that this is equal to
    #     the (conventional) level of self when self is a congruence subgroup.
    #     This can fail if self is odd, but the actual level is at most twice the
    #     generalised level. See the paper by Kiming, Schuett and Verrill for
    #     more examples.
       
    #     EXAMPLE::

    #         sage: Gamma0(18).generalised_level()
    #         18
    #         sage: sage.modular.arithgroup.congroup_generic.CongruenceSubgroup(5).index()
    #         Traceback (most recent call last):
    #         ...
    #         NotImplementedError

    #     In the following example, the actual level is twice the generalised
    #     level. This is the group `G_2` from Example 17 of K-S-V.
        
    #     ::

    #         sage: G = CongruenceSubgroup(8, [ [1,1,0,1], [3,-1,4,-1] ])
    #         sage: G.level()
    #         8
    #         sage: G.generalised_level()
    #         4
    #     """
    #     return arith.lcm([self.cusp_width(c) for c in self.cusps()])

    # def projective_index(self):
    #     r"""
    #     Return the index of the image of self in `{\rm PSL}_2(\ZZ)`. This is equal
    #     to the index of self if self contains -1, and half of this otherwise.

    #     This is equal to the degree of the natural map from the modular curve
    #     of self to the `j`-line.
        
    #     EXAMPLE::

    #         sage: Gamma0(5).projective_index()
    #         6
    #         sage: Gamma1(5).projective_index()
    #         12
    #     """
        
    #     if self.is_even():
    #         return self.index()
    #     else:
    #         return self.index() // 2

    def is_congruence(self):
        r"""
        Return True if self is a congruence subgroup.

        EXAMPLE::

            sage: Gamma0(5).is_congruence()
            True
            sage: sage.modular.arithgroup.arithgroup_generic.ArithmeticSubgroup().is_congruence()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """

        raise NotImplementedError


    @cached_method
    def generators(self):
        r"""
        Return a list of generators for this congruence subgroup. The result is cached.

        
        EXAMPLE::


        """
        raise NotImplementedError

    def gens(self, *args, **kwds):
        r"""
        Return a tuple of generators for this congruence subgroup.

        The generators need not be minimal. For arguments, see :meth:`~generators`.

        EXAMPLES::

        """
        return tuple(self.generators(*args, **kwds))

    def gen(self, i):
        r"""
        Return the i-th generator of self, i.e. the i-th element of the
        tuple self.gens().

        EXAMPLES::

        """
        return self.generators()[i]

    def ngens(self):
        r"""
        Return the size of the minimal generating set of self returned by
        :meth:`generators`.

        EXAMPLES::


        """
        return len(self.generators())

    def ncusps(self):
        r"""
        Return the number of cusps of this arithmetic subgroup. This is
        provided as a separate function since for dimension formulae in even
        weight all we need to know is the number of cusps, and this can be
        calculated very quickly, while enumerating all cusps is much slower.

        EXAMPLES::

            sage: sage.modular.arithgroup.arithgroup_generic.ArithmeticSubgroup.ncusps(Gamma0(7))
            2
        """

        return ZZ(len(self.cusps()))


    def dimension_modular_forms(self, k=2):
        r"""
        Return the dimension of the space of weight k modular forms for this
        group.

        EXAMPLE::

        """
        raise NotImplementedError

    def dimension_cusp_forms(self, k=2):
        r"""
        Return the dimension of the space of weight k cusp forms for this
        group.

        EXAMPLE::

        """
        raise NotImplementedError
        
    def dimension_eis(self, weight=(2,)):
        r"""
        Return the dimension of the space of  Eisenstein series of the given weight for
        this group, which is a subspace of the space of modular forms
        complementary to the space of cusp forms.

        INPUT:
        
        - ``k`` - an integer (default 2).

        EXAMPLES::


        """

        raise NotImplementedError

    
    def sturm_bound(self, weight=(2,)):
        r"""
        Returns the Sturm bound for modular forms of the given weight and level
        this subgroup.
        
        INPUT:
        
        -  ``weight`` - an tuple of integers `\geq 2` (default: 2)
        
        EXAMPLES::

        """
        
        raise NotImplementedError

