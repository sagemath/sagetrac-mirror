r"""
Testing origamis
"""

import random
from sage.misc.misc import cputime

from origami_dense import Origami_dense_pyx
from origami import Origami

from sage.combinat.permutation import Permutations

from copy import copy

def random_origami(n,connected=True):
    r"""
    Returns a random origami on n squares.

    INPUT:

    - ``n`` - integer - the number of squares
    
    - ``connected`` - boolean (default: True) - whether or not the output
      origami should be connected.
    """
    G = Permutations(n)
    r = [i-1 for i in G.random_element()]
    u = [i-1 for i in G.random_element()]
    o = Origami(r,u,check=False,as_tuple=True)
    if connected:
       while not o.is_connected():
           r = [i-1 for i in G.random_element()]
           u = [i-1 for i in G.random_element()]
           o = Origami(r,u,check=False,as_tuple=True)
    return o

class Test:
    r"""
    Origami testing class

    Beware: some of the tests fail for not connected graphs

    INPUT:

    - ``nb_squares`` - the number of squares for the origami on which the tests
      are performed

    - ``connected`` - whether or not consider only connected origamis
    """
    def __init__(self,min_nb_squares=6,max_nb_squares=20):
        self.min_nb_squares = min_nb_squares
        self.max_nb_squares = max_nb_squares

    def __repr__(self):
        return "Origami test for squares between %d and %d"%(self.min_nb_squares,self.max_nb_squares)

    def _do(self, name):
        """
        Perform the test 'test_name', where name is specified as an
        argument. This function exists to avoid a call to eval.

        EXAMPLES:
            sage: sage.modular.modsym.tests.Test()._do("random") 
            test_random
            ...
        """
        Test.__dict__["test_%s"%name](self)

    def random(self, seconds=0):
        """
        Perform random tests for a given number of seconds, or
        indefinitely if seconds is not specified.

        EXAMPLES:
            sage: sage.dynamics.flat_surfaces.origamis.tests.Test().random(1)
            test_random
            ...
        """
        self.test("random", seconds)
        
    def test(self, name, seconds=0):
        """
        Repeatedly run 'test_name', where name is passed as an
        argument. If seconds is nonzero, run for that many seconds. If
        seconds is 0, run indefinitely.
        """
        seconds = float(seconds)
        total = cputime()
        n = 1
        while seconds == 0 or cputime(total) < seconds:
            s = "** test_origami: number %s"%n
            if seconds > 0:
                s += " (will stop after about %s seconds)"%seconds
            t = cputime()
            self._do(name)
            print "\ttime=%s\telapsed=%s"%(cputime(t),cputime(total))
            n += 1

    def _get_random_origami(self):
        r"""
        Compute a random origami with the good parameters
        """
        return random_origami(random.randint(self.min_nb_squares,self.max_nb_squares))

    #
    # individual tests with randomize input
    #

    def test_random(self):
        """
        Do a random test from all the possible tests.
        """
        tests = [a for a in Test.__dict__.keys() if a[:5] == "test_" and a != "test_random"]
        name = random.choice(tests)
        print "Doing random test"
        Test.__dict__[name](self)

    def test_standard_form(self):
        r"""
        Test that the standard form is standard.
        """
        from sage.all import SymmetricGroup
        o = self._get_random_origami()
        o._set_standard_form()
        G = SymmetricGroup(o.nb_squares())

        for _ in xrange(10):
            p = G.random_element()

            oo = copy(o)
            oo.relabel(p)
            oo_std = oo.to_standard_form()

            if o_std != oo_std:
                raise AssertionError, "no equality between relabeled\n%s\n- - - - - - - -\n%s" %(o,oo)

    def test_comparison(self):
        r"""
        Test whether the comparison of dense origami is well implemented.
        """
        o1 = self._get_random_origami()
        o2 = self._get_random_origami()

        i1 = o1.r_tuple()
        j1 = o1.u_tuple()
        i2 = o2.r_tuple()
        j2 = o2.u_tuple()

        print "test comparison between (%s,%s) and (%s,%s)" %(i1,j1,i2,j2)

        oo1 = Origami_dense_pyx(i1,j1)
        oo2 = Origami_dense_pyx(i2,j2)
        if (oo1 == oo2) != ((i1,j1) == (i2,j2)):
            raise AssertionError,  "equality between (%s,%s) and (%s,%s)" %(i1,j1,i2,j2)
        if (oo1 != oo2) != ((i1,j1) != (i2,j2)):
             raise AssertionError, "difference between (%s,%s) and (%s,%s)" %(i1,j1,i2,j2)
        if (oo1 < oo2) != ((i1,j1) < (i2,j2)):
             raise AssertionError, "less than between (%s,%s) and (%s,%s)" %(i1,j1,i2,j2)
        if (oo1 <= oo2) != ((i1,j1) <= (i2,j2)):
            raise AssertionError, "leq than between (%s,%s) and (%s,%s)" %(i1,j1,i2,j2)
        if (oo1 > oo2) != ((i1,j1) > (i2,j2)):
            raise AssertionError, "greater than between (%s,%s) and (%s,%s)" %(i1,j1,i2,j2)
        if (oo1 >= oo2) != ((i1,j1) >= (i2,j2)):
            raise AssertionError, "geq than between (%s,%s) and (%s,%s)" %(i1,j1,i2,j2)

    def test_orbit_and_action(self):
        r"""
        Test whether the functions horizontal_twist and vertical_twist are
        compatible with the graph obtained by the method gl2z_edges
        """
        o = self._get_random_origami()

        l, i = o.gl2z_edges()

        for oo in l:
            if oo != oo.horizontal_twist(1).horizontal_twist(-1):
                raise AssertionError, "horiz. twist back and forth not identity (%s,%s)" %(oo.r_tuple(),oo.u_tuple())
            if oo != oo.horizontal_twist(-1).horizontal_twist(1):
                raise AssertionError, "horiz. twist forth and back forth not identity (%s,%s)" %(oo.r_tuple(),oo.u_tuple())
            if oo.horizontal_twist(2) != oo.horizontal_twist(1).horizontal_twist(1):
                raise AssertionError, "multiple positive horiz. twist (%s,%s)" %(oo.r_tuple,oo.u_tuple)
            if oo.horizontal_twist(-2) != oo.horizontal_twist(-1).horizontal_twist(-1):
                raise AssertionError, "multiple negative horiz. twist (%s,%s)" %(oo.r_tuple(),oo.u_tuple())

            ooo = oo.horizontal_twist()
            ooo._set_standard_form()
            if ooo != l[oo]:
                raise AssertionError, "horiz. twist with (%s,%s)" %(oo.r_tuple(),oo.u_tuple())

            if oo != oo.vertical_twist(1).vertical_twist(-1):
                raise AssertionError, "vert. twist back and forth not identity (%s,%s)" %(oo.r_tuple(),oo.u_tuple())
            if oo != oo.vertical_twist(-1).vertical_twist(1):
                raise AssertionError, "vert. twist forth and back forth not identity (%s,%s)" %(oo.r_tuple(),oo.u_tuple())
            if oo.vertical_twist(2) != oo.vertical_twist(1).vertical_twist(1):
                raise AssertionError, "multiple positive vert. twist (%s,%s)" %(oo.r_tuple,oo.u_tuple)
            if oo.vertical_twist(-2) != oo.vertical_twist(-1).vertical_twist(-1):
                raise AssertionError, "multiple negative vert. twist (%s,%s)" %(oo.r_tuple(),oo.u_tuple())

            ooo = oo.vertical_twist()
            ooo._set_standard_form()
            if ooo != i[l[i[oo]]]:
                raise AssertionError, "vert. twist with (%s,%s)" %(oo.r_tuple(),oo.u_tuple())

    def test_cylinder_coordinates(self):
        r"""
        Test cylinder coordinates and the reconstruction of an origami from a
        cylinder diagrams and the 3-tuple ``(lengths,heights,twists)``.
        """
        o = self._get_random_origami()
        o._set_standard_form()

        c,lengths,heights,twists = o.cylinder_diagram(True)
        oo = c.cylcoord_to_origami(lengths,heights,twists)
        oo._set_standard_form()

        if o != oo:
            raise AssertionError, "cylinder coordinates on (%s,%s)" %(o.r_tuple(),o.u_tuple())

    def test_stratum_component(self):
        r"""
        Test that the method stratum component is constant along the GL2(Z)
        orbit of an origami.
        """
        if self.max_nb_squares <= 5:
            return
        o = self._get_random_origami()
        s = o.stratum()
        while s.is_connected():
            o = self._get_random_origami()
            s = o.stratum()

        cc = o.stratum_component()
        oo = o
        for _ in xrange(20):
            if random.randint(0,2): # proba 2/3: we apply a twist
                oo = Origami(oo.r(), oo.u()*oo.r())
            else: # proba 1/3: we exchange
                oo = Origami(oo.u(), oo.r())
            if oo.stratum_component() != cc:
                raise AssertionError, "the two following origamis are in the same GL2(Z) orbit but do not have the same stratum component r1=%s, u1=%s\nr2=%s, u2=%s"%(str(o.r()),str(o.u()),str(oo.r()),str(oo.u()))

    #
    # exhaustive tests
    #

    def test_comparisons(self):
        from itertools import product,tee,permutations

        for i1,j1,i2,j2 in product(*tee(permutations(range(n)),4)):
            o1 = Origami_dense_pyx(i1,j1)
            o2 = Origami_dense_pyx(i2,j2)
            if (o1 == o2) != ((i1,j1) == (i2,j2)):
                raise AssertionError,  "equality between (%s,%s) and (%s,%s)" %(i1,j1,i2,j2)
            if (o1 != o2) != ((i1,j1) != (i2,j2)):
                 raise AssertionError, "difference between (%s,%s) and (%s,%s)" %(i1,j1,i2,j2)
            if (o1 < o2) != ((i1,j1) < (i2,j2)):
                 raise AssertionError, "less than between (%s,%s) and (%s,%s)" %(i1,j1,i2,j2)
            if (o1 <= o2) != ((i1,j1) <= (i2,j2)):
                raise AssertionError, "leq than between (%s,%s) and (%s,%s)" %(i1,j1,i2,j2)
            if (o1 > o2) != ((i1,j1) > (i2,j2)):
                raise AssertionError, "greater than between (%s,%s) and (%s,%s)" %(i1,j1,i2,j2)
            if (o1 >= o2) != ((i1,j1) >= (i2,j2)):
                raise AssertionError, "geq than between (%s,%s) and (%s,%s)" %(i1,j1,i2,j2)


class Test_Origami:
    r"""
    Test various functionnality of origamis
    """
    
    def test_projective_standard_form(self, o):
        from sage.all import SymmetricGroup
        G = SymmetricGroup(o.nb_squares())

        o_std = o.projective_standard_form()
        o_inv_std = o.inverse().projective_standard_form()

        if o_std != o_inv_std:
            raise AssertionError, "standard forms not projectively equal for inverses\n%s" %(o)

        for _ in xrange(10):
            p = G.random_element()

            oo = copy(o)
            oo.relabel(p)
            oo_std = o.projective_standard_form()
            oo_inv_std = oo.inverse().projective_standard_form()

            if oo_std != oo_inv_std:
                raise AssertionError, "standard forms not projectively equal for inverses\n%s" %(oo)

            if o_std != oo_std:
                raise AssertionError, "no equality between relabeled\n%s\n- - - - - - - -\n%s" %(o,oo)

    def test_cylcoord(self, o):
        r"""
        Test cylinder coordinates
        """
        oo = o.to_standard_form()
        c,l,h,t = o.cylinder_diagram(data = True)
        ooo = c.cylcoord_to_origami(l,h,t)
        ooo._set_standard_form()
        if oo != ooo:
            raise AssertionError, "cylcoord do not give self back on \n%s" %str(o)

    def test_PSL2Z_action(self, o):
        r"""
        Tests that the projective action of S is such that S^2=id and that the S
        action is such that S^2=-id
        """
        oo = o.S_action().S_action()
        o_inv = o.inverse()

        o_inv_std_proj = o_inv.projective_standard_form()
        o_std_proj = o.projective_standard_form()
        oo_std_proj = oo.projective_standard_form()

        if o_inv_std_proj != o_std_proj:
            raise AssertionError, "-id does not act projectively as identity on \n%s" %str(o)

        if o_std_proj != oo_std_proj:
            raise AssertionError, "S^2 does not act projectively as identity on\n%s" %str(o)

        oo = o.S_action().L_action().S_action().L_action().S_action().L_action()
        oo_std_proj = oo.projective_standard_form()

        if o_std_proj != oo_std_proj:
            raise AssertionError, "(SL)^3 does not act projectively as identity on\n%s" %str(o)

        oo = o.S_action().R_action().S_action().R_action().S_action().R_action()
        oo_std_proj = oo.projective_standard_form()

        if o_std_proj != oo_std_proj:
            raise AssertionError, "(SR)^3 does not act projectively as identity on\n%s" %str(o)

    def test_SL2Z_action(self, o):
        r"""
        Tests that the projective action of S is such that S^2=id and that the S
        action is such that S^2=-id
        """
        o_std = o.to_standard_form()
        o_inv = o.inverse()
        o_inv_std = o_inv.to_standard_form()

        # test S^2
        oo = o.S_action().S_action()
        oo_std = oo.to_standard_form()

        if o_inv_std != oo_std:
            raise AssertionError, "S^2 does not act as -id on \n%s" %str(o)

        # test (SL)^3 = -1
        oo = o.S_action().L_action().S_action().L_action().S_action().L_action()
        oo_std = oo.to_standard_form()

        if o_inv_std != oo_std:
            raise AssertionError, "(SL)^3 does not act as -id on \n%s" %str(o)

        # test (LS)^3 = -1
        oo = o.L_action().S_action().L_action().S_action().L_action().S_action()
        oo_std = oo.to_standard_form()

        if o_inv_std != oo_std:
            raise AssertionError, "(LS)^3 does not act as -id on \n%s" %str(o)

        # test (SR)^3 = 1
        oo = o.S_action().R_action().S_action().R_action().S_action().R_action()
        oo_std = oo.to_standard_form()

        if o_std != oo_std:
            raise AssertionError, "(SR)^3 does not act as id on \n%s" %str(o)

        # test (RS)^3 = 1
        oo = o.R_action().S_action().R_action().S_action().R_action().S_action()
        oo_std = oo.to_standard_form()

        if o_std != oo_std:
            raise AssertionError, "(RS)^3 does not act as -id on \n%s" %str(o)

    def test_simple_closed_curves(self, o):
        edges = OrigamiEdges(o)
        Z = edges.cycle_space()
        basis = edges.basis_of_simple_closed_curves()

        if len(basis) != Z.dimension():
            raise AssertionError, "size error %d simple closed curve for dim(Z1) = %d" %(
                    len(basis), Z.dimension())

        for v in basis:
            if v not in Z:
                raise AssertionError, "vector %s not in Z1" %str(v)

        for i in xrange(len(basis)):
            v = basis[i]
            if edges.intersection(v,v) != 0:
                raise AssertionError, "auto intersection for %s" %str(v)
            for j in xrange(i):
                w = basis[j]
                i_vw = edges.intersection(v,w)
                i_wv = edges.intersection(w,v)
                if i_vw != i_wv:
                    raise AssertionError, "anti symmetry error for\nv=%s\nw=%s\ni(v,w)=%d and i(w,v)=%d" %(str(v),str(w),i_vw,i_wv)

    def test_chain_complex_dimensions(self, o):
        r"""
        Test the dimensions of the (reduced) chain complex of the origami o
        """
        CC = o.chain_complex()
        n = o.nb_squares()
        nv = o.nb_vertices(True)
        
        dims = [(1,1,1), (nv,nv-1,nv-1), (2*n,n-1,2*n-nv+1), (n,0,1)]

        for i in xrange(4):
            C = CC.chain_space(i-1)
            B = CC.border_space(i-1)
            Z = CC.cycle_space(i-1)
            dim = dims[i]

            if not B.is_subspace(Z) or not Z.is_subspace(C):
                raise AssertionError, "bad inclusion in degree %d" %(i-1)

            if C.dimension() != dim[0] or B.dimension() != dim[1] or Z.dimension() != dim[2]:
                raise AssertionError, "wrong dimensions in degree %d dim(C) = %d dim(B) = %d dim(Z) = %d for\n%s" %(
                          i-1, C.dimension(), B.dimension(), Z.dimension(), str(o))

    def test_vertices(self, o):
        r"""
        Test the function of vertices on o
        """
        for v in o.vertices(True):
            edge_in, edge_out = v.adjacent_edge_indices()
            pos_in, pos_out = v.edge_positions()
            for i in xrange(len(edge_in)):
                j = 2*i - i%2
                if pos_in[edge_in[i]] != j:
                    raise AssertionError, "bad pos for edge %d" %edge_in[i]
                if pos_out[edge_out[i]] != j+2:
                    raise AssertionError, "bad pos for edge %d" %edge_out[i]


