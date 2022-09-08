from sage.modules.free_quadratic_module_integer_symmetric import FreeQuadraticModule_integer_symmetric

class LatticeWithIsometry(FreeQuadraticModule_integer_symmetric):
    r"""
    """

    def __init__(self, Lattice, iso, order=None, gramE=None, magmaRep=None, basisE=None):
        r"""
        """
        self.L = Lattice
        assert True==isinstance(Lattice,FreeQuadraticModule_integer_symmetric)
        self.iso = iso#Lattice.orthogonal_group([iso.change_ring(ZZ)]).gen(0)
        if order is None:
            order = iso.multiplicative_order()
        self.order = order
        self._gramE = gramE
        self._magmaRep = magmaRep
        self._Oq_equiv = None
        self._basisE = basisE


    def __repr__(self):
        return "Lattice with an isometry of order " + str(self.order)+".\n"+self.L.__repr__()

    def lattice(self):
        r"""
        """
        return self.L

    def isometry(self):
        r"""
        """
        return self.iso

    def sublattice_from_poly(self, poly):
        r"""
        """
        R.<x> = ZZ[]
        poly = R(poly)
        L_new = self._L.sublattice(poly(self._iso.matrix().change_ring(ZZ)).kernel().basis_matrix())
        iso_new = self._iso.matrix()
        return LatticeWithIsometry(L_new, iso_new)

    def Oq_equiv(self):
        r"""
        Image of the equivariant part of the orthogonal group.
        """
        if self._Oq_equiv is not None:
            return self._Oq_equiv
        if self.order==1 or self.order==2:
            return self.L.image_in_Oq()
        elif prod(self.L.signature_pair()) == 0:
            # E = self._gramE.base_ring()
            # gensG = self._magmaRep.AutomorphismGroup().Generators()
            # gensG = [g.ChangeRing(E).sage() for g in gensG]
            # a = E.gen().matrix()
            # gensG = [block_matrix(g.nrows(),g.ncols(),[i.polynomial().subs(x=a) for i in g.list() ]) for g in gensG]
            # G = self.L.orthogonal_group(gens=gensG)
            # Oq_equiv = Oq.subgroup([Oq(g) for g in G.gens()])
            if self.order==13:
                G = [self.L.orthogonal_group(gens=[self.iso]).gen(0)]
            else:
                G = self.L.orthogonal_group()
                f = G(self.iso)
                G = [G(g) for g in G.gap().Centralizer(f.gap()).GeneratorsOfGroup()]
            Oq = self.L.discriminant_group().orthogonal_group()
            Oq_equiv = Oq.subgroup([Oq(g) for g in G])
        elif self.L.rank() == euler_phi(self.order):
            O_equiv = self.L.O([-self.iso**0, self.iso])
            Oq_equiv = self.L.q().O(O_equiv.gens())
        else:
            # print(self.L, self.iso, self._gramE,self._magmaRep)
            Oq_equiv = Oq_equiv_herm(self.L, self.iso, self._gramE, self._magmaRep)
        self._Oq_equiv = Oq_equiv
        #self._gramE = None
        #self._magmaRep = None
        return Oq_equiv
