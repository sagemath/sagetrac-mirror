from sage.modules.module import Module
from sage.structure.factory import UniqueFactory
from distributions import Distributions
from sage.modular.dirichlet import DirichletCharacter
from sage.modular.arithgroup.all import Gamma0
from sage.rings.integer import Integer
from sage.rings.rational_field import QQ
from modsym import ModularSymbolElement, ModSymAction
from fund_domain import manin_relations

class PSModularSymbols_constructor(UniqueFactory):
    def create_key(self, group=None, weight=None, sign=0, base_ring=None, p=None, prec_cap=None, coefficients=None):
        if isinstance(group, (int, Integer)):
            group = Gamma0(group)
        if base_ring is None and p is None:
            base_ring = QQ
        if coefficients is None:
            if p is not None and prec_cap is None:
                prec_cap = 20
            if isinstance(group, DirichletCharacter):
                character = group.minimize_base_ring()
                group = Gamma0(character.modulus())
                character = (character, None)
            else:
                character = None
            if weight is None: raise ValueError("you must specify a weight or coefficient module")
            k = weight - 2
            coefficients = Distributions(k, p, prec_cap, base_ring, character)
        return (group, coefficients, sign)

    def create_object(self, version, key):
        return PSModularSymbolSpace(*key)

PSModularSymbols = PSModularSymbols_constructor('PSModularSymbols')

class PSModularSymbolSpace(Module):
    """
    A class for spaces of modular symbols that use Glenn Stevens'
    conventions.

    There are two main differences between the modular symbols in this
    directory and the ones in sage.modular.modsym.

    - There is a shift in the weight: weight `k=0` here corresponds to
      weight `k=2` there.

    - There is a duality: these modular symbols are functions from
      `Div^0(P^1(\QQ))`, the others are formal linear combinations of
      such elements.

    INPUT:

    - ``V`` -- the coefficient module, which should have a right action of `M_2(\ZZ)`

    - ``domain`` -- a set or None, giving the domain 
    """
    Element = ModularSymbolElement
    def __init__(self, group, coefficients, sign):
        self._group = group
        self._coefficients = coefficients
        self._sign = sign
        self._manin_relations = manin_relations(group.level()) # should distingish between Gamma0 and Gamma1...
        act = ModSymAction(self)
        self._populate_coercion_lists_(action_list=[act])

    def ngens(self):
        r"""
        Returns the number of generators defining self.

        EXAMPLES:

        ::
        sage: E = EllipticCurve('11a')
        sage: from sage.modular.overconvergent.pollack.modsym_symk import form_modsym_from_elliptic_curve
        sage: phi = form_modsym_from_elliptic_curve(E); phi
        [-1/5, 3/2, -1/2]
        sage: phi.ngens()
        3
        sage: E = EllipticCurve('37a')
        sage: from sage.modular.overconvergent.pollack.modsym_symk import form_modsym_from_elliptic_curve
        sage: phi = form_modsym_from_elliptic_curve(E); phi
        [0, 1, 0, 0, 0, -1, 1, 0, 0]
        sage: phi.ngens()
        9
        """
        return len(self._manin_relations.generator_indices())

    def ncoset_reps(self):
        r"""
        Returns the number of coset representatives defining the full_data of self

        OUTPUT:
        The number of coset representatives stored in the manin relations. (Just the size
        of P^1(Z/NZ))

        EXAMPLES:

        ::

        sage: E = EllipticCurve('11a')
        sage: from sage.modular.overconvergent.pollack.modsym_symk import form_modsym_from_elliptic_curve
        sage: phi = form_modsym_from_elliptic_curve(E); phi
        [-1/5, 3/2, -1/2]
        sage: phi.ncoset_reps()
        12
        sage: E = EllipticCurve('37a')
        sage: phi = form_modsym_from_elliptic_curve(E); phi
        [0, 1, 0, 0, 0, -1, 1, 0, 0]
        sage: phi.ncoset_reps()
        38
        """
        return len(self._manin_relations.coset_reps())

    def level(self):
        r"""
        Returns the level `N` when self is of level `Gamma_0(N)`.

        INPUT:
            none

        OUTPUT:
        
        The level `N` when self is of level `Gamma_0(N)`

        EXAMPLES:

        ::

        sage: E = EllipticCurve('11a')
        sage: from sage.modular.overconvergent.pollack.modsym_symk import form_modsym_from_elliptic_curve 
        sage: phi = form_modsym_from_elliptic_curve(E); phi
        [-1/5, 3/2, -1/2]
        sage: phi.level()
        11
        """

        return self._manin_relations.level()

