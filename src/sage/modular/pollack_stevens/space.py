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
