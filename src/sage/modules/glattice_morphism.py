from __future__ import print_function, absolute_import

import sage.modules.glattice
from sage.modules.free_module_morphism import FreeModuleMorphism
from sage.structure.element import is_Matrix
from inspect import isfunction
from . import matrix_morphism
from sage.matrix.constructor import matrix
from sage.matrix.matrix_space import MatrixSpace

class GLatticeMorphism(FreeModuleMorphism):
    def _repr_(self):
        r = "Lattice morphism defined by the matrix\n{!r}\nDomain: {}\nCodomain: {}"
        return r.format(self.matrix(), self.domain(), self.codomain())

    def kernel(self):
    	K = FreeModuleMorphism.kernel(self)
    	return self.domain().sublattice(K.basis())


class GLatticeHomspace(sage.categories.homset.HomsetWithBase):
    def __call__(self, A, check=True):
        from . import free_module_morphism
        C = self.codomain()
        D = self.domain()
        if not is_Matrix(A):
            # Compute the matrix of the morphism that sends the
            # generators of the domain to the elements of A.
            try:
                if isfunction(A):
                    v = [C(A(g)) for g in self.domain().gens()]
                else:
                    v = [C(a) for a in A]
                A = matrix([C.coordinates(a) for a in v], ncols=C.rank())
            except TypeError:
                # Let us hope that FreeModuleMorphism knows to handle
                # that case
                pass	
        for e in D.gens():
            for g in D.group().gens():
                if A*(D._action(g)(e)) != C._action(g)(A*e):
                    raise TypeError("The morphism does not preserve the action of the group")
        return GLatticeMorphism(self, A)

    def _matrix_space(self):
        try:
            return self.__matrix_space
        except AttributeError:
            R = self.codomain().base_ring()
            M = MatrixSpace(R, self.domain().rank(), self.codomain().rank())
            self.__matrix_space = M
            return M
