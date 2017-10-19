"""
A standardized interface to a homotopy continuation solver.
"""

from abc import abstractmethod, ABCMeta
from sage.rings.all import Integer
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.misc.flatten import flatten
from sage.structure.all import SageObject


class HomotopyContinuationEngine(SageObject):
    __metaclass__ = ABCMeta

    @abstractmethod
    def __init__(self, numthreads=1, prec=53, digits=None, useadaptiveprec=False):
        self.__numthreads = numthreads

    def mixed_volume(self, polynomialsystem, stable=False):
        raise NotImplementedError
    def bezout_bound(self, polynomialsystem):
        raise NotImplementedError

    @abstractmethod
    def numerical_irreducible_decomp(self, polynomialsystem, topdim=None):
        raise NotImplementedError

    @abstractmethod
    def track_paths(self, homotopy, parameterstartvalue, \
        parameterendvalue, startsolutions):
        raise NotImplementedError

    @abstractmethod
    def zero_dim_solve(self, system):
        """
        Assumes the system is 0-dimensional.
        INPUT:
        - ``prec`` -- precision in bits
        - ``digits`` -- precision in decimal digits (only used if
          ``prec`` is not given)
        """
        raise NotImplementedError

class PHCpackEngine(HomotopyContinuationEngine):
    def __init__(self, numthreads=1, prec=53, digits=None, useadaptiveprec=False):
        HomotopyContinuationEngine.__init__(self)
        try:
            import phcpy
            self.phcpy = phcpy
        except ImportError:
            raise RuntimeError("phcpy not found")
        try:
            assert numthreads >= 1
        except Exception:
            raise ValueError("Invalid thread count: " + str(numthreads))
        # __phcpyvartoinputvardict records the correspondence between the
        # original variables given by the user, and our own self-created
        # variables, which we use to avoid the issue of phcpy not allowing
        # certain variable names such as e, i, etc. Also, makes the strsol2dict
        # work if the user chose m, t, etc. as a variable name.
        self.__phcpyvartoinputvardict = {}
        self.__originalring = None

        self.set_precision(prec, digits, useadaptiveprec)

    def __repr__(self):
        return "Interface to " + self.phcpy.py2c_PHCpack_version_string()
    def __str__(self):
        return self.__repr__()

    def __syst_to_phcpy_strs(self, polynomialsystem):
        inputring = polynomialsystem.ring
        self.__originalring = inputring
        genindexlist = range(inputring.ngens())
        newvarnames = ['x'+str(i) for i in genindexlist]
        self.__phcpyvartoinputvardict = {newvarnames[i]: inputring.gen(i) for i in genindexlist}
        newring = PolynomialRing(inputring.base_ring(), names=newvarnames)
        subbedpolys = [
            p.subs({inputring.gen(i): newring.gen(i) for i in genindexlist})
            for p in polynomialsystem.polys]
        return [str(p) + ';' for p in subbedpolys]

    def __phcpy_sols_to_numerical_pts(self, phcpysols):
        toreturn = []
        basering = self.__originalring.base_ring()
        for phcpysol in phcpysols:
            soldict = self.phcpy.solutions.strsol2dict(phcpysol)
            pointdict = {self.__phcpyvartoinputvardict[key] : basering(soldict[key])
                         for key in self.__phcpyvartoinputvardict.keys()}
            toreturn.append(NumericalPoint(pointdict,
                                           ring=self.__originalring,
                                           multiplicity=soldict['m'],
                                           rco=soldict['rco'],
                                           err=soldict['err'],
                                           res=soldict['res']))
        return toreturn

    @classmethod
    def __call_phcpy_function(cls, func):
        try:
            return func()
        except Exception:
            raise RuntimeError("Error running phcpy")

    def set_precision(self, prec=53, digits=None, useadaptiveprec=False):
        r"""
        Set the precision of PHCpack computations done with this engine object.

        Rounds up to the nearest precision of double, double-double, or quad-double.

        INPUT:
        - ``prec`` -- precision in bits
        - ``digits`` -- precision in decimal digits (only used if
          ``prec`` is not given)
        - ``useadaptiveprec`` -- whether to use adaptive multiprecision. Currently not supported.
        """
        if useadaptiveprec:
            raise NotImplementedError("Adaptive precision not supported")

        if prec is not None:
            if prec <= 53:
                phcpyprec = 'd'
            elif prec <= 106:
                phcpyprec = 'dd'
            elif prec <= 212:
                phcpyprec = 'qd'
            else:
                raise ValueError(
                    "phcpy cannot work in greater than 212 bits of precision.")
        elif digits is not None:
            if digits <= 16:
                phcpyprec = 'd'
            elif digits <= 32:
                phcpyprec = 'dd'
            elif digits <= 64:
                phcpyprec = 'qd'
            else:
                raise ValueError(
                    "phcpy cannot work in greater than 64 digits of precision.")
        else:
            phcpyprec = 'd'

        self.__phcpy_prec = phcpyprec

    def mixed_volume(self, system, stable=False):
        func = lambda: Integer(self.phcpy.solver.mixed_volume(
            self.__syst_to_phcpy_strs(system), stable))
        return self.__call_phcpy_function(func)

    def bezout_bound(self, polynomialSystem):
        func = lambda: Integer(self.phcpy.solver.total_degree(
            self.__syst_to_phcpy_strs(polynomialSystem)))
        return self.__call_phcpy_function(func)

    def numerical_irreducible_decomp(self, polynomialSystem, topDim=None):
        numvars = polynomialSystem.ring.nvars()
        if topDim is None:
            topDim = numvars - 1
        exponents = flatten([p.exponents(as_ETuples=False)
                             for p in polynomialSystem.polys])

        islaurent = False
        for exponent in exponents:
            if exponent < 0:
                islaurent = True
                break

        self.phcpy.factor.solve(numvars,
                                dim=topDim,
                                pols=polynomialSystem.polys,
                                islaurent=islaurent,
                                precision=self.__phcpy_prec,
                                tasks=self.__numthreads - 1,
                                verbose=self.__verbose)

    def track_paths(self, homotopy, parameterStartValue, parameterEndValue, startSolutions):
        pass

    def zero_dim_solve(self, polynomialSystem):
        # assumes 0-dim, doesn't check
        syststrs = self.__syst_to_phcpy_strs(polynomialSystem)
        func = lambda: self.phcpy.solver.solve(syststrs,
                                               verbose=self.__verbose,
                                               tasks=self.__numthreads - 1,
                                               precision=self.__phcpyprec)
        sols = self.__call_phcpy_function(func)
        return self.__phcpy_sols_to_numerical_pts(sols)

