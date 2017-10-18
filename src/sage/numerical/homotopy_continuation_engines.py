"""
A standardized interface to a homotopy continuation solver.
"""

import abc

# TODO: Correct classes

class HomotopyContinuationEngine(SageObject):
    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def __init__(self, numThreads=1):
        raise NotImplementedError

    def mixed_volume(self, polynomialSystem, stable=False):
        raise NotImplementedError
    def bezout_bound(self, polynomialSystem):
        raise NotImplementedError

    @abc.abstractmethod
    def track_paths(self, homotopy, parameterStartValue, \
        parameterEndValue, startSolutions):
        raise NotImplementedError

    @abc.abstractmethod
    def zero_dim_solve(self, system, prec=None, digits=None, adaptivePrec=False):
        """
        Assumes the system is 0-dimensional.
        INPUT:
        - ``prec`` -- precision in bits
        - ``digits`` -- precision in decimal digits (only used if
          ``prec`` is not given)
        """
        raise NotImplementedError

class PHCpackEngine(HomotopyContinuationEngine):
    def __init__(self, numThreads=1):
        try:
            import phcpy
            self.phcpy = phcpy
        except ImportError:
            raise RuntimeError("phcpy not found")
        try: assert numThreads >= 1
        except Exception: raise ValueError("Invalid thread count: " + str(numThreads))
        self.numThreads = numThreads
        # This records the correspondence between the original variables given by the user,
        # and our own self-created variables, which we use to avoid the issue of phcpy not
        # allowing certain variable names such as e, i, etc. Also, makes the strsol2dict
        # work if the user chose m, t, etc. as a variable name.
        self.__phcpyVarToInputVarDict = {}

    def mixed_volume(self, system, stable=False):
        f = lambda: Integer(self.phcpy.solver.mixed_volume(self.__syst_to_phcpy_strs(system), stable))
        return self.__call_phcpy_function(f)
    def bezout_bound(self, polynomialSystem):
        f = lambda: Integer(self.phcpy.solver.total_degree(self.__syst_to_phcpy_strs(polynomialSystem)))
        return self.__call_phcpy_function(f)

    def __repr__(self): return "Interface to " + self.phcpy.py2c_PHCpack_version_string()
    def __str__(self): return self.__repr__()

    def __syst_to_phcpy_strs(self, polynomialSystem):
        self.__phcpyVarToInputVarDict
        return map(lambda p:str(p)+';', polynomialSystem.polys)

    def __phcpy_sols_to_NumericalPoints(self, phcpySols, ring):
        toReturn = []
        for phcpySol in phcpySols:
            solDict = self.phcpy.solutions.strsol2dict(phcpySol)
            pointDict = {variable:solDict[str(variable)] for variable in ring.gens()}
            #NumericalPoint(pointDict, ring=ring, multiplicity = None, condition_number=None):

    def __call_phcpy_function(self, func):
        try: return func()
        except Exception: raise RuntimeError("Error running phcpy")

    def track_paths(self, homotopy, parameterStartValue, \
        parameterEndValue, startSolutions):
        pass

    def zero_dim_solve(self, system, prec=None, digits=None, adaptivePrec=False):
        if adaptivePrec: raise NotImplementedError("adaptive precision not currently supported")

        if not (prec is None):
            if prec <= 53: phcpyPrec = 'd'
            elif prec <= 106: phcpyPrec = 'dd'
            elif prec <= 212: phcpyPrec = 'qd'
            else: raise ValueError("Cannot work in greater than 212 bits of precision with PHCpack")
        elif not (digits is None):
            if digits <= 16: phcpyPrec = 'd'
            elif digits <= 32: phcpyPrec = 'dd'
            elif digits <= 64: phcpyPrec = 'qd'
            else: raise ValueError("Cannot work in greater than 64 digits of precision with PHCpack")
        else: # prec, digits both none
            phcpyPrec = 'd'
        systStrs = self.__syst_to_phcpy_strs(system)
        f = lambda: self.phcpy.solver.solve(systStrs, verbose=False, tasks=self.numThreads - 1, precision=phcpyPrec)
        sols = self.__call_phcpy_function(f)
        return self.__phcpy_sols_to_NumericalPoints(sols)


