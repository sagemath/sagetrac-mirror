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

    @abc.abstractmethod
    def mixed_volume(self, polynomialSystem):
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

    def __repr__(self): return "Interface to " + self.phcpy.py2c_PHCpack_version_string()
    def __str__(self): return self.__repr__()

    def __syst_to_phcpy_strs(self, polynomialSystem):
        return map(lambda p:str(p)+';', polynomialSystem.polys)

#def __init__(self, coords, ring=None, multiplicity=None, condition_number=None):
    def __phcpy_sols_to_NumericalPoints(self, phcpySols):
        pass
        """
        import phcpy.solutions
        toReturn = []
        for phcpySol in phcpySols:
            solDict = self.phcpy.solutions.strsol2dict(phcpySol)
            NumericalPoint(
        """

    def __call_phcpy_function(self, func):
        try: return func()
        except Exception: raise RuntimeError("Error running phcpy")

    def track_paths(self, homotopy, parameterStartValue, \
        parameterEndValue, startSolutions):
        pass

    def zero_dim_solve(self, system, prec=None, digits=None, adaptivePrec=False):
        import phcpy.solver ##### kill w/ new version of phcpy
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

        f = lambda: self.phcpy.solver.solve(self.__syst_to_phcpy_strs(system), verbose = False, tasks = self.numThreads - 1, precision = phcpyPrec)
        sols = self.__call_phcpy_function(f)
        return self.__phcpy_sols_to_NumericalPoints(sols)

    def mixed_volume(self, system):
        import phcpy.solver ##### kill w/ new version of phcpy
        f = lambda: self.phcpy.solver.mixed_volume(self.__syst_to_phcpy_strs(system))
        return Integer(self.__call_phcpy_function(f))
