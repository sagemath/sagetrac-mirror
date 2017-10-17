"""
A standardized interface to a homotopy continuation solver.
"""

from abc import ABCMeta

# TODO: Correct classes

class HomotopyContinuationEngine(SageObject):
    def __init__(self, numThreads):
        raise NotImplementedError

    def mixed_volume(self, polynomialSystem):
        raise NotImplementedError

    def track_paths(self, homotopy, parameterStartValue, \
        parameterEndValue, startSolutions):
        raise NotImplementedError

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
    def __init__(self,numThreads):
        try:
            import phcpy
        except ImportError:
            raise RuntimeError("phcpy not found")
        self.numThreads = numThreads
        return None

    def __repr__(self): return phcpy.py2c_PHCpack_version_string()
    def __str__(self): return self.__repr__()

    def syst_to_phcpy_strs(self, polynomialSystem):
        pass

    def phcpy_sols_to_NomericalPoints(self, phcpySols):
        pass

    def call_phcpy_function(func):
        try: return func()
        except Exception: raise RuntimeError("Error running phcpy")

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

        f = lambda: phcpy.solver.solve(syst_to_phcpy_strs(system), verbose = False, precision = phcpyPrec)
        sols = call_phcpy_function(f)
        return phcpy_sols_to_NumericalPoints(sols)

    def mixed_volume(self, system):
        f = lambda: phcpy.solver.mixed_volume(syst_to_phcpy_strs(system))
        return Integer(call_phcpy_function(f))
