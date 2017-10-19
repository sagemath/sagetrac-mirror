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
    def zero_dim_solve(self, system, prec=None, digits=None, useAdaptivePrec=False):
        """
        Assumes the system is 0-dimensional.
        INPUT:
        - ``prec`` -- precision in bits
        - ``digits`` -- precision in decimal digits (only used if
          ``prec`` is not given)
        """
        raise NotImplementedError

class PHCpackEngine(HomotopyContinuationEngine):
    def __init__(self, numThreads=1, prec=53, digits=None, useAdaptivePrec=False, verbose=False):
        try:
            import phcpy
            self.phcpy = phcpy
        except ImportError:
            raise RuntimeError("phcpy not found")
        try: assert numThreads >= 1
        except Exception: raise ValueError("Invalid thread count: " + str(numThreads))
        self.__numThreads = numThreads
        # __phcpyVarToInputVarDict records the correspondence between the
        # original variables given by the user, and our own self-created
        # variables, which we use to avoid the issue of phcpy not allowing
        # certain variable names such as e, i, etc. Also, makes the strsol2dict
        # work if the user chose m, t, etc. as a variable name.
        self.__phcpyVarToInputVarDict = {}
        self.__originalRing = None
        self.__verbose = verbose
        self.set_prec(prec, digits, useAdaptivePrec)

    def __repr__(self): return "Interface to " + self.phcpy.py2c_PHCpack_version_string()
    def __str__(self): return self.__repr__()

    def __syst_to_phcpy_strs(self, polynomialSystem):
        inputRing = polynomialSystem.ring
        self.__originalRing = inputRing
        genIndexList = range(inputRing.ngens())
        newVarNames = ['x'+str(i) for i in genIndexList]
        self.__phcpyVarToInputVarDict = {newVarNames[i]: inputRing.gen(i) for i in genIndexList}
        newRing = PolynomialRing(inputRing.base_ring(), names=newVarNames)
        subbedPolys = [p.subs({inputRing.gen(i): newRing.gen(i) for i in genIndexList}) for p in polynomialSystem.polys]
        return map(lambda p:str(p)+';', subbedPolys)

    def __phcpy_sols_to_NumericalPoints(self, phcpySols):
        toReturn = []
        baseRing = self.__originalRing.base_ring()
        for phcpySol in phcpySols:
            solDict = self.phcpy.solutions.strsol2dict(phcpySol)
            pointDict = {self.__phcpyVarToInputVarDict[key] : baseRing(solDict[key]) for key in self.__phcpyVarToInputVarDict.keys()}
            newPt = NumericalPoint(pointDict, ring=self.__originalRing, multiplicity=solDict['m'], rco=solDict['rco'], err=solDict['err'], res=solDict['res'])
            toReturn.append(newPt)
        return toReturn

    def __call_phcpy_function(self, func):
        try: return func()
        except Exception: raise RuntimeError("Error running phcpy")
        
    def set_precision(prec=53, digits=None, useAdaptivePrec=False):
        r"""
        Set the precision of PHCpack computations done with this engine object.

        Rounds up to the nearest precision of double, double-double, or quad-double.

        INPUT:
        - ``prec`` -- precision in bits
        - ``digits`` -- precision in decimal digits (only used if
          ``prec`` is not given)
        - ``useAdaptivePrec`` -- whether to use adaptive multiprecision. Currently not supported.
        """
        if useAdaptivePrec: raise NotImplementedError("adaptive precision not currently supported")
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
        self.__phcpy_prec = phcpyPrec

    def mixed_volume(self, system, stable=False):
        f = lambda: Integer(self.phcpy.solver.mixed_volume(self.__syst_to_phcpy_strs(system), stable))
        return self.__call_phcpy_function(f)
    def bezout_bound(self, polynomialSystem):
        f = lambda: Integer(self.phcpy.solver.total_degree(self.__syst_to_phcpy_strs(polynomialSystem)))
        return self.__call_phcpy_function(f)


    def track_paths(self, homotopy, parameterStartValue, \
        parameterEndValue, startSolutions):
        pass

    def zero_dim_solve(self, system):
        # assumes 0-dim, doesn't check
        systStrs = self.__syst_to_phcpy_strs(system)
        f = lambda: self.phcpy.solver.solve(systStrs, verbose=self.__verbose, tasks=self.__numThreads - 1, precision=self.__phcpyPrec)
        sols = self.__call_phcpy_function(f)
        return self.__phcpy_sols_to_NumericalPoints(sols)

