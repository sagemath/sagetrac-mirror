"""
A standardized interface to a homotopy continuation solver.
"""

from abc import abstractmethod, ABCMeta
from sage.rings.all import Integer
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.misc.flatten import flatten
from sage.structure.all import SageObject
from sage.numerical.polynomial_homotopy_types import NumericalPoint

class HomotopyContinuationEngine(SageObject):
    """
    A base class for homotopy continuation engines. Any engine should inherit
    from this class.
    """
    __metaclass__ = ABCMeta

    @abstractmethod
    def __init__(self, numthreads=1, prec=53, digits=None, useadaptiveprec=False):
        pass

    @abstractmethod
    def __repr__(self):
        return "Homotopy engine base class instance."

    @abstractmethod
    def __str__(self):
        return self.__repr__()

    @abstractmethod
    def mixed_volume(self, polynomial_system, stable=False):
        """
        Compute mixed volume
        """
        raise NotImplementedError

    @abstractmethod
    def bezout_bound(self, polynomial_system):
        """
        Compute Bezout bound
        """
        raise NotImplementedError

    @abstractmethod
    def numerical_irreducible_decomp(self, polynomial_system, topdim=None):
        """
        Find a numerical irreducible decomposition
        """
        raise NotImplementedError

    @abstractmethod
    def track_paths(self, homotopy, parameter_start_value, \
        parameter_end_value, startsolutions):
        """
        Track along a homotopy
        """
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
    """
    PHCpack engine class for performing homotopy continuation.
    """
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
        # __phcpy_var_to_input_var_dict records the correspondence between the
        # original variables given by the user, and our own self-created
        # variables, which we use to avoid the issue of phcpy not allowing
        # certain variable names such as e, i, etc. Also, makes the strsol2dict
        # work if the user chose m, t, etc. as a variable name.
        self.__phcpy_var_to_input_var_dict = {}
        self.__original_ring = None
        self.__numthreads = numthreads

        self.set_precision(prec, digits, useadaptiveprec)

    def __repr__(self):
        return "Interface to " + self.phcpy.py2c_PHCpack_version_string()

    def __str__(self):
        return self.__repr__()

    def __syst_to_phcpy_strs(self, polynomial_system):
        input_ring = polynomial_system.ring()
        self.__original_ring = input_ring
        gen_index_list = range(input_ring.ngens())
        new_var_names = ['x'+str(i) for i in gen_index_list]
        self.__phcpy_var_to_input_var_dict = {new_var_names[i]: input_ring.gen(i) for i in gen_index_list}
        new_ring = PolynomialRing(input_ring.base_ring(), names=new_var_names)
        subbed_polys = [
            p.subs({input_ring.gen(i): new_ring.gen(i) for i in gen_index_list})
            for p in polynomial_system.polynomials()]
        return [str(p) + ';' for p in subbed_polys]

    def __phcpy_sols_to_numerical_pts(self, phcpy_sols):
        to_return = []
        base_ring = self.__original_ring.base_ring()
        for phcpy_sol in phcpy_sols:
            sol_dict = self.phcpy.solutions().strsol2dict(phcpy_sol)
            try:
                point_dict = {self.__phcpy_var_to_input_var_dict[key] : base_ring(sol_dict[key])
                             for key in self.__phcpy_var_to_input_var_dict.keys()}
            except TypeError:
                base_ring = ComplexField(prec = self.__prec)
                point_dict = {self.__phcpy_var_to_input_var_dict[key] : base_ring(sol_dict[key])
                             for key in self.__phcpy_var_to_input_var_dict.keys()}
            to_return.append(NumericalPoint(point_dict,
                                           ring=self.__original_ring,
                                           multiplicity=sol_dict['m'],
                                           rco=sol_dict['rco'],
                                           err=sol_dict['err'],
                                           res=sol_dict['res']))
        return to_return

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
                phcpy_prec = 'd'
            elif prec <= 106:
                phcpy_prec = 'dd'
            elif prec <= 212:
                phcpy_prec = 'qd'
            else:
                raise ValueError(
                    "phcpy cannot work in greater than 212 bits of precision.")
        elif digits is not None:
            if digits <= 16:
                phcpy_prec = 'd'
            elif digits <= 32:
                phcpy_prec = 'dd'
            elif digits <= 64:
                phcpy_prec = 'qd'
            else:
                raise ValueError(
                    "phcpy cannot work in greater than 64 digits of precision.")
        else:
            phcpy_prec = 'd'

        if phcpy_prec == 'd':
            self.__digits = 16
            self.__prec = 53
        elif phcpy_prec == 'dd':
            self.__digits = 32
            self.__prec = 106
        elif phcpy_prec == 'qd':
            self.__digits = 64
            self.__prec = 212
        else:
            raise ValueError("Invalid phcpy precision type: %s"%phcpy_prec)
            
        self.__phcpy_prec = phcpy_prec

    def mixed_volume(self, system, stable=False):
        func = lambda: Integer(self.phcpy.solver.mixed_volume(
            self.__syst_to_phcpy_strs(system), stable))
        return self.__call_phcpy_function(func)

    def bezout_bound(self, polynomialSystem):
        func = lambda: Integer(self.phcpy.solver.total_degree(
            self.__syst_to_phcpy_strs(polynomialSystem)))
        return self.__call_phcpy_function(func)

    def numerical_irreducible_decomp(self, polynomialSystem, topDim=None):
        numvars = polynomialSystem.ring().nvars()
        if topDim is None:
            topDim = numvars - 1
        exponents = flatten([p.exponents(as_ETuples=False)
                             for p in polynomialSystem.polynomials()])

        islaurent = False
        for exponent in exponents:
            if exponent < 0:
                islaurent = True
                break

        sols = self.phcpy.factor.solve(numvars,
                                       dim=topDim,
                                       pols=polynomialSystem.polynomials(),
                                       islaurent=islaurent,
                                       precision=self.__phcpy_prec,
                                       tasks=self.__numthreads - 1,
                                       verbose=False)
        print sols

    def track_paths(self, homotopy, parameterStartValue, parameterEndValue, startSolutions):
        pass

    def zero_dim_solve(self, polynomialSystem):
        # assumes 0-dim, doesn't check
        syststrs = self.__syst_to_phcpy_strs(polynomialSystem)
        func = lambda: self.phcpy.solver.solve(syststrs,
                                               verbose=False,
                                               tasks=self.__numthreads - 1,
                                               precision=self.__phcpy_prec)
        sols = self.__call_phcpy_function(func)
        return self.__phcpy_sols_to_numerical_pts(sols)

