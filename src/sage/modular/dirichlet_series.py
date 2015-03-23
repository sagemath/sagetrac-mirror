"""
Class for formal manipulation of Dirichlet series
=================================================
"""

#*****************************************************************************
#       Copyright (C) 2011 Jonathan Hanke
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.sage_object import SageObject
from sage.rings.infinity import Infinity
from sage.rings.integer_ring import ZZ
from sage.misc.defaults import series_precision
from sage.symbolic.expression import Expression
from sage.functions.transcendental import zeta

class DirichletSeries(SageObject):
    """
    Class for manipulation of Dirichlet series over a ring R.
    """
    def __init__(self, arg, is_exact_flag=False, default_var_string='s',
                 precision=series_precision(), base_ring=ZZ):
        """
        INPUT:

            dirichlet_series(coeff_list) -- initialize a Dirichlet series
                    with coeffs in R from a list of coefficients

            dirichlet_series(zeta(s)^2, precision=b) -- Initialize a Dirichlet series
                    from an arithmetic function, with b terms (so + O((b+1)^s))

            dirichlet_series(object) -- create a Dirichlet series from an object,
                    which should call the .dirichlet_series() method of that object.

        We will be able to set the formal variable later if we want... but for 
        now it defaults to s.  We'll need to handle this when we do substitutions...
        
        EXAMPLES::

            sage: s,t = var('s,t')
            sage: dirichlet_series(1)
            1 + O(21^(-s))
            sage: dirichlet_series(zeta(s))
            1 + 1/(2^s) + 1/(3^s) + 1/(4^s) + 1/(5^s) + 1/(6^s) + 1/(7^s) + 1/(8^s) + 1/(9^s) + 1/(10^s) + 1/(11^s) + 1/(12^s) + 1/(13^s) + 1/(14^s) + 1/(15^s) + 1/(16^s) + 1/(17^s) + 1/(18^s) + 1/(19^s) + 1/(20^s) + O(21^(-s))
            sage: dirichlet_series(zeta(t), precision=10)
            1 + 1/(2^t) + 1/(3^t) + 1/(4^t) + 1/(5^t) + 1/(6^t) + 1/(7^t) + 1/(8^t) + 1/(9^t) + 1/(10^t) + O(11^(-t))
            sage: dirichlet_series(zeta(s), precision=2)
            1 + 1/(2^s) + O(3^(-s))

        TESTS::

            sage: dirichlet_series(zeta(5))
            Traceback (most recent call last):
            ...
            ValueError: Generating function must have exactly one variable
            sage: dirichlet_series(zeta(s)*zeta(t))
            Traceback (most recent call last):
            ...
            ValueError: Generating function must have exactly one variable
        """
        if not hasattr(base_ring, 'is_ring') or not base_ring.is_ring():
            raise TypeError, "The base_ring argument must be a ring!"
        
        tmp_coeff_list = None
        defining_function = None

        if not isinstance(arg, list):
            if arg == 1:  # identity
                self.__coefficient_list = [0,1] + [0] * (precision-1)
            elif isinstance(arg, Expression):
                vars = arg.variables()
                if len(vars) != 1:
                    raise ValueError('Generating function must have exactly one variable')
                if arg.operator() == zeta and arg.operands()[0].is_symbol():
                    self.__coefficient_list = [1] * (1 + precision)
                else:
                    raise NotImplementedError

            ## Initialize the Internal variables
            self.__base_ring = base_ring
            self.__creation_function_expression = arg
            self.__number_of_coeffs = len(self.__coefficient_list) - 1
            self.__is_exact = is_exact_flag
            self.__s_variable = default_var_string

        else:
            tmp_coeff_list = [base_ring(x) for x in arg]
    
            ## Initialize the Internal variables
            self.__base_ring = base_ring
            self.__creation_function_expression = []
            self.__coefficient_list = ['X'] + tmp_coeff_list
            self.__number_of_coeffs = len(self.__coefficient_list) - 1
            self.__is_exact = is_exact_flag
            self.__s_variable = default_var_string

    def base_ring(self):
        """
        Returns the base ring of the Dirichlet series.
        """
        return self.__base_ring

        
    def precision(self):
        """
        Returns the computed precision of the Dirichlet series.
        """
        return self.__number_of_coeffs
    

    def extend_precision_to(self, n):
        """
        Extend the precision of the Dirichlet series to be O((n+1)^s).
        """
        raise NotImplementedError

    def is_Eulerian(self):
        """
        Verify that the series is Eulerian up to the given precision, and cache the result.
        """
        raise NotImplementedError

    def is_exact(self):
        """
        Return if the series is an exact (finite) Dirichlet series
        """
        return self.__is_exact

    def nonzero_coeffs(self):
        """
        Return a list of the indices n (up to the precision) where the n-th coefficient is non-zero.
        """
        nonzero_list = [n  for n in range(1, self.precision() + 1)  if self[n] != 0]
        return nonzero_list

    def __getitem__(self, n):
        """
        Return the value of the n-th coefficient.
        """
        if (n >=1) and (n <= self.precision()):                 ## This check is probably too slow!
            return self.__coefficient_list[n]
        elif self.is_exact() and (n > self.precision()):    ## Allow arbitrary evaluation precision of exact series! =)
            return self.base_ring()(0)
        else:
            raise TypeError, "The coeffecient you requested (n = " + str(n) + ") is out of the precomputed range."

    def __repr__(self, n_max=Infinity):
        """
        Print first n coefficients of the Dirichlet series, using its internally specified variable.
        """
        out_str = ""
        n_prec = min(self.precision(), n_max)
        
        ## Add the first coefficient
        if n_prec >= 1:
            out_str += str(self.__coefficient_list[1])

        ## Add all other non-zero coefficients up to the desired precision
        for i in range(2, n_prec + 1):
            if (self.__coefficient_list[i] != 0):
                out_str += " + " + str(self.__coefficient_list[i]) + "/(" + str(i) + "^" + self.__s_variable + ")"
            
        ## Add the error term, if the result is not exact.
        if not self.__is_exact:
            out_str += " + O(" + str(n_prec + 1) + "^(-" + self.__s_variable + "))"
            
        ## Return the output string
        return out_str

    def __add__(self, other):
        """
        Form the sum of two Dirichlet series.

        EXAMPLES::

            sage: D1 = dirichlet_series([1,1,1,1]);
            sage: D1 + D1
            2 + 2/(2^s) + 2/(3^s) + 2/(4^s) + O(5^(-s))
            sage: D1 + D1 + D1
            3 + 3/(2^s) + 3/(3^s) + 3/(4^s) + O(5^(-s))
            sage: s = var('s')
            sage: D2 = dirichlet_series(zeta(s), precision=10); D2
            1 + 1/(2^s) + 1/(3^s) + 1/(4^s) + 1/(5^s) + 1/(6^s) + 1/(7^s) + 1/(8^s) + 1/(9^s) + 1/(10^s) + O(11^(-s))
            sage: D1 + D2
            2 + 2/(2^s) + 2/(3^s) + 2/(4^s) + O(5^(-s))

        """
        if self.base_ring() != other.base_ring():
            raise NotImplementedError, "For now the base rings must be the same!  TO DO: Add ring coersion!"

        ## Assume that these are just given by lists of coefficients for now!
        sum_coeff_list = [self[i] + other[i]  for i in range(1, min(self.precision(), other.precision()) + 1)]
        D_sum = dirichlet_series(sum_coeff_list)

        return D_sum

    def __sub__(self, other):
        """
        Form the sum of two Dirichlet series.

        EXAMPLES::

            sage: s = var('s')
            sage: D1 = dirichlet_series([1,1,1,1]);
            sage: D1 - D1
            0 + O(5^(-s))
            sage: D2 = dirichlet_series(zeta(s), precision=10); D2
            1 + 1/(2^s) + 1/(3^s) + 1/(4^s) + 1/(5^s) + 1/(6^s) + 1/(7^s) + 1/(8^s) + 1/(9^s) + 1/(10^s) + O(11^(-s))
            sage: D1 - D2
            0 + O(5^(-s))

        """
        return self + (other * (-1))

    def __mul__(self, other):
        """
        Define the product of two Dirichlet series, or of a Dirichlet series and a number.

        EXAMPLES::

            sage: s = var('s')
            sage: D1 = dirichlet_series([1,1,1,1]);
            sage: D1 = dirichlet_series([1,1,1,1]); D1
            1 + 1/(2^s) + 1/(3^s) + 1/(4^s) + O(5^(-s))
            sage: D1 * 3
            3 + 3/(2^s) + 3/(3^s) + 3/(4^s) + O(5^(-s))
            sage: D1 = dirichlet_series([1,1,1,1]); D1
            1 + 1/(2^s) + 1/(3^s) + 1/(4^s) + O(5^(-s))
            sage: D1 * D1
            1 + 2/(2^s) + 2/(3^s) + 3/(4^s) + O(5^(-s))
            sage: D2 = dirichlet_series(zeta(s), precision=8); D2
            1 + 1/(2^s) + 1/(3^s) + 1/(4^s) + 1/(5^s) + 1/(6^s) + 1/(7^s) + 1/(8^s) + O(9^(-s))
            sage: D22 = D2 * D2
            new_coeff_list = [1, 2, 2, 3, 2, 4, 2, 4]
            sage: D22
            1 + 2/(2^s) + 2/(3^s) + 3/(4^s) + 2/(5^s) + 4/(6^s) + 2/(7^s) + 4/(8^s) + O(9^(-s))

        """
        from sage.rings.arith import divisors
        R = self.base_ring()
        try:
            ## Scalar multiplication
            scale_factor = R(other)
            return dirichlet_series([scale_factor * self[i]  for i in range(1, self.precision() + 1)])

        except:
            ## Dirichlet Multiplication
            if not self.is_exact() and not other.is_exact():
                new_precision = min(self.precision(), other.precision())
            elif self.is_exact() and other.is_exact():
                new_precision = max(self.precision(), other.precision())     ## If both are exact, we can take the precision of th
            elif self.is_exact():
                new_precision = other.precision()   
            elif other.is_exact():
                new_precision = self.precision()                
                
            new_coeff_list = [sum([self[d] * other[n/d]  for d in divisors(n)])  for n in range(1, new_precision + 1)]
            new_exact_flag = self.is_exact() and other.is_exact()
            return dirichlet_series(new_coeff_list, is_exact_flag=new_exact_flag)

    def __pow__(self, n):
        """
        Take any integer power of the current Dirichlet series.

        EXAMPLES::

            sage: s = var('s')
            sage: D2 = dirichlet_series(zeta(s), precision=8); D2
            1 + 1/(2^s) + 1/(3^s) + 1/(4^s) + 1/(5^s) + 1/(6^s) + 1/(7^s) + 1/(8^s) + O(9^(-s))
            sage: D2^2
            1 + 2/(2^s) + 2/(3^s) + 3/(4^s) + 2/(5^s) + 4/(6^s) + 2/(7^s) + 4/(8^s) + O(9^(-s))
            sage: D2^5
            1 + 5/(2^s) + 5/(3^s) + 15/(4^s) + 5/(5^s) + 25/(6^s) + 5/(7^s) + 35/(8^s) + O(9^(-s))
            sage: D2^0
            1 + 1/(2^s) + 1/(3^s) + 1/(4^s) + 1/(5^s) + 1/(6^s) + 1/(7^s) + 1/(8^s) + O(9^(-s))
            sage: D2^(-1)
            1 + -1/(2^s) + -1/(3^s) + 0/(4^s) + -1/(5^s) + 1/(6^s) + -1/(7^s) + 0/(8^s) + O(9^(-s))

        """
        R = self.base_ring()
        Identity_series = dirichlet_series(1, precision=self.precision(), is_exact_flag=True)
        
        ## SANITY CHECK: Check that n is an integer
        if not n in ZZ:
            raise TypeError, "The power must be an integer!"

        ## Identity (n=0)
        if n == 0:
            return Identity_series
        
        ## Use powers of self (n>0)
        elif n > 0:
            tmp_new_series = self
            for i in range(n-1):
                tmp_new_series = tmp_new_series * self  

        ## Use powers of the inverse (n<0)
        elif n < 0:
            inv_series = self.inverse()
            tmp_new_series = inv_series
            for i in range(abs(n)-1):
                tmp_new_series = tmp_new_series * inv_series  
                

        ## Return the series
        return tmp_new_series

    def inverse(self):
        """
        Compute the inverse Dirichlet series under Dirichlet multiplication.

        sage: s = var('s')
        sage: D2 = dirichlet_series(zeta(s), precision=8); D2
        1 + 1/(2^s) + 1/(3^s) + 1/(4^s) + 1/(5^s) + 1/(6^s) + 1/(7^s) + 1/(8^s) + O(9^(-s))
        sage: D2.inverse()
        1 + -1/(2^s) + -1/(3^s) + 0/(4^s) + -1/(5^s) + 1/(6^s) + -1/(7^s) + 0/(8^s) + O(9^(-s))
        sage: D2 * D2.inverse()
        1 + 0/(2^s) + 0/(3^s) + 0/(4^s) + 0/(5^s) + 0/(6^s) + 0/(7^s) + 0/(8^s) + O(9^(-s))
        """
        from sage.rings.arith import divisors
        R = self.base_ring()

        ## Check if the first coefficient is invertible in R
        a1 = self[1]
        a1_inv = a1^(-1)
        if not a1_inv in R:
            raise RuntimeError, "The leading term is not invertible in R, so the Dirichlet series is not invertible."
        
        extended_inv_coeff_list = ["X", a1_inv] 
        for n in range(2, self.precision() + 1):
            extended_inv_coeff_list.append(-a1_inv * sum([self[n/d] * extended_inv_coeff_list[d]  for d in divisors(n)[:-1]]))
#        print "extended_inv_coeff_list = " + str(extended_inv_coeff_list)

        return dirichlet_series(extended_inv_coeff_list[1:])

    def list(self):
        """
        Returns the list of coefficients of the Dirichlet series, where the zeroth coefficient is 'X'.
        """
        return self.__coefficient_list

    def scale_variable_by(self, a):
        """
        Scale the variable s by any positive integer a >= 1.  This takes the a-th power 
        of all coefficient indices.
        """
        ## Check that a is a positive integer
        
        ## Make the new coefficients
        old_prec = self.precision()
        new_prec = old_prec^a
        new_coeff_list = ["X"] + [0  for i in range(1, new_prec + 1)]
        for i in range(1, old_prec + 1):
            new_coeff_list[i^a] = self[i]

        ## Return the new Dirichlet series
        R = self.base_ring()
        return dirichlet_series(new_coeff_list[1:], is_exact_flag=self.is_exact())
        
    def shift_variable_by(self, a):
        """
        Shift the variable s by any non-positive integer a <= 0, or by any integer if 
        the base ring contains QQ.  This multiplies each coefficient by the a-th power 
        of its index.
        """
        ## Perform some checks

        ## Make the new coefficients
        old_prec = self.precision()
        new_coeff_list = [self[n] * (n^(-a)) for n in range(1, old_prec + 1)]

        ## Return the new Dirichlet series
        R = self.base_ring()
        return dirichlet_series(new_coeff_list, is_exact_flag=self.is_exact())

dirichlet_series = DirichletSeries

##############################################################################################################
## Some simple constructors:
## =========================

def zeta__series(n):
    """
    Returns the Dirichlet series of the Riemann zeta function, with precision n.
    """
    return dirichlet_series(QQ, "zeta", n)

def L__series(chi, n):
    """
    Returns the L-series of a quadratic Dirichlet character chi, with precision n.
    """
    return dirichlet_series(QQ, chi, n)
