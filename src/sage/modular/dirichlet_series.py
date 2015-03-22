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

import types



###################################################################
## Class for formal manipulation of Dirichlet series over a ring ## 
###################################################################


class DirichletSeries():
    """
    Class for manipulation of Dirichlet series over a ring R.
    """


    def __init__(self, R, B=None, C=None, is_exact_flag=False, default_var_string='s'):
        """

        INPUT:
            DirichletSeries(R, coeff_list) -- initialize a Dirichlet series 
                    with coeffs in R from a list of coefficients

            DirichletSeries(R, f, b) -- Initialize a Dirichlet series over R 
                    from an arithmetic function, with b terms (so + O((b+1)^s))

            DirichletSeries(object) -- create a Dirichlet series from an object,
                    which should call the .dirichlet_series()method of that object.

        We will be able to set the formal variable later if we want... but for 
        now it defaults to s.  We'll need to handle this when we do substitutions...
        
        
        EXAMPLES:
            sage: DirichletSeries(ZZ, "zeta", 20)
            1 + 1/(2^s) + 1/(3^s) + 1/(4^s) + 1/(5^s) + 1/(6^s) + 1/(7^s) + 1/(8^s) + 1/(9^s) + 1/(10^s) + 1/(11^s) + 1/(12^s) + 1/(13^s) + 1/(14^s) + 1/(15^s) + 1/(16^s) + 1/(17^s) + 1/(18^s) + 1/(19^s) + 1/(20^s) + O(21^(-s))
            sage: DirichletSeries(ZZ, "zeta", 10)
            1 + 1/(2^s) + 1/(3^s) + 1/(4^s) + 1/(5^s) + 1/(6^s) + 1/(7^s) + 1/(8^s) + 1/(9^s) + 1/(10^s) + O(11^(-s))
            sage: DirichletSeries(ZZ, "zeta", 2)
            1 + 1/(2^s) + O(3^(-s))
        
            
            
        """
        ## SANITY CHECK:  Check that R is a ring
        ## -------------------------------------
        if not is_Ring(R):
            raise TypeError, "The first argument must be a ring!"
        

        ##  Set some default variables
        tmp_coeff_list = None
        defining_function = None


        ## Parse the construction from an arithmetic function:
        ## ---------------------------------------------------
        if C != None:

            ## Allow functions to be passed
            if isinstance(B, types.FunctionType):
                try:
                    f = B(1)  ## Make sure we can evaluate the function at some natural number
                    f = B  ## Set the function
                except:
                    raise TypeError, "The given arithmetic function f(n) doesn't make sense at n=1."

            ## Allow special names for the arithmetic function (if its a string)
            elif isinstance(B, str):
                if B == "zeta":
                    f = lambda x: 1
                elif (B == "moebius") or (B == "mu"):
                    f = moebius
                elif B == "identity":
                    f = lambda x: kronecker_delta(x,1)

            ## Error if no allowed function initialization is detected
            else:
                raise TypeError, "Something is wrong with the arithmetic function argument..."


            ## Initialize the Internal variables
            self.__base_ring = deepcopy(R)
            self.__creation_function_expression = [f]
            self.__coefficient_list = ['X'] + [f(n)  for n in range(1, C+1)]
            self.__number_of_coeffs = len(self.__coefficient_list) - 1
            self.__is_exact = False
            self.__s_variable = deepcopy(default_var_string)
            

                    
    
        ## Parse the construction from a list:
        ## -----------------------------------
        elif isinstance(B, list): #and isinstance(C, bool):
            tmp_coeff_list = [R(x)  for x in B]
    
            ## Initialize the Internal variables
            self.__base_ring = deepcopy(R)
            self.__creation_function_expression = []
            self.__coefficient_list = ['X'] + tmp_coeff_list
            self.__number_of_coeffs = len(self.__coefficient_list) - 1
            self.__is_exact = is_exact_flag
            self.__s_variable = deepcopy(default_var_string)


        
        
    def base_ring(self):
        """
        Returns the base ring of the Dirichlet series.
        """
        return deepcopy(self.__base_ring)

        
    def precision(self):
        """
        Returns the computed precision of the Dirichlet series.
        """
        return self.__number_of_coeffs
    

    def extend_precision_to(self, n):
        """
        Extend the precision of the Dirichlet series to be O((n+1)^s).
        """
        pass


    def is_Eulerian(self):
        """
        Verify that the series is Eulerian up to the given precision, and cache the result.
        """
        pass


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
            return deepcopy(self.__coefficient_list[n])
        elif self.is_exact() and (n > self.precision()):    ## Allow arbitrary evaluation precision of exact series! =)
            return self.base_ring()(0)
        else:
            raise TypeError, "The coeffecient you requested (n = " + str(n) + ") is out of the precomputed range."


#    def __call__(self, n):



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

        EXAMPLES:
            sage: D1 = DirichletSeries(ZZ, [1,1,1,1]);
            sage: 
            sage: D1 + D1
            2 + 2/(2^s) + 2/(3^s) + 2/(4^s) + O(5^(-s))
            sage: D1 + D1 + D1
            3 + 3/(2^s) + 3/(3^s) + 3/(4^s) + O(5^(-s))

            sage: D2 = DirichletSeries(ZZ, "zeta", 10); D2
            1 + 1/(2^s) + 1/(3^s) + 1/(4^s) + 1/(5^s) + 1/(6^s) + 1/(7^s) + 1/(8^s) + 1/(9^s) + 1/(10^s) + O(11^(-s))
            sage: D1 + D2
            2 + 2/(2^s) + 2/(3^s) + 2/(4^s) + O(5^(-s))

        """

        ## SANITY CHECK -- the rings are the same?
        if self.base_ring() != other.base_ring():
            raise NotImplementedError, "For now the base rings must be the same!  TO DO: Add ring coersion!"


        ## Assume that these are just given by lists of coefficients for now!
        sum_coeff_list = [self[i] + other[i]  for i in range(1, min(self.precision(), other.precision()) + 1)]
        D_sum = DirichletSeries(self.base_ring(), sum_coeff_list)

        return D_sum

        #####################
        ## CONVERSION RULES:
        ## -----------------
        ## Preserve the form of the local variable if it's the same in both, otherwise use the self variable.
        
        ## Take the min of the precisions if they're both made from lists (no defining functions), or bath from defining functions.
        
        ## Take the biggest precision if there is one list and one defining function.





    def __sub__(self, other):
        """
        Form the sum of two Dirichlet series.

        EXAMPLES:
            sage: D1 = DirichletSeries(ZZ, [1,1,1,1]);
            sage: D1 - D1

            sage: D2 = DirichletSeries(ZZ, "zeta", 10); D2
            1 + 1/(2^s) + 1/(3^s) + 1/(4^s) + 1/(5^s) + 1/(6^s) + 1/(7^s) + 1/(8^s) + 1/(9^s) + 1/(10^s) + O(11^(-s))
            sage: D1 - D2

        """
        return self + (other * (-1))



    def __mul__(self, other):
        """
        Define the product of two Dirichlet series, or of a Dirichlet series and a number.

        EXAMPLES:
            sage: D1 = DirichletSeries(ZZ, [1,1,1,1]);
            sage: D1 = DirichletSeries(ZZ, [1,1,1,1]); D1
            1 + 1/(2^s) + 1/(3^s) + 1/(4^s) + O(5^(-s))
            sage: D1 * 3
            3 + 3/(2^s) + 3/(3^s) + 3/(4^s) + O(5^(-s))
            sage: 3 * D1 
            ---------------------------------------------------------------------------
            TypeError                                 Traceback (most recent call last)

            /Users/jonhanke/Dropbox/SAGE/sage-4.6/<ipython console> in <module>()

            /Users/jonhanke/Dropbox/SAGE/sage-4.6/local/lib/python2.6/site-packages/sage/structure/element.so in sage.structure.element.RingElement.__mul__ (sage/structure/element.c:11399)()

            /Users/jonhanke/Dropbox/SAGE/sage-4.6/local/lib/python2.6/site-packages/sage/structure/coerce.so in sage.structure.coerce.CoercionModel_cache_maps.bin_op (sage/structure/coerce.c:6995)()

            TypeError: unsupported operand parent(s) for '*': 'Integer Ring' and '<type 'instance'>'
            sage:


            sage: D1 = DirichletSeries(ZZ, [1,1,1,1]); D1
            1 + 1/(2^s) + 1/(3^s) + 1/(4^s) + O(5^(-s))
            sage: D1 * D1
            1 + 2/(2^s) + 2/(3^s) + 3/(4^s) + O(5^(-s))
            sage: 
            sage: 
            sage: D2 = DirichletSeries(ZZ, "zeta", 8); D2
            1 + 1/(2^s) + 1/(3^s) + 1/(4^s) + 1/(5^s) + 1/(6^s) + 1/(7^s) + 1/(8^s) + O(9^(-s))
            sage: 
            sage: 
            sage: 
            sage: D22 = D2 * D2
            new_coeff_list = [1, 2, 2, 3, 2, 4, 2, 4]
            sage: D22
            1 + 2/(2^s) + 2/(3^s) + 3/(4^s) + 2/(5^s) + 4/(6^s) + 2/(7^s) + 4/(8^s) + O(9^(-s))


        """
        R = self.base_ring()
        try:
            ## Scalar multiplication
            scale_factor = R(other)
            return DirichletSeries(R, [scale_factor * self[i]  for i in range(1, self.precision() + 1)])

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
            return DirichletSeries(R, new_coeff_list, is_exact_flag=new_exact_flag)


    def __pow__(self, n):
        """
        Take any integer power of the current Dirichlet series.

        EXAMPLES:
            sage: D2 = DirichletSeries(ZZ, "zeta", 8); D2
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
        Identity_series = DirichletSeries(R, "identity", self.precision(), is_exact_flag=True)
        
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

        ## TO DO: Use binary to break the powers down to minimize the number of operations!




    def inverse(self):
        """
        Compute the inverse Dirichlet series under Dirichlet multiplication.

        sage: D2 = DirichletSeries(ZZ, "zeta", 8); D2
        1 + 1/(2^s) + 1/(3^s) + 1/(4^s) + 1/(5^s) + 1/(6^s) + 1/(7^s) + 1/(8^s) + O(9^(-s))
        sage: 
        sage: 
        sage: D2.inverse()
        1 + -1/(2^s) + -1/(3^s) + 0/(4^s) + -1/(5^s) + 1/(6^s) + -1/(7^s) + 0/(8^s) + O(9^(-s))
        sage: D2 * D2.inverse()
        1 + 0/(2^s) + 0/(3^s) + 0/(4^s) + 0/(5^s) + 0/(6^s) + 0/(7^s) + 0/(8^s) + O(9^(-s))
        """
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

        return DirichletSeries(R, extended_inv_coeff_list[1:])



#    def is_approx_identity(self):


    def list_of_coefficients(self):
        """
        Returns the list of coefficients of the Dirichlet series, where the zeroth coefficient is 'X'.
        """
        return deepcopy(self.__coefficient_list)
        


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
        return DirichletSeries(R, new_coeff_list[1:], is_exact_flag=self.is_exact())

        
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
        return DirichletSeries(R, new_coeff_list, is_exact_flag=self.is_exact())



##############################################################################################################
## Some simple constructors:
## =========================

def zeta__series(n):
    """
    Returns the Dirichlet series of the Riemann zeta function, with precision n.
    """
    return DirichletSeries(QQ, "zeta", n)


def L__series(chi, n):
    """
    Returns the L-series of a quadratic Dirichlet character chi, with precision n.
    """
    return DirichletSeries(QQ, chi, n)




