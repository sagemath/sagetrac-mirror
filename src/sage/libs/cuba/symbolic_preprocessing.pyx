"""
Symbolic preprocessor for the cuba integration in cuba.pyx
"""


#*****************************************************************************
#       Copyright (C) 2013 Svetlin Tassev <tassev@astro.princeton.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

########################################################################

from sage.misc.functional import integrate
from sage.symbolic.assumptions import assume
from multiprocessing import Pool,TimeoutError
from time import time
from sage.symbolic.expression import Expression

def _integ(toint):
    """
    Do a 1-d symbolic integration using integrate(). 
    Take a tuple containing ``(integrand,integration variable, 
    list of free variables)`` .The result is either an Expression or 
    a None. If the result contains an integrate() (i.e. Maxima did 
    not succeed) it returns a None, as it is useless for numerics, 
    and this allows for automatically terminating further 
    integration attempts for this particilar variable order by its 
    caller ``_multidim_analytical_integration_unit_cube`` . See that 
    function for more. This is tested by all the examples in the 
    rest of this file which use symbolic integration. Here is one 
    test for the sake of it::
    
        sage: from sage.libs.cuba.symbolic_preprocessing import _integ
        sage: y,z = var('y z')
        sage: _integ((x^2*y*sin(z),x,[y,z]))
        1/3*y*sin(z)
    
    
    """
    for v in toint[2]: # The analytical integration is done over the unit hypercube, 
                       # so make that explicit to Sage:
        assume(v > 0)
        assume(v < 1)

    l = len(toint[0].variables())
    
    
    r = integrate(toint[0],toint[1],0,1)
    if not isinstance(r,Expression):
        return r
    
    if l > len(r.variables()): # This check ensures that r does not contain 
                               # integrate() expressions. 
                               # Those are useless for numerics, and moreover, 
                               # if passed through the 
                               # pools, result sporadically in segfaults!!!
            return r
    else:
        return None



########################################################################
    
    

def  _multidim_analytical_integration_unit_cube(func, timelimit, dimlimit=0):
    """ 
    Return the analytical multidimensional integral of a 
    function over the unit hypercube. 
    
    May return only a partially integrated function, if the integral 
    cannot be integrated in all integration variables, either in 
    principle (by Maxima), or within the ``timelimit`` and 
    ``dimlimit`` .
    
    This function should be considered as a preprocessor for 
    ``multidim_integration``.
    
    INPUT:
    
    -   ``func`` -- the integrand, which must be an Expression. The 
        integral is assumed to be over all ``func.variables()`` of 
        ``func``, and therefore the result ultimately should be numeric. 
        However, the analytical integral is not always possible to be 
        reduced to a number. In such cases, the result is a function 
        which is to be integrated over the unit hypercube. In many 
        cases, the resulting function depends on fewer variables than 
        ``func``, and therefore the dimensionality of the integral is 
        reduced.
    
    -   ``timelimit`` -- a number. It gives the number of seconds 
        integration should be attempted before a result is returned. 
   
    -   ``dimlimit``  -- an integer (default:0). It gives the 
        minimum number of integration variables to which the integrator 
        should try to reduce the integrand. See ``multidim_integration``.
   
    OUTPUT: 
    
    A tuple containing:
    
    ``(ndim,result)``
    
    where ``result`` is either a number (for which ``ndim`` =0), in 
    which case the integral has been fully performed; or an 
    Expression, which is to be integrated over the unit hypercube in 
    all ``result.variables()`` . The latter result is obtained when 
    only some of the integrals are possible to do analytically by 
    Sage.
    
    
    .. WARNING::
    
        -   The second warning in the Warning section of 
            ``multidim_integration`` applies to this function as well.
        -   This function should mostly be used as a low-level 
            version of ``multidim_integration``. Use that function 
            instead as most debugging checks are done there. Also, note 
            that the examples in ``multidim_integration`` test this 
            function as well.
   
    ALGORITHM:
      
    Take a function and try to integrate it symbolically on the 
    unit hypercube within timelimit by running many forks in 
    parallel. Each fork has a different integration variable order. 
    We do that because whether an integration succeeds or not and how fast 
    depends on that order. Compare:
    
    integrate(integrate(2000*y^2*sin(2000*x*y),(y,0,1)),(x,0,1)) 
    
    and
    
    integrate(integrate(2000*y^2*sin(2000*x*y),(x,0,1)),(y,0,1)) 
    
    
    EXAMPLES::
        
        sage: from sage.libs.cuba.symbolic_preprocessing import _multidim_analytical_integration_unit_cube
        sage: y,z,w=var('y z w')
        sage: f=2000*y^2*sin(2000*x*y)
        sage: # May fail for very slow computers within set timelimit:
        sage: _multidim_analytical_integration_unit_cube(f,10) # long time (1sec). 
        (0, -1/4000000*cos(2000) - 1/2000*sin(2000) + 2000001/4000000)
        sage: # For the next example, one should get one of the answers below depending on
        sage: # which proceeds faster -- the integral in x or in y. The result below is 
        sage: # twice as fast as the one above as we do not reduce the integral to zero 
        sage: # integration variables. 
        sage: _multidim_analytical_integration_unit_cube(f,10,dimlimit=1) # random output long time
        (1, -y^2*(cos(2000*y)/y - 1/y))
        (1, -1/2000000*((2000000*x^2 - 1)*cos(2000*x) - 2000*x*sin(2000*x))/x^3 - 1/2000000/x^3)
        sage: f=20000*y^2*sin(20000*x*y*z)
        sage: # Output below may have x instead of z, depending on which fork 
        sage: # finishes first.
        sage: _multidim_analytical_integration_unit_cube(f,10) # random output long time
        (1, 1/400000000*(200000000*z^2 - 20000*z*sin(20000*z) - cos(20000*z))/z^3 + 1/400000000/z^3)
        sage: # The calculation below is 3x faster than the above example, as the  
        sage: # integration terminates immediately when the integrand is reduced to 
        sage: # one integration variable. 
        sage: _multidim_analytical_integration_unit_cube(f,10,dimlimit=1) # random output, long time
        (1, 1/400000000*(200000000*x^2 - 20000*x*sin(20000*x) - cos(20000*x))/x^3 + 1/400000000/x^3)
        sage: f=20000*y^2*sin(20000*x^x*y^y*z^z)
        sage: # The result below cannot be integrated in any of its variables:
        sage: _multidim_analytical_integration_unit_cube(f,10) # long time
        (3, 20000*y^2*sin(20000*x^x*y^y*z^z))
        sage: f=20000*y^2*sin(20000*x*y^y*z^z) 
        sage: # If we set dimlimit=2 below, the computation is 5 times faster as further
        sage: # reduction of the integral (below that dimlimit) cannot be done.
        sage: # May fail for very slow computers within the timelimit:
        sage: _multidim_analytical_integration_unit_cube(f,10) # long time 
        (2, -(y^(-y)*z^(-z)*cos(20000*y^y*z^z) - y^(-y)*z^(-z))*y^2)
        sage: f=2000*y^2*sin(2000*x*y*sin(z))
        sage: _multidim_analytical_integration_unit_cube(f,10)  # long time . random output
        (1, -1/2000000*(2000000*(2*cos(2*z) - 1)*cos(4*z) - 1000000*cos(4*z)^2 - 4000000*cos(2*z)^2 + 1000*(cos(4*z) - 2*cos(2*z) + 1)*cos(3*z + 2000*sin(z)) - (cos(4*z) - 2*cos(2*z) + 1)*cos(2*z + 2000*sin(z)) - 1000*(cos(4*z) - 2*cos(2*z) + 1)*cos(z + 2000*sin(z)) + 1000*(cos(4*z) - 2*cos(2*z) + 1)*cos(-z + 2000*sin(z)) - (cos(4*z) - 2*cos(2*z) + 1)*cos(-2*z + 2000*sin(z)) - 1000*(cos(4*z) - 2*cos(2*z) + 1)*cos(-3*z + 2000*sin(z)) - 1000000*sin(4*z)^2 + 4000000*sin(4*z)*sin(2*z) - 4000000*sin(2*z)^2 + 1000*(sin(4*z) - 2*sin(2*z))*sin(3*z + 2000*sin(z)) - (sin(4*z) - 2*sin(2*z))*sin(2*z + 2000*sin(z)) - 1000*(sin(4*z) - 2*sin(2*z))*sin(z + 2000*sin(z)) - 1000*(sin(4*z) - 2*sin(2*z))*sin(-z + 2000*sin(z)) + (sin(4*z) - 2*sin(2*z))*sin(-2*z + 2000*sin(z)) + 1000*(sin(4*z) - 2*sin(2*z))*sin(-3*z + 2000*sin(z)) + 4000000*cos(2*z) - 1000000)/(cos(4*z)^2*sin(z) + 4*cos(2*z)^2*sin(z) + sin(4*z)^2*sin(z) - 4*sin(4*z)*sin(2*z)*sin(z) + 4*sin(2*z)^2*sin(z) - 2*(2*cos(2*z)*sin(z) - sin(z))*cos(4*z) - 4*cos(2*z)*sin(z) + sin(z)) - 1/1000000*(cos(4*z)*cos(2*z) - 2*cos(2*z)^2 + sin(4*z)*sin(2*z) - 2*sin(2*z)^2 + cos(2*z))/(cos(4*z)^2*sin(z) + 4*cos(2*z)^2*sin(z) + sin(4*z)^2*sin(z) - 4*sin(4*z)*sin(2*z)*sin(z) + 4*sin(2*z)^2*sin(z) - 2*(2*cos(2*z)*sin(z) - sin(z))*cos(4*z) - 4*cos(2*z)*sin(z) + sin(z)))


    .. NOTE::
    
        -   When modifying the code, make sure that the function can 
            return numerical answers (if any) as they become available, 
            not waiting for all forks to finish. See, for example, the 
            first example above. Also, make sure all forks are killed in 
            that case before the function returns. And test extensively 
            to make sure that no segfaults appear randomly especially by 
            issuing one and the same command in rapid succession.
        
        -   When modifying the code, make sure that your 
            modifications do not result in speed regressions.
    
    
    
    AUTHORS:
    
    - Svetlin Tassev (2013-12)
    
        
    """
    
    
    
    # _name__ == '__main__':

     
    result=func
    vars = list(func.variables())
    n_vars = len(vars)
    
    it = []
    pool = []
    ns = []
    v_arr = []

    start = time()

    try:
        pool.append( Pool(processes=n_vars) )
    
        it.append( pool[0].imap_unordered(_integ, zip([func]*n_vars,vars,[vars]*n_vars)) )
        ns.append( n_vars )
        v_arr.append( n_vars )

        nthreads = 1
        next_thread = 0
        nmax = n_vars
        
        while (nthreads > 0) and ((time() - start) < timelimit):
            if ns[next_thread] > 0:
                try:
                    func = it[next_thread].next(timeout = 0.01)
                    
                    
                    ns[next_thread] -= 1
                    if func != None:
                        
                        if not isinstance(func, Expression):
                            result = func.n()
                            nthreads = 0
                            n_vars = 0
                            break
                            
                        vars = list(func.variables())
                        n_vars1 = len(vars)



                        if n_vars1 <= dimlimit:
                            result = func
                            n_vars = n_vars1
                            nthreads = 0
                            break
                        
                        nthreads += 1
                        if nmax > n_vars1:
                            nmax = n_vars1
                            n_vars = n_vars1
                            result = func
                           
                        
                        v_arr.append( n_vars1 )
                        pool.append( Pool(processes=n_vars1) )
                        it.append( pool[nthreads-1].imap_unordered(_integ, zip([func]*n_vars1,vars,[vars]*n_vars1)) )
                        ns.append( n_vars1 )
                        
                    next_thread = (next_thread+1) % nthreads
                except (TimeoutError):
                    next_thread = (next_thread+1) % nthreads
            else:
                pool.pop(next_thread)
                it.pop(next_thread)
                ns.pop(next_thread)
                v_arr.pop(next_thread)
                nthreads -= 1 
                next_thread = next_thread % nthreads
    finally:
        for i in pool:
            if not isinstance(i,int):
                 i.terminate()
                 i.join()
    return n_vars, result

