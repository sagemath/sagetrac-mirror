"""
Multidimensional Numerical Integration

Return the multidimensional numerical integral of a function.

    EXAMPLES::

    Let us compute the integral of `x^2\sin(y)` for `x` in the interval
    `(0,3)` and `y` in the interval `(\pi,0.1+3x)` ::

        sage: from sage.libs.cuba.cuba import *
        sage: y = var('y')
        sage: r = multidim_integration(x^2*sin(y),(x,0,3),(y,pi,0.1+3*x),verbose=0) 
        sage: print "result = "+str(r[0])+"+-"+str(r[1]) # abs tol 3e-2
        result = -9.30923867705+-3.85710195604e-05

    See ``multidim_integration`` and 
    ``multidim_integration_unit_cube`` for more examples.
    
    AUTHORS:

    - Svetlin Tassev (2013-12): initial version

"""

#*****************************************************************************
#       Copyright (C) 2013 Svetlin Tassev <tassev@astro.princeton.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************



include 'sage/ext/cdefs.pxi'
include 'sage/ext/interrupt.pxi'
include 'cuba.pxi'
from types import BuiltinFunctionType,FunctionType
from sage.symbolic.expression import Expression
from sage.ext.fast_eval cimport FastDoubleFunc
from sage.ext.fast_callable import fast_callable
from sage.calculus.functional import diff
from sage.matrix.constructor import matrix
from sage.misc.functional import det
from sage.symbolic.ring import SR
from sage.libs.cuba.symbolic_preprocessing import _multidim_analytical_integration_unit_cube
from sage.libs.cuba.cython_compilation import _compile_with_cython_and_then_integrate
from libc.stdlib cimport malloc, free
from sage.gsl.integration import numerical_integral
import inspect


# This function will be evaluated thousands, possibly even millions of times. 
# To optimize it, we split it into three copies:
cdef int _compiled_functional(const int *ndim, const double x[], 
                              const int *ncomp, double f[], 
                              void *userdata):
    """
    A wrapper for Cuba. See ``multidim_integration_unit_cube``.
    """
    t=[]
    for i in range(ndim[0]):
        t.append(x[i])
    value = (<FastDoubleFunc>userdata)(*(t))
    f[0] = value
    return 0

cdef int _compiled_fast_functional_vec(const int *ndim, const double x[], 
         const int *ncomp, double f[], void *userdata):
    """
    A wrapper for Cuba. See ``multidim_integration_unit_cube``.
    """
    value = (<FastDoubleFunc>userdata)._call_c(x) 
    #value = (<FastDoubleFunc>userdata)(x) 
    f[0] = value
    return 0
    
cdef int _compiled_functional_vec(const int *ndim, const double x[], 
         const int *ncomp, double f[], void *userdata):
    """
    A wrapper for Cuba. See ``multidim_integration_unit_cube``.
    """
    t=[]
    for i in range(ndim[0]):
        t.append(x[i])
    value = (<FastDoubleFunc>userdata)(t) 
    f[0] = value
    return 0
# End of three copies




cdef _jacobian_and_substitution(ranges,new_vars,simplify):
    """
    Take integration intervals and the new variable 
    names that live on the unit hypercube and compute the Jacobian 
    of the transformation between the integration variables and the 
    hypercube variables, as well as the transformation itself 
    between them. Passing a simplify value of True, will try to 
    simplify the results. See ``multidim_integration``.
    """
    

    X = []
    v = []
    tosub = []
    transform = []
    for i in range(len(ranges)):
            x = ranges[i][0]
            v.append(x)
            a = ranges[i][1]# a and b can explicitly depend on the integration 
                            # variables. we eliminate them below.
            b = ranges[i][2]
            try:
                X.append((x-a)/(b-a))
            except(ZeroDivisionError):
                return 0,0  # NOTE: The function returns zero for zero 
                # intervals (i.e. when in fact 1/jac=0). So, in this case 
                # the integrator must return 0. When in fact jac=0, we raise 
                # an exception below.
            tosub.append((x,new_vars[i]*(b-a)+a))
            transform.append(new_vars[i]*(b-a)+a)
    jac = det(matrix([[diff(xx,vv) for xx in X] for vv in v]))
    
    if (jac == 0): #So, when in fact jac=0, we raise an exception
        raise Exception('Error: Zero jacobian! Check your integration range.')
    
    tosub = dict(tosub)
    
    # This will eliminate variables for all integration ranges which make any sense. 
    # An example which does    make sense: \int^1_0 dx \int^1_x dy f(x,y)
    # An example which doesn't make sense: \int^y_0 dx \int^1_x dy f(x,y)
    for i in range(len(ranges)):
        tmp = []
        for t in transform:
            t1 = t.subs(tosub)
            if (simplify):
                t1 = t1.simplify_full()
            tmp.append(t1)
        transform = tmp
        
        jac = jac.subs(tosub)
        if (simplify):
            jac = jac.simplify_full()

            
    
    for i in transform:
        for vv in v:
            if i.has(vv):
                raise  Exception('Error: Cannot eliminate variables! Check your integration range.')
    
    if (jac == 0): #So, when in fact jac=0, we raise an exception
        raise Exception('Error: Zero jacobian! Check your integration range.')
        
    return jac,dict(zip(v,transform))








    
    
    
    
    
    
    
    
    
    
    
    
    
    




def multidim_integration(func, *ranges, compile_with_cython=True,
                         preprocess_string=True, only_inspect=False, 
                         verbose=1, cuba_verbose=0, replace_dictionary=[], 
                         additional_imports="", symbolic=False,
                         symbolic_dimension_limit=1, symbolic_time_limit=10,
                         simplify=False,
                         **kwargs):
    r'''
    Calculate the numerical multidimensional integral of a function along
    with an error estimate and convergence information.
    
    To speed up the execution, the integrator makes an honest effort 
    to compile the function before the actual integration. The 
    integrator offers many options to explore and modify the code, 
    which is to be compiled. As an option, symbolic integration can 
    be attempted before the numerical integration, which can be 
    especially useful when the numerical integration is slow to 
    converge to a requested precision. The result is always 
    converted to a floating point number.
        
        
    INPUT:
    
    -   ``func`` -- the integrand. It can be:
    
        - an Expression
        - a FunctionType, BuiltinFunctionType, FastDoubleFunc, a fast_callable (and 
          probably other similar types)
        - a constant
        - a string.  In that case, the name of the integrand is 
          always assumed to be "the_supplied_integrand", and the 
          integration variables in ``ranges`` must follow the order of 
          the arguments of the function. (see the Example and Warnings 
          sections)
    -   ``ranges`` -- the integration intervals in each variable as 
        tuples, e.g. (x,0,10),(y,x^2,200). Each interval is of the form 
        (variable name, lower bound, upper bound). It can only include 
        expressions or numbers. Infinite intervals (using +/-Infinity) 
        are not supported at this time. The functional interdependence 
        of the intervals is sorted out internally, and so for the 
        example above one can equivalently pass (y,x^2,200),(x,0,10). 
        The integration variable names must match the function argument 
        names. Only for Expressions, the latter can be a subset of the 
        former. For certain functions (e.g. BuiltinFunctionType), for 
        which one cannot extract information about their argument names, 
        the integrator assumes that the order of the integration 
        variables in ``ranges`` matches the order of the function 
        arguments. In such cases, the integrator prints a warning for 
        ``verbose`` >=1.
    -   ``kwargs`` -- arguments to be passed on to 
        ``multidim_integration_unit_cube``. Those include the 
        integration method to be used, integrator switches, etc. See 
        ``multidim_integration_unit_cube``.
    -   ``verbose`` -- integer between 0 and 3 (default:1). Sets the code
        verbosity level. All critical warnings are displayed when 
        ``verbose`` >=1. The calculation proceeds silently for ``verbose`` =0.
    -   ``cuba_verbose`` -- integer (default:0). This is passed on to
        ``multidim_integration_unit_cube`` as a ``verbose`` parameter.
        See ``multidim_integration_unit_cube``.
    -   ``compile_with_cython`` -- boolean (default: True). Has an effect
        if the integrand is not an Expression. Tries to compile the 
        integrand transformed to the unit hypercube variables (needed by 
        Cuba) using ``cython`` if ``_fast_float_`` comilation is 
        not available. Information about 
        the compilation process is available if ``verbose`` >=2. If set 
        to False, the integrator proceeds using the fall-back slow 
        function evaluation. Setting it to False is useful for 
        debugging, or if compilation takes too much time compared to the 
        slow integral evaluation.
    -   ``symbolic`` -- boolean (default: False). When set to True, 
        the code first tries to calculate the integral analytically. If 
        a numerical result is found, then the function returns it. 
        Otherwise, the function tries to reduce the dimensionality of 
        the integral by evaluating as many of the integrals analytically 
        as possible within the ``symbolic_time_limit`` . Available only 
        for ``func`` that is an Expression. 
    -   ``symbolic_time_limit`` -- number (default:10; but for fast 
        machines setting this to more than 1 may not result in 
        improvements) of seconds to try with the analytical evaluation 
        before proceeding with the numerical evaluation of the 
        integrand. Has an effect only if ``symbolic`` =True. Note that 
        this time limit is only approximate.
    -   ``symbolic_dimension_limit`` -- an integer (default:1). The 
        lowest dimension (number of integration variables) to which the 
        symbolic integrator should attempt to reduce the integral before 
        proceeding with the numerical integration. Has an effect only if 
        ``symbolic`` =True. When the dimension is 1, then the default 
        integration method is 'GSL' which is really fast at getting high 
        precision answers. So, the value of 1 is a good compromise 
        between speed and accuracy. For wildly oscillating integrands, 
        however, one may want to try the integral with 
        ``symbolic_dimension_limit`` = 0.
    -   ``preprocess_string`` -- boolean (default: True). Applied 
        only when ``func`` is a string. It tells the integrator whether 
        or not to preprocess ``func``. If set to True, it converts ^ to 
        ** (using the sage ``preparse``), all numbers to floats (so, 
        that fractions such as 1/3 are not rounded to zero), and 
        converts sage functions to equivalent C math.h names, imported 
        with ``libc.math``. To see the result of the preprocessing, one 
        can inspect the final code to be compiled by ``cython`` by 
        setting ``only_inspect`` =True.
    -   ``simplify`` -- boolean (default: False). Rarely useful. If 
        set to True, the code tries to symolically simplify the 
        transformation (and Jacobian of that transformation) from the 
        supplied integration variables to variables on the unit 
        hypercube. That transformation is required by the Cuba library. 
        Slows down the preprocessing, but may result in a speed-up once 
        the function is passed to the integrator.
    -   ``only_inspect`` -- boolean (default: False). Rarely useful. 
        Prints out the code for the integrand transformed to the unit 
        hypercube variables. Does not procede to compile the code or 
        evaluate the integral. Useful for debugging.
    -   ``replace_dictionary`` -- a list of tuples (default:[]) 
        containing two strings, e.g. ``('pi','M_PI')``. Rarely useful. 
        It is used if a cython compilation of the integrand transformed 
        to the unit cube is attempted. Each tuple tells the integrator 
        how to replace typical sage names to names available to cython. 
        Keep in mind that the code passed on to cython for compilation 
        of the integrand imports all ``libc.math`` definitions. 
        ``replace_dictionary`` tuples are added to the default list 
        which includes replacement rules from sage names to available 
        ``libc.math`` names (e.g. ``('arcsin','asin')``). Only useful if 
        a compilation of the integrand fails due to undefined functions, 
        and if the fall-back slow integration is too slow. Usually used 
        in combination with ``additional_imports``.
    -   ``additional_imports`` -- a string (default:''). Rarely 
        useful. Used if a cython compilation of the integrand 
        transformed to the unit cube is attempted. It is inserted in the 
        cython code, and is useful if the integration intervals or the 
        integrand depend on functions not imported by default in cython. 
        Only useful if a compilation of the integrand fails due to 
        undefined functions, and if the fall-back slow integration is 
        too slow. Usually used in combination with ``replace_dictionary``.
    
    OUTPUT:
        
    -   If ``only_inspect`` =False (default), then the function returns 
        a tuple containing: 
        
        ``(integration result, integration error, chi2 probability,``
        ``number of integrand evaluations, fail code, actual method, actual ndim)``
        
        ``integration result`` and ``integration error`` are floats,
        giving the integral result and error estimate.
        
        ``actual method`` gives back the actual integration method 
        used for the final result. Can be one of the following 
        strings: 'ZeroInterval' (for zero integration intervals for 
        which the integrator returns zero), 'Constant' (for constant 
        integrands over the unit hypercube; see the Algorithm 
        section), 'Symbolic' (for results obtained using symblic 
        integration), 'GSL', 'Vegas', 'Cuhre', 'Divonne', 'Suave' 
        (the last five are as generated by 
        ``multidim_integration_unit_cube`` ). 
        
        When ``actual method`` is 'ZeroInterval', 'Constant' or 
        'Symbolic', the code returns zeros for ``integration error`` 
        , ``chi2 probability``, ``number of integrand evaluations`` 
        , ``fail code`` and ``actual ndim`` . Otherwise, those are 
        output by ``multidim_integration_unit_cube``. See  
        ``multidim_integration_unit_cube``. 
    
    -   If ``only_inspect`` =True, the function returns an empty string.
    
    
        
    
    EXAMPLES:
    
    Let us integrate `x\times y+z` over the interval `(0,2)` in all variables::
    
        sage: from sage.libs.cuba.cuba import *
        sage: y,z = var('y z')
        sage: multidim_integration(x*y+z,(x,0,2),(y,0,2),(z,0,2))[0] # rel tol 3e-3
        16.000000000000004
    
    Now let us integrate `xy/z^2` for `x` in the interval `(0,1)`, for `y`
    in ther interval `(x,2x)` and for `z` in the interval `(2y,5.43x^2)`::
    
        sage: from sage.libs.cuba.cuba import *
        sage: y,z = var('y z')
        sage: multidim_integration(x*y/z^2,(x,0,1),(y,x,2*x) # rel tol 3e-3
        ....:                                     ,(z,2*y,5.43*x^2))[0]  
        0.028548568652693552
        
    We can evaluate the above integral by defining the integrand in various other ways::
    
        sage: from sage.libs.cuba.cuba import *
        sage: y,z = var('y z')
        sage: f(x,y,z) = x*y/z^2
        sage: multidim_integration(f,(x,0,1),(y,x,2*x) # rel tol 3e-3
        ....:                               ,(z,2*y,5.43*x^2))[0]  # FAST
        0.028548568652693552
        sage: ff = fast_callable(f,vars=(x,y,z),domain=float) 
        sage: multidim_integration(ff,(x,0,1),(y,x,2*x) # long time rel tol 3e-3
        ....:                               ,(z,2*y,5.43*x^2),verbose=0)[0] # SLOW
        0.028548568652693552
        sage: ff = f._fast_float_(x,y,z) 
        sage: multidim_integration(ff,(x,0,1),(y,x,2*x) # long time rel tol 3e-3
        ....:                               ,(z,2*y,5.43*x^2),verbose=0)[0] # SLOW
        0.028548568652693552
        sage: y,z = var('y z')
        sage: def f(x,y,z):
        ....:       return x*y/z^2
        sage: multidim_integration(f,(x,0,1),(y,x,2*x),(z,2*y,5.43*x^2)  # rel tol 3e-3
        ....:                               ,compile_with_cython=False)[0] # INTERMEDIATE
        0.02854856865269355
        sage: cython("""       
        ....: cpdef f(x,y,z):  
        ....:     return x*y/z**2""") 
        sage: multidim_integration(f,(x,0,1),(y,x,2*x)   # long time rel tol 3e-3
        ....:                               ,(z,2*y,5.43*x^2),verbose=0)[0]   # SLOW
        0.028548568652693566
        sage: f=(                              
        ....: """                                
        ....: cpdef the_supplied_integrand(x,y,z): 
        ....:     return x*y/z^2   # here we can safely use ^ instead of ** since 
        ....:                      # string input integrands are preparsed by sage       
        ....: """
        ....: )
        sage: multidim_integration(f,(x,0,1),(y,x,2*x) # long time, rel tol 3e-3
        ....:                               ,(z,2*y,5.43*x^2),verbose=0)[0]   # SLOW
        0.028548568652693566
    
    The comments above indicate which function definitions result in 
    faster and which in slower execution times. Not all of them 
    might seem reasonable, so let us explain those results in more 
    detail. 
    
    The integrator behind ``multidim_integration`` is 
    ``multidim_integration_unit_cube``, which in turn calls the Cuba 
    library [Cuba]_. Cuba itself integrates a provided integrand 
    over the unit hypercube over all function variables. That means 
    that ``multidim_integration`` has to transform the integration 
    ranges to the unit hypercube first. So, even in cases where the 
    function is already compiled, the integrator has to do some more 
    work to compile the variable transformation (the default behavior
    set using ``compile_with_cython`` =True). If the integrand is 
    an expression then the integrand together with the integration 
    variable transformation are compiled using ``fast_callable`` which
    is fast. In all other cases, however, the compilation (if 
    requested, which is default) proceeds using cython, which 
    requires time -- in this case some substantial fraction of a 
    second. For this particular case compilation is the biggest 
    bottle neck, as the integral requires very few integrand 
    evaluations to converge. So, switching off cython copmilation 
    (``compile_with_cython`` =False) makes this particular integral 
    to be evaluated about 3 times faster. 
    
    However, the above integral was trivial. For integrals requiring 
    tens of thousands or even millions (or more!) of evaluations, 
    compilation is crucial to get a result in any reasonable time.
    So, compare::
    
        sage: from sage.libs.cuba.cuba import *
        sage: y,z = var('y z')
        sage: multidim_integration(x*sin(y)              # rel tol 1e-3
        ....:       ,(y,0,200*x),(x,0,1),maxeval=3e5)[0] # fraction of a sec
        0.5043733652048799
        sage: def f(x,y):
        ....:  return x*sin(y)
        sage: multidim_integration(f                     # rel tol 1e-3 long time
        ....:       ,(y,0,200*x),(x,0,1),maxeval=3e5)[0] # ~2 sec
        0.50435733652048821
        sage: # This example takes twice longer:
        sage: multidim_integration(f                     # rel tol 1e-3 long time
        ....:       ,(y,0,200*x),(x,0,1),maxeval=3e5,compile_with_cython=False)[0] 
        0.5043733652048821
        
    Note that the last two examples above can be sped up 
    substantially if the integrand was compiled beforehand using 
    cython or was provided as a string (see the examples further up).
    
    For integrands which are even slower to converge, such as the 
    wildly oscillating integrand below, one should consider using 
    the ``symbolic`` option::
        
        sage: from sage.libs.cuba.cuba import *
        sage: y,z = var('y z')
        sage: multidim_integration(x*sin(y),(y,0,2000*x) # abs tol 1e-8 long time ~1s
        ....:         ,(x,0,1),symbolic=True,symbolic_dimension_limit=0) 
        (0.499535322112679, 0.0, 0.0, 0, 0, 'Symbolic', 0)
        sage: # Integration variable order is irrelevant for speed
        sage: # as all permutations are tried in parallel:
        sage: multidim_integration(x*sin(y),(x,0,1)         # abs tol 1e-8 long time ~1s
        ....:           ,(y,0,2000*x),symbolic=True,symbolic_dimension_limit=0) 
        (0.499535322112679, 0.0, 0.0, 0, 0, 'Symbolic', 0)
        sage: # With the default symbolic_dimension_limit=1, 
        sage: # the integration proceeds twice faster:
        sage: multidim_integration(x*sin(y),(y,0,2000*x),   # abs tol 1e-8 long time ~1s
        ....:               (x,0,1),symbolic=True,symbolic_dimension_limit=1)
        (0.49953532211267804, 5.603590508274081e-15, 0, 50000, 0, 'GSL', 1)
    
    If one of the integration intervals is a point, then the result is always a zero::
    
        sage: multidim_integration(sqrt(y),(x,0,1),(y,-2*x,-2*x)) # abs tol 1e-8
        (0.0, 0.0, 0.0, 0, 0, 'ZeroInterval', 0)
    
    Let us find the value of `\pi` and compare the performance of the different methods::
        
        sage: from sage.libs.cuba.cuba import *
        sage: y,z = var('y z')
        sage: def g(x,y): 
        ....:     if (x^2+y^2<1.0):
        ....:         return 1.0
        ....:     else:
        ....:         return 0.0
        sage: multidim_integration(g,(x,-1,1) # rel tol 3e-3 long time ~1s
        ....:                       ,(y,-1,1),maxeval=100000,eps_rel=1e-3)[0]  
        3.1417868104610505
        sage: multidim_integration(g,(x,-1,1) # rel tol 3e-3 long time ;not compiled, so 5s
        ....:                       ,(y,-1,1),maxeval=100000
        ....:                       ,eps_rel=1e-3,compile_with_cython=False)[0]   
        3.1417868104610505
        sage: cython("""
        ....: cpdef g(x,y): 
        ....:     if (x**2+y**2 < 1.0):
        ....:         return 1.0
        ....:     else:
        ....:         return 0.0""")
        sage: multidim_integration(g,(x,-1,1) # long time (~1s) rel tol 3e-6 
        ....:                       ,(y,-1,1),verbose=0,eps_rel=1e-6
        ....:                       ,maxeval=1e7,method='Divonne')[0]
        3.1415969838455995
        sage: # Note that this is again slow but for 1000x smaller error! 
        sage: # Could not manage to make Mathematica converge to this precision 
        sage: # within any reasonable time without symbolic processing.
        3.1415969838455995
        sage: # Or in this case one can move the integrand to the integration intervals.
        sage: # Here we can also successfully use symbolic=True
        sage: multidim_integration(1,         #     abs tol 1e-8 (takes 0.01 sec!)
        ....:                       (x,-sqrt(1-y^2),sqrt(1-y^2)),(y,-1,1),eps_rel=1e-12)
        (3.141592653589801, 1.7424572844538184e-12, 0, 50000, 0, 'GSL', 1)
    
    The last calculation is blazingly fast for an insane precision, 
    especially when compared to the preceeding calculations. Let us 
    look under the hood and check why that is::
        
        sage: from sage.libs.cuba.cuba import *
        sage: y = var('y')
        sage: multidim_integration(1,(x,-sqrt(1-y^2),sqrt(1-y^2))
        ....:                       ,(y,-1,1),only_inspect=True)
        The function to be integrated over the unit hypercube is:
        -------------------
        4.00000000000000*sqrt(-(2*X[1] - 1)^2 + 1)
        -------------------
        Only code inspection requested.
        ''
        sage: def g(x,y): 
        ....:     if (x^2+y^2<1.0):
        ....:         return 1.0
        ....:     else:
        ....:         return 0.0
        sage: # g below is renamed to the_supplied_integrand, 
        sage: # and it is never in fact called when inspection is requested.
        sage: # Therefore, we could have defined it as a dummy.
        sage: multidim_integration(g,(x,-1,1),(y,-1,1),only_inspect=True)
        The source code for the compiled function to be integrated over
        the unit hypercube is as follows:
        -------------------------
        <BLANKLINE>
        cargs -O3
        from libc.math cimport *   
        <BLANKLINE>
        from sage.libs.cuba.cuba import multidim_integration_unit_cube
        from sage.libs.cuba.cython_compilation import the_supplied_integrand
        cpdef double func1(X):
            jac1=(4.0)
            vals1=[float(2)*float(X[0]) - float(1), float(2)*float(X[1]) - float(1)]
            return the_supplied_integrand(*vals1)*jac1
        <BLANKLINE>
        x=multidim_integration_unit_cube(func1,x_is_vector=True,ndim=2,verbose=0,**{})
        -------------------------
        Only source code inspection requested.
        ''
    
    and so, the integrand is in fact 1-dimensional for the fast 
    calculation (the integrator for ``ndim`` =1 defaults to the GSL 
    integrator called using ``numerical_integral``); while the 
    integral over the (p|c)ython function `g` is 2-dimensional (in 
    the two vector components of ``X``). 
    
    
    Let us do a 3-dimensional Gaussian integral in a way which shows 
    how one can pass parameters to a function to be integrated:: 
    
        sage: from sage.libs.cuba.cuba import *
        sage: def g(a,b,c):
        ....:     x,y,z = var('x y z')
        ....:     f(x,y,z) = exp(-x^2/2/a^2-y^2/2/b^2-z^2/2/c^2)/(2*pi)^(3/2)/(a*b*c)
        ....:     return multidim_integration(f,(x,-100,100),(y,-100,100)
        ....:                                  ,(z,-100,100),method='GaussianLike')[0]
        sage: g(2.3,1.6,6.4) # rel tol 3e-3
        1.0000116384056714
        
    So, the Gaussian above is properly normalized. Now let us  cut a 
    hole in the integrand and set it to zero for values which are 
    above `0.2f(0,0,0)`. Below is the fastest version of the code I 
    came up with. It uses a string input function. One can do it 
    relatively fast with a cython function (declared using cpdef) as 
    well::
    
        sage: from sage.libs.cuba.cuba import *
        sage: x,y,z = var('x,y,z')
        sage: fHole = ("""
        ....: a=2.3
        ....: b=1.6
        ....: c=6.4
        ....: cpdef f(x,y,z):         # preprocess_string =True by default,
        ....:                         # so using ^ instead of ** is allowed.
        ....:    return exp(-x^2/2/a^2-y^2/2/b^2-z^2/2/c^2)/(2*pi)^(3/2)/(a*b*c)
        ....: v0=f(0,0,0)
        ....: cpdef the_supplied_integrand(x,y,z):
        ....:    v=f(x,y,z)
        ....:    if (v<v0*0.2):
        ....:        return v
        ....:    else:
        ....:        return 0.0
        ....: """)
        sage: # Mathematica 9 never converged for this calculation (using 
        sage: # Boole[] to make the hole) to the correct result given
        sage: # the requested precision. Moreover, if SymbolicProcessing
        sage: # is not explicitly off there, the result it gets is
        sage: # consistently off by 50%.
        sage: multidim_integration(fHole,(x,-100,100)  # rel tol 3e-4 long time (~1s)
        ....:                           ,(y,-100,100),(z,-100,100),verbose=0,
        ....:       method = 'GaussianLike',eps_rel=1e-4,maxeval=1e7)[0] 
        0.35892236759872764
        
    For smaller precisions, for the example above one can in fact 
    use an expression by defining 
    ``fHole(x,y,z)=f(x,y,z)*exp(-(f(x,y,z)/(v0*0.2))**n)``, with 
    ``n`` being large, say 250 (larger precision requires larger 
    ``n``, which at some point crashes the code). But tricks like 
    this are beyond this code.
    
    
    Here is a more involved example requiring an external function 
    import done using a string function. Note that we set 
    ``preprocess_string=False`` since the function string 
    preprocessor converts all integers in the code to floats, which 
    in turn crashes the legendre calculation as its first argument 
    must be an int. As the preprocessor is off, we have to be 
    careful to write down (1.0/3.0) and not (1/3), which is zero 
    according to cython::
    
        sage: from sage.libs.cuba.cuba import *
        sage: y,z=var('y z')
        sage: func = (
        ....: """   # import the fast legendreP from the gsl library:
        ....: cdef extern from "gsl/gsl_sf_legendre.h":      
        ....:     double  gsl_sf_legendre_Pl(int l, double x)
        ....: cpdef the_supplied_integrand(x,y,z):    # we set ``preprocess_string`` =False
        ....:     return 1.0/3.0*gsl_sf_legendre_Pl(4,cos(x*y*z))  # so, write 1.0/3.0
        ....: """
        ....: )
        sage: multidim_integration(func,(x,0,1),(y,x,2*x) # long time (1s) rel tol 5e-6
        ....:                          ,(z,2*y,5*y),eps_rel=1e-6,maxeval=1e6,verbose=0
        ....:                          ,preprocess_string=False)[0] 
        0.07031205742221526
    
    For the given precision, the calculation above (which is slow 
    mostly due to compilation time) is about 40 times faster than in 
    Mathematica where it requires tweaking with the integration 
    parameters to converge.
    
    Let us do an integral which Sage can integrate symbolically in 
    some variables but not all. Let us request a verbose output to 
    observe what is happening::
    
        sage: multidim_integration(x*sin(y*z)     # long time ~2s, abs tol 3e-3
        ....:                       ,(x,0,1),(y,0,2000*x),(z,0,10)
        ....:                       ,symbolic=True,verbose=2,eps_rel=1e-5)  
        Trying symbolic integration.
        Symbolic calculation reduced dimensionality of integral to: 1
        Attempting to compile function using fast_callable ...
        Success!
        (4.990345235599539, 0.00013703927547439887, 0, 50000, 0, 'GSL', 1)
    
    If we did not use ``symbolic`` =True, the above calculation 
    would not converge within any reasonable time. Try the above 
    example with ``verbose`` =3 to see the actual final function, 
    which is integrated.
    
    .. WARNING:: 
        
        -   Pay attention to the warning in the description of 
            ``ranges`` about integration variable order. Let us repeat 
            it here for completeness. For certain functions (e.g. 
            BuiltinFunctionType), for which one cannot extract 
            information about their argument names, the integrator 
            assumes that the order of the integration variables in 
            ``ranges`` matches the order of the function arguments. In 
            such cases, the integrator prints a warning for ``verbose`` 
            >=1. Here is an example where argument placement in 
            ``ranges`` does not matter::
            
                sage: from sage.libs.cuba.cuba import *
                sage: y=var('y')
                sage: multidim_integration(x,(x,0,y),(y,0,1))[0] # rel tol 3e-3
                0.16666666666666677
                sage: multidim_integration(x,(y,0,1),(x,0,y))[0] # rel tol 3e-3
                0.16666666666666677
                sage: def f(x,y): 
                ....:       return x*y-x
                sage: multidim_integration(f,(x,0,y),(y,0,1))[0] # rel tol 3e-3
                -0.04166666666666667
                sage: multidim_integration(f,(y,0,1),(x,0,y))[0] # rel tol 3e-3
                -0.04166666666666667

            But for BuiltinFunctionType it is critical::
            
                sage: from sage.libs.cuba.cuba import *
                sage: y=var('y')
                sage: cython(
                ....: """cpdef f(x,y): 
                ....:    return x*y-x
                ....: """)
                sage: multidim_integration(f,(x,0,y),(y,0,1),verbose=0)[0] # rel tol 3e-3
                -0.04166666666666667
                sage: multidim_integration(f,(y,0,1),(x,0,y))[0] # rel tol 3e-3
                WARNING! Not enough available information for this function type! 
                Therefore, you must make sure that:
                1) the variables in the integration ranges follow the order 
                   of the function argument list.
                2) the number of integration variables and function arguments 
                   is one and the same!
                Example: For a function with arguments (x,a,z), the ranges must 
                         be (x,xmin,xmax),(a,amin,amax),(z,zmin,zmax)
                -0.2083333333333334

            Note that the first calculation gives the correct 
            result. There the order of the function arguments and 
            the integration variables matches. In the second 
            example, the order does not match and the result is 
            wrong. We did not set ``verbose`` =0 in the second 
            example on purpose to show the warning message.

        
        -   The Sage preparser can sometimes lead to unexpected 
            results. To avoid this, when passing a FunctionType 
            `f(x,y,...)` as an argument, pass it as ``f``, and not as 
            ``f(x,y,...)``. Let us see the problem in an example::
            
                sage: from sage.libs.cuba.cuba import *
                sage: y=var('y')
                sage: def f(x,y): 
                ....:       if x<y:
                ....:           return 0
                ....:       else:
                ....:           return 1
                sage: # This will be slow  (~1s) due to non-cython definition of f:
                sage: multidim_integration(f,(x,0,1),(y,0,1))[0]  # rel tol 5e-3 long time
                0.5000942725825704
                sage: multidim_integration(f(x,y),(x,0,1),(y,0,1))[0] # rel tol 1e-8
                1.00000000000000
                sage: f(x,y)
                1
                sage: f(1,0)
                1
                sage: f(0,1)
                0
            
            The exact answer is 0.5, and so the first calculation 
            gets it right, while the second calculation fails. The 
            reason is that when one calls the integrator with 
            `f(x,y)`, then Sage preparses the expression. That 
            returns `f(x,y)=1` even though for numerical arguments 
            `f` returns the expected results (see the output above). 
            The reason for the strange behavior of `f` is that, 
            quite unfortunately, at this time Sage returns False for 
            unknown sybolic comparisons such as `(x<y)` in `f`. So, 
            in the second calculation above, the value of `f(x,y)=1` 
            is passed on to ``multidim_integration``, which in turn 
            produces the correct result but for the unintended 
            function `f(x,y)=1`.
            
            Long story short, one should always make sure that what 
            Sage passes on to the integrator is in fact what one 
            wants to be integrated. One can use one of the many 
            debugging options for that (such as setting 
            ``only_inspect`` =True).
            
        -   One can pass the integrand function as a string. In that 
            case, the name of the integrand is always assumed to be 
            "the_supplied_integrand", and the integration variables in 
            ``ranges`` must follow the order of the arguments of the 
            function. The string definition is then inserted into a 
            cython code as a blob to be compiled without any checks, so 
            in case of errors, one may want to set ``only_inspect`` 
            =True to inspect the code of the function to be integrated. 
            For faster results one should probably define the function 
            using cpdef instead of def. Because the string is directly 
            inserted in a cython source code, one can even import 
            function definitions with the same string (see the examples).
            When ``preprocess_string`` =True (default), the string is 
            preprocessed for convenience to convert all integers to 
            floats, ^ to ** , arcsin to asin, etc. However, this may 
            crash the code as in the example using the Legendre 
            polynomials from GSL above. In such cases one should try 
            with ``preprocess_string`` =False (keeping in mind that now 
            ^ is not automatically converted to ** and that things like 
            1/3 are rounded to zero by cython), and/or set 
            ``only_inspect`` =True to identify the problem.
        
        
    
    TESTS:
    
    Make sure the integral over a constant integrand works::
    
        sage: from sage.libs.cuba.cuba import *
        sage: y,z = var('y z')
        sage: multidim_integration(1,(x,0,2),(y,0,2),(z,0,2)) # abs tol 1e-8
        (8.00000000000000, 0.0, 0.0, 0, 0, 'Constant', 0)
        sage: # This also reduces to a constant after the transformation of 
        sage: # variables to the unit cube. To see that explicitly, run 
        sage: # with only_inspect=True
        sage: multidim_integration(1,(x,0,2),(y,0,2),(z,x,2+x)) # abs tol 1e-8 
        (8.00000000000000, 0.0, 0.0, 0, 0, 'Constant', 0)
    
    Make sure the 1-dim integral also works whether defined through 
    an Expression or a FunctionType::
    
        sage: def g(x):
        ....:     return x
        sage: multidim_integration(g,(x,0,2))[0] # rel tol 3e-3
        2.0
        sage: multidim_integration(x,(x,0,2))[0] # rel tol 3e-3
        2.0
    
    Make sure this works::
    
        sage: multidim_integration(x*sin(x^x*y) # abs tol 3e-8 long time
        ....:                                ,(y,0,x),(x,0,1),symbolic=True) 
        (0.10224975449206862, 1.1352003169936412e-15, 0, 50000, 0, 'GSL', 1)
        
    ... and this ...::
    
        sage: from sage.libs.cuba.cuba import *
        sage: x,y,z = var('x,y,z')
        sage: a=2.3
        sage: b=1.6
        sage: c=6.4
        sage: def f(x,y,z):        
        ....:     return exp(-x^2/2/a^2-y^2/2/b^2-z^2/2/c^2)/(2*pi)^(3/2)/(a*b*c)
        sage: v0=f(0,0,0)
        sage: def fHole(x,y,z):        
        ....:     v=f(x,y,z)  # defining a local variable was leading to issues
        ....:     if (v<v0*0.2):
        ....:         return v
        ....:     else:
        ....:         return 0.0
        sage: multidim_integration(fHole,(x,-100,100)  # rel tol 3e-2 long time (~1s)
        ....:                           ,(y,-100,100),(z,-100,100),verbose=0,
        ....:       method = 'GaussianLike',eps_rel=1e-2,maxeval=1e7)[0] 
        0.3600773268303449
        
    ALGORITHM:
    
    Here is a sketch of the algorithm.
    
    1. The code transforms the variables of integration to new 
    variables over the unit hypercube. It computes the Jacobian of 
    the transformation as well.
    
    2. If the integrand is constant over the unit hypercube, the 
    constant is returned as it is the result of the integration. If 
    not, the evaluation continues on.
    
    3. When symbolic integration is requested for Expression 
    integrands, it is tried on the transformed function, since then 
    it is trivial to handle the assumptions.
    
    4. If symbolic integration is not requested (or is unavailable, 
    or is not successful within the set time limit to reduce the 
    integral to a number), then the code attempts to compile the 
    transformed integrand either using ``fast_callable`` (for 
    Expressions) or using ``cython`` (for functions and strings).
    
    5. The code then passes the integrand on to 
    ``multidim_integration_unit_cube``, which in turn uses either 
    the GSL ``numerical_integral`` (for one dimensional integrals 
    when ``method`` is set appropriately) or the Cuba library for 
    multidimensional numerical integration [Cuba]_. See 
    ``multidim_integration_unit_cube``.
    
    .. TODO::
    
        -   Add support for infinite intervals. 
        -   Make possible to 
            view the integration space partitioning similar to partview, 
            which comes with Cuba.
    
    .. NOTE::
    
        -   When modifying the code, make sure that your 
            modifications do not result in speed regressions.
        
    
    REFERENCES:
    
    ..  [Cuba] T. Hahn. Cuba - a library for multidimensional numerical integration.
        Computer Physics Communications, Volume 168, Issue 2, p. 78-95. (2005).
        Available on http://arxiv.org/abs/hep-ph/0404043
        
    AUTHORS:
    
    - Svetlin Tassev (2013-12)
        
   '''
    
    ranges = list(ranges)
    
    r_vars = [v[0] for v in ranges]
    num_vars = len(r_vars)
    
    
    
    
    # Make sure this is a float and not some symbolic number
    # Needed for when we check for constant integrand; or for passing
    # to the compiled function, when the jacobian is not a constant.
    try:
         func = func.n()
    except (AttributeError,TypeError):
         pass
    
    if isinstance(func, str) and (verbose>=1):
                print "WARNING! Not enough available information for this function type!"
                print "Therefore, you must make sure that:"
                print "1) the variables in the integration ranges follow the order"
                print "   of the function argument list."
                print "2) the number of integration variables and function arguments"
                print "   is one and the same!"
                print "Example: For a function with arguments (x,a,z), the ranges must"
                print "         be (x,xmin,xmax),(a,amin,amax),(z,zmin,zmax)"

    
    if callable(func):
        if isinstance(func, FunctionType):
            #pass_vars=func.func_code.co_varnames # This is a list of strings
            pass_vars = inspect.getargspec(func)[0]
            if (len(pass_vars) != num_vars):
                raise Exception('Number of function arguments and number of integration variables must match for FunctionType!')
            # To compare to pass_vars, we need to convert r_vars to strings:
            str_vars=tuple([str(r) for r in r_vars]) 
            if not(set(str_vars) == set(pass_vars)):
                print "Your integrand depends on the following variables:"
                print pass_vars
                print "While the integration ranges depend on these variables"
                print str_vars
                raise Exception('Wrong integration variables!')
            #rearrange ranges array to match function argument order:
            if (str_vars != pass_vars):
                str_vars = list(str_vars)
                for j in range(len(pass_vars)):
                    for i in range(len(str_vars)):
                        if (str_vars[i]==pass_vars[j]):
                            tmp=ranges[j]
                            ranges[j]=ranges[i]
                            ranges[i]=tmp
                            tmp=str_vars[j]
                            str_vars[j]=str_vars[i]
                            str_vars[i]=tmp
                r_vars=[v[0] for v in ranges]
        elif isinstance(func, FastDoubleFunc) or  hasattr(func,'get_orig_args'):
            if (verbose>=1):
                print "WARNING! I will assume that the variables in the integration ranges"
                print "follow the order of the function argument list."
                print "Example: For a function with arguments (x,a,z), the ranges must"
                print "         be (x,xmin,xmax),(a,amin,amax),(z,zmin,zmax)"
            if hasattr(func,'get_orig_args'):
                num_vars_temp = func.get_orig_args().get('args')
            else:
                num_vars_temp = func.nargs
            if (num_vars_temp != num_vars):
                raise Exception('Number of function arguments and number of integration variables must match for FastDoubleFunc and fast_callable!')
            pass_vars=r_vars
        elif hasattr(func, 'arguments'):
            pass_vars = func.arguments()
        elif hasattr(func, 'variables'):
            pass_vars = func.variables()
        else:
            if (verbose>=1):
                print "WARNING! Not enough available information for this function type!"
                print "Therefore, you must make sure that:"
                print "1) the variables in the integration ranges follow the order"
                print "   of the function argument list."
                print "2) the number of integration variables and function arguments"
                print "   is one and the same!"
                print "Example: For a function with arguments (x,a,z), the ranges must"
                print "         be (x,xmin,xmax),(a,amin,amax),(z,zmin,zmax)"
            pass_vars = r_vars
            
            #try it out:
            #try:
            #    func(*[0 for i in range(num_vars)])
            #except(TypeError):
            #    raise
            #except:
            #    pass
        if not isinstance(func, FunctionType):
            if not (set(pass_vars)<=set(r_vars)):
                print "The integrand depends on the following variables:"
                print pass_vars
                print "While the integration ranges depend on these variables"
                print r_vars
                print "The former should be a subset of the latter!"
                raise Exception('Function and integration variables do not match!')
    
   
    
    new_vars=[]
    for i in range(len(ranges)):                  
        new_vars.append(SR.symbol('new_vars_BEgIn423N'+str(i)+'eN4324QxdD')) 
                                #  new_vars are supposed to be random and,
                                #  therefore, safe, as I'm doing string replace below
    
    # Need to map integration range onto the unit hypercube for Cuba
    jac,tosub=_jacobian_and_substitution(ranges,new_vars,simplify)
    if jac==0: # Note that _jacobian_and_substitution returning zero is a special 
               # case of zero intervals. So, the integral is zero.
        if only_inspect or (verbose>=2):
            print "One or more of the integration intervals is zero. Thus, the integral is zero."
            if only_inspect:
                print "Only code inspection requested."
                return ""
        return 0.0,0.0,0.0,0,0,'ZeroInterval',0
    jac=1.0/jac
    
    
    #handle constant case
    if (not callable(func)) and (not isinstance(func,str)):
        try:
            fjac = (func*jac).n()
            if only_inspect or (verbose>=2):
                print "The integrand is constant over the unit hypercube."
                if only_inspect:
                    print "The integral is "+str(fjac)
                    print "Only code inspection requested."
                    return ""
            return fjac,0.0,0.0,0,0,'Constant',0
        except (AttributeError,TypeError):
            pass
    
    # For expressions, it is straightforward:
    if ((hasattr(func,'subs') and isinstance(func,Expression)) or 
        ((not callable(func)) and (not isinstance(func,str)))):
        if (not callable(func)):
            func=func*jac
        else:
            func=func.subs(tosub)*jac
            
        if (verbose>=3) or (only_inspect):
            print "The function to be integrated over the unit hypercube is:"
            print "-------------------"
            tmp=str(func).replace('new_vars_BEgIn423N','X[').replace('eN4324QxdD',']')
            #tmp=repl("(.*)\|-->","",tmp)
            print tmp
            print "-------------------"

        
        if (not callable(func)):
            if only_inspect or (verbose>=2):
                print "The integrand is constant over the unit hypercube."
                if only_inspect:
                    print "The integral is "+str(func)
                    print "Only code inspection requested."
                    return ""
            return func,0.0,0.0,0,0,'Constant',0 # func was already multiplied by Jacobian above
        
        # If the true ndim is smaller due to the var transformation, it is 
        # taken care of with this since we use an Expression
        # If the code fails here, file a bug report.
        vars=func.variables()
        num_vars=len(vars)
        func=func.function(*vars) # fixing argument list to new vars.
        
        
        
############# symbolic block
        if (symbolic) and num_vars > symbolic_dimension_limit:
            if verbose >= 2:
                print "Trying symbolic integration."
            X = _multidim_analytical_integration_unit_cube(func, symbolic_time_limit, 
                                                    dimlimit=symbolic_dimension_limit)
            if X[0] == 0:
                if (verbose >= 2) or only_inspect:
                    print "Symbolic calculation was successful!"
                    print "The integrand is constant over the unit hypercube."
                    if only_inspect:
                        print "The integral is "+str(X[1])
                        print "Only code inspection requested."
                        return ""
                return X[1].n(), 0.0, 0.0, 0, 0, 'Symbolic', 0
            else:
                vars = list(X[1].variables())
                num_vars=X[0]
                func = fast_callable(X[1],vars=vars,domain=float)  
                if (verbose >= 2) or only_inspect:
                    print "Symbolic calculation reduced dimensionality of integral to: "+str(X[0])
                if (verbose >= 3) or only_inspect:
                    print "The resulting function over the unit hypercube is:"
                    print "----------------------"
                    print str(X[1]).replace('new_vars_BEgIn423N','X[').replace('eN4324QxdD',']')
                    print "----------------------" 
                
########## end symbolic block

        if only_inspect:
            print "Only code inspection requested."
            return ""
                                        
########## compiling Expressions with fast_callable is straightforward and fast:      
        if verbose >= 2:
            print "Attempting to compile function using fast_callable ..."
        if hasattr(func,'_fast_callable_'):
            func = fast_callable(func,vars=vars,domain=float) # fix argument list to new vars
        if verbose >= 2:
            print "Success!"

        return multidim_integration_unit_cube(func,ndim=num_vars,verbose=cuba_verbose,**kwargs)
        # here ndim is set only tentatively the ...unit_cube function will check it. 
        # In reality, it may be a smaller number due to the irrelevance of some of the vars.
        # Although, above this should already be fixed, unless the code is modified ...

########## if not Expression, try cython compilation:


    # For functions it is more work as we first should make an honest effort to compile them,
    # or otherwise it will be too slow:
    
    vals=[v.subs(tosub) for v in r_vars]
    
    
    if (compile_with_cython):
        return _compile_with_cython_and_then_integrate(func, vals,
                       jac,num_vars, additional_imports, replace_dictionary,
                       preprocess_string,
                       verbose, cuba_verbose, only_inspect,**kwargs) 

########## slow integration:

    if (verbose>=2):
        print "Continuing with slow (no cython) evaluation..." 
    
    if isinstance(func,str):
        raise Exception("Slow evaluation does not work for string integrands. Use compile_with_cython=True")

    def func1(X):
        plugin=dict(zip(new_vars,X))
        jac1=jac.subs(plugin)
        if callable(func):
            vals1=[v.subs(plugin) for v in vals]
            return func(*vals1)*jac1
        else:
            return func*jac1
    return multidim_integration_unit_cube(func1, x_is_vector=True, ndim=num_vars, 
                                          verbose=cuba_verbose, **kwargs)
        
        







def multidim_integration_unit_cube(func, x_is_vector=False, ndim=0, 
             eps_rel=1e-3, eps_abs=1e-12, maxeval=50000, verbose=0,
             method='Default', GSL_algorithm = 'qag',GSL_rule=6, seed=1, mineval=0,
             pseudorandomlevel=0, accumulate_points=True, additional_flags=0, 
             retain_statefile=False,statefile='',
             C_key=0,
             V_nstart=1000, V_nincrease=500, 
             V_nbatch=1000, V_gridno=0,
             D_key1=47, D_key2=1, D_key3=1, D_maxpass=5,
             D_border=1.0e-14, D_maxchisq=10, 
             D_mindeviation=0.25,
             D_ldxgiven=0, D_ngiven=0, D_xgiven=[],
             S_nnew=1000,S_flatness=50
    ): 
    ''' 
    Return the numerical multidimensional integral of a 
    function over the unit hypercube along with an error estimate 
    and convergence information.
   
    This function should mostly be used as a low-level version of 
    ``multidim_integration``. If unsure, use that function instead.
   
    INPUT:
    
    -   ``func`` -- the integrand to be integrated in ``ndim`` 
        dimensions. No parameters are allowed as the function must 
        return real numbers for numerical arguments. Argument order is 
        irrelevant as the integration intervals are `(0,1)` for all 
        variables.
    -   ``ndim`` -- integer (default:0). The number of dimensions over which 
        to compute the integral. The integrator tries to figure it out 
        automatically (it also fixes wrong input ``ndim`` values). 
        However, user input may be requested in case that fails.
    -   ``x_is_vector`` -- boolean (default:False). Whether or not 
        the ``ndim`` arguments are separate variables, or are the 
        ``ndim`` components of one vector argument.
    -   ``eps_rel`` -- a number (default:1e-3) setting the relative error tolerance
    -   ``eps_abs`` -- a number (default:1e-12) setting the absolute error tolerance
    -   ``maxeval`` -- an integer (default:50000). Maximum number of function evaluations
    -   ``verbose`` -- an integer between 0 and 3 (default:0). Sets 
        the verbosity level for Cuba (see below).
    -   ``method`` -- a string (default:'Default'). Tells which 
        integration method to use. Available options are:
        
        -   'Default' -- For ``ndim`` =1, it sets the method to 
            'GSL' for ``x_is_vector`` =False and to 'Vegas' for 
            ``x_is_vector`` =True. For higher dimensions it sets the 
            method to 'Cuhre'.
        -   'GSL' -- available for ``ndim`` =1 and ``x_is_vector`` =False. 
            It uses GSL integration through ``numerical_integral``.
        -   Meta methods:
        
            -   'Oscillatory' -- Suitable for oscillatory 
                integrands. Chosen from the available Cuba methods 
                according to  the results on p.23 of [Cuba]_
            -   'ProductPeak' -- Suitable for integrands with several peaks. Chosen as above.
            -   'CornerPeak' --  Suitable for integrands peaking at 
                a corner of the unit cube. Chosen as above.
            -   'GaussianLike' --  Suitable for Gaussian-like integrands. Chosen as above.
            -   'C0continuous' --  Suitable for C0-continuous integrands. Chosen as above.
            -   'Discontinuous' --  Suitable for discontinuous integrands. Chosen as above.
        -   Cuba methods. See [CubaUpdated]_ for detailed descriptions of these methods.
        
            -   'Vegas'
            -   'Cuhre'
            -   'Divonne'
            -   'Suave'
            
    -   ``GSL_algorithm`` and ``GSL_rule`` are directly passed on to 
        ``numerical_integral``. See ``numerical_integral`` for the 
        description of these parameters (without the 'GSL' prefix).
    
    -   ``mineval`` -- an integer (default:0). The minimum number of 
        integrand evaluations for Cuba methods. Corresponds to the Cuba 
        argument of the same name. See [CubaUpdated]_ for the details.
    
    -   ``seed`` -- an integer (default:1). Corresponds to the 
        ``seed`` argument in Cuba. See [CubaUpdated]_ for the details.
    
    -   The Cuba ``flags`` integer argument is set using the 
        arguments of ``multidim_integration_unit_cube`` according to:
        
        ``flags = verbose+4*bool(accumulate_points)+256*pseudorandomlevel + additional_flags``
    
        The names of the ``multidim_integration_unit_cube`` 
        arguments on the rhs above are more or less self-explanatory.
        See [CubaUpdated]_ for the details, however. The defaults 
        are as follows: ``verbose`` =0, ``accumulate_points`` = 
        True, ``additional_flags`` =0.
        
    -   ``retain_statefile`` -- boolean (default:False). Whether to 
        retain a statefile. See [CubaUpdated]_ for the details.
    
    -   ``statefile`` -- string (default:''). The name of the 
        statefile. See [CubaUpdated]_ for the details.
    
    -   The following arguments with their defaults are passed to 
        the corresponding Cuba integrator method. The names (stripped 
        from the underscore prefix) are identical to the ones in 
        [CubaUpdated]_.  See [CubaUpdated]_ for the details.
        
        -   Arguments for Cuhre: ``C_key`` =0
        -   Arguments for Vegas: ``V_nstart`` =1000, ``V_nincrease`` =500, 
            ``V_nbatch`` =1000, ``V_gridno`` =0, ``D_key1`` =47, 
            ``D_key2`` =1, ``D_key3`` =1, ``D_maxpass`` =5, ``D_border`` =1.0e-14, 
            ``D_maxchisq`` =10, ``D_mindeviation`` =0.25
        -   Arguments for Divonne: ``D_ldxgiven`` =0, ``D_ngiven`` =0, 
            ``D_xgiven`` =[]
        -   Arguments for Suave: ``S_nnew`` =1000, ``S_flatness`` =50
    
    OUTPUT:
    
    -   For ``method`` ='GSL', the output is a tuple giving:
        
        
        ``(integration result, integration error, 0 , maxeval , 0, 'GSL', 1 )``
        
        where ``integration result`` and ``integration error`` are 
        floats as computed by ``numerical_integral``, and 
        ``maxeval`` is the maximum number of integrand evaluations 
        requested.
        
    -   For the Cuba methods, the output is a tuple containing:
    
        ``( integration result, integration error, chi2 probability,`` 
        ``number of integrand evaluations, fail code, actual method,``
        ``actual ndim)``
    
        -   ``integration result`` and ``integration error`` are the 
            numerical result and error estimate for the integral 
            obtained using Cuba.
        
        -   ``chi2 probability`` is a float returning the `\chi^2` 
            probability that the returned integration error is not a 
            reliable estimate of the true error. A value <<1 indicates a 
            good error estimate. See the Cuba library documentation 
            [CubaUpdated]_ for more details.
   
        -   ``number of integrand evaluations`` is an integer 
            returning the actual number of integrand evaluations the 
            code made.
   
        -   ``fail code`` is an integer returning the value of the 
            ``fail`` argument of Cuba. Possible values (quoting from the 
            Cuba manual [CubaUpdated]_) are:
            
            -   ``fail`` =0, the desired accuracy was reached
            -   ``fail`` =-1, the dimension of the integral is 
                unavailable for the chosen combination of integrator 
                method and switches
            -   ``fail`` >0, the accuracy goal was not met within 
                the allowed maximum number of integrand evaluations. 
                While 'Vegas', 'Suave', and 'Cuhre' simply return 1, 
                'Divonne' can estimate the number of points by which 
                ``maxeval`` needs to be increased to reach the desired 
                accuracy and returns this value.

        -   ``actual method`` is a string returning the actual Cuba 
            method used for the integration
   
        -   ``actual ndim`` is an integer returning the actual 
            dimensionality of the integral evaluated by Cuba
   
    -   For the meta methods, the result is one of the above.
   
    EXAMPLES:
    
    Let us integrate `x\\times y` on the unit hypercube::
    
        sage: from sage.libs.cuba.cuba import *
        sage: y=var('y')
        sage: multidim_integration_unit_cube(x*y)[0] # rel tol 3e-3
        0.24999999999999994
   
    For more examples see ``multidim_integration`` instead.
   
    TESTS:
   
    Let us try a ``fast_callable`` compiled function as well::
    
        sage: from sage.libs.cuba.cuba import *
        sage: y=var('y')
        sage: f(x,y)=x*y^2
        sage: f=fast_callable(f,vars=(x,y),domain=float)
        sage: multidim_integration_unit_cube(f)[0] # rel tol 3e-3
        0.16666666666666677
   
   
    The integral of a constant is the same constant::
    
        sage: from sage.libs.cuba.cuba import *
        sage: multidim_integration_unit_cube(3) # rel tol 1e-8
        (3.0, 0.0, 0.0, 0, 0, 'Constant', 0)
        
    Let us check that the different methods agree::
        
        sage: from sage.libs.cuba.cuba import *
        sage: def g(x,y):
        ....:     return y*sin(x)
        sage: multidim_integration_unit_cube(g,method='Cuhre')[0] # rel tol 3e-3
        0.2298488470659301
        sage: multidim_integration_unit_cube(g,method='Vegas')[0] # rel tol 3e-3
        0.22969997926596447
        sage: multidim_integration_unit_cube(g,method='Suave')[0] # rel tol 3e-3 long time
        0.22797366049605247
        sage: multidim_integration_unit_cube(g,method='Divonne')[0] # rel tol 3e-3
        0.22986187577090006

    Let us try passing a cython compiled function::

        sage: from sage.libs.cuba.cuba import *
        sage: cython("""
        ....: cpdef g(x,y): 
        ....:     if (x**2+y**2<1.0):
        ....:         return  4.0/3.1415926
        ....:     else:
        ....:         return 0.0""")
        sage: multidim_integration_unit_cube(g,ndim=2)[0] # rel tol 3e-3
        0.9998436417774216
    
    Let us try the 1-dim case::
        
        sage: from sage.libs.cuba.cuba import *
        sage: # This uses the default 'GSL' method:
        sage: multidim_integration_unit_cube(x)[0] # rel tol 3e-3
        0.5
        sage: multidim_integration_unit_cube(x,method='Vegas')[0] # rel tol 3e-3
        0.4999535857167853
        
    .. WARNING::
    
        -   The second warning in the Warning section of 
            ``multidim_integration`` applies to this function as well.
        -   This function should mostly be used as a low-level 
            version of ``multidim_integration``. Use that function 
            instead as most debugging checks are done there. Also, note 
            that the examples in ``multidim_integration`` test this 
            function as well.
   
    ALGORITHM:
    
    Uses the Cuba library for multidimensional numerical integration [Cuba]_.
    
    For 1-dimensional integration, it optionally supports calling 
    the GSL (GNU Scientific Library) integrator implemented using a 
    call to ``numerical_integral``.
    
    
    .. NOTE::
    
        When modifying the code, make sure that your modifications 
        do not result in speed regressions.
    
    REFERENCES:
    
    ..  [CubaUpdated] T. Hahn. Cuba - a library for multidimensional 
        numerical integration. The updated Cuba reference manual 
        distributed with the Cuba source code on 
        http://www.feynarts.de/cuba/
    
    
    AUTHORS:
    
    - Svetlin Tassev (2013-12)
   
   
    '''
    if not callable(func):
        return <double>func,0.0,0.0,0,0,'Constant',0
    
    cdef integrand_t F
    cdef void * userdata
    
    if not x_is_vector: #assume f(x,y,z,...)
            if  hasattr(func, 'arguments'):
                vars = func.arguments()
                ndim = len(vars)
            elif hasattr(func, 'variables'):
                vars = func.variables()
                ndim = len(vars)
            elif hasattr(func, 'func_code'):
                #vars = func.func_code.co_varnames
                vars = inspect.getargspec(func)[0]
                ndim = len(vars)
            elif (ndim>0):
                vars = []
                for i in range(ndim):                  
                    vars.append(SR.symbol('new_vars68bfni8ht37_'+str(i))) 
                                         # new_vars is some random name
                                         
            if hasattr(func, '_fast_callable_'):  
                func = fast_callable(func,vars=vars,domain=float)
                
            if isinstance(func, FastDoubleFunc) :
                ndim = func.nargs
            if hasattr(func,'get_orig_args'):
                ndim = func.get_orig_args().get('args')
            
            if ndim == 0:
                    raise Exception("Cannot infer number of function arguments. Please specify ndim.")
                    
            if (method == 'GSL' or method == 'Default') and ndim == 1:
                if verbose >= 1:
                    print "Will be using the GSL integrator."
                
                r = numerical_integral(func, 0, 1, eps_abs = eps_abs, 
                    eps_rel = eps_rel, rule=GSL_rule, max_points=maxeval,
                    algorithm = GSL_algorithm)
                return r[0],r[1],0,maxeval,0,'GSL',1
                    
            F = _compiled_functional
            userdata = <void *>func
    else:
            if ndim>0:
                if hasattr(func,'_fast_callable_'): 
                    func = fast_callable(func,vars='x',domain=float)
                if isinstance(func, FastDoubleFunc) or  hasattr(func,'get_orig_args'):
                    F = _compiled_fast_functional_vec
                    userdata = <void *>func
                else:
                    F = _compiled_functional_vec
                    userdata = <void *>func
            else:
                raise Exception("Cannot infer dimensionality of x vector. Please specify ndim.")
    
    
    
    
    cdef double prob
    cdef double result
    cdef double err
    cdef int fail
    cdef int flags
    cdef int finalInt
    cdef int neval
    cdef char *statefilePass
    if (retain_statefile):
        statefilePass = statefile
    else:
        statefilePass = NULL 
        
    if (accumulate_points):
        finalInt = 0
    else:
        finalInt = 4
        
    flags = verbose+finalInt+256*pseudorandomlevel+additional_flags
    cdef int ncomp = 1
    
    if method == 'Default':
            method = 'Cuhre' # GSL is chosen higher up for ndim=1 for non-vector arguments
    
    if method == 'Oscillatory':
        if ndim > 1:
            method = 'Cuhre'
        else:
            method = 'Vegas' 
    if method == 'ProductPeak':
        method = 'Vegas' 
    if method == 'CornerPeak':
        if ndim <= 8 and ndim>1:
            method = 'Cuhre'
        else:
            method = 'Vegas' 
    if method == 'GaussianLike':
        if ndim > 1:
            method = 'Divonne'
        else:
            method = 'Vegas'
    if method == 'C0continuous':
        method = 'Vegas'
    if method == 'Discontinuous':
        if (ndim <= 8) and (ndim > 1):
            method = 'Cuhre'
        elif ndim > 8:
            method = 'Vegas' 
        elif ndim == 1:
            method = 'Vegas'
    
    if ndim == 1:
        if method == 'Cuhre' or method == 'Divonne':
                if (verbose >= 1):
                    print "Method unavailable for ndim=1. Switching method to Vegas"
                method='Vegas'
    
    foundMethod=False
    
    if verbose>=1:
        print "Using method = "+method
    
    if method=='Vegas':
            foundMethod=True
            sig_on()
            Vegas(ndim,ncomp,F,userdata,eps_rel,eps_abs, flags, seed, 
                 mineval, maxeval,
                 V_nstart, V_nincrease, V_nbatch, V_gridno, statefilePass, 
                 &neval, &fail, &result,  &err, &prob)
            sig_off()
    
    cdef int nregions
    
    if method=='Suave':
            foundMethod=True
            sig_on()
            Suave(ndim,ncomp,F,userdata, eps_rel, eps_abs, flags, seed, 
                  mineval, maxeval,
                  S_nnew, S_flatness, statefilePass, &nregions, 
                  &neval, &fail, &result,  &err, &prob)
            sig_off()
    
    
    
    if method=='Cuhre':
            foundMethod=True
            sig_on()
            Cuhre(ndim,ncomp,F,userdata,eps_rel,eps_abs, flags, 
                  mineval, maxeval,
                  C_key, statefilePass, &nregions, 
                  &neval,&fail, &result,  &err, &prob)
            sig_off()
    
    
    
    cdef int nxg=D_ldxgiven*D_ngiven
    cdef double *xgiven
    
    if method=='Divonne':
            foundMethod=True
            xgiven = <double *> malloc(nxg * sizeof(double))
            for i in range(D_ldxgiven*D_ngiven):
                xgiven[i] = D_xgiven[i]
            sig_on()
            Divonne(ndim, ncomp, F, userdata, eps_rel, eps_abs, flags, seed, 
                    mineval, maxeval,
                    D_key1, D_key2, D_key3, D_maxpass, D_border, D_maxchisq,
                    D_mindeviation,
                    D_ngiven, D_ldxgiven, xgiven, 0, NULL, statefilePass,
                    &nregions, &neval, &fail, &result, &err, &prob)
            sig_off()
            free(xgiven) 
    
    if not foundMethod:
        raise Exception('The chosen integration method does not exist!')
    
    return result, err, prob, neval,fail,method, ndim
    













