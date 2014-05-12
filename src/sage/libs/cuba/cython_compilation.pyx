"""
Integrand compilation using cython for the cuba integration in cuba.pyx
"""


#*****************************************************************************
#       Copyright (C) 2013 Svetlin Tassev <tassev@astro.princeton.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************


from sage.misc.preparser import preparse
from sage.misc.cython_c import cython
from re import sub as repl




def _THEFUNCTION45624634563457584g5435b6hbv564c(*args):
    """
    This global function is initialized as a dummy. It is used (as an ugly way) 
    to pass the integrand to the compiled cython code in
    ``_compile_with_cython_and_then_integrate``.
    The name is pretty random to avoid possibilities of conflicts. See
    ``_compile_with_cython_and_then_integrate``. 
    
    Since it is overwritten in the calls to the integrators in the examples 
    below, this definition must be placed above them, if you want the 
    test below to work. 
    
    
    TESTS::
    
        The function must be initialized to return a scalar and take 
        arbitrary number of arguments:
        
            sage: from sage.libs.cuba.cython_compilation import _THEFUNCTION45624634563457584g5435b6hbv564c
            sage: x,y,a,z=var('x y a z')
            sage: _THEFUNCTION45624634563457584g5435b6hbv564c(x,y,a,x,z)
            0.0
        
    """
    return 0.0 



def _compile_with_cython_and_then_integrate(func, vals, jac, num_vars, 
                      additional_imports, replace_dictionary, preprocess_string,
                      verbose, 
                      cuba_verbose, only_inspect,**kwargs):
    '''
    Attempt to compile a function with cython
    and then to integrate it. The input and output are exactly the same 
    as the ones of ``multidim_integration``. Also most tests of that
    function actually test this one. So, see that function for more.
    
    TESTS::
    
    These tests look like a mess because of the required random
    integration variable names, which must not conflict anything 
    a sane person would choose.
    
    Test with a Python function::
    
        sage: from sage.libs.cuba.cython_compilation import _compile_with_cython_and_then_integrate
        sage: new_vars_BEgIn423N0eN4324QxdD,new_vars_BEgIn423N1eN4324QxdD=var('new_vars_BEgIn423N0eN4324QxdD,new_vars_BEgIn423N1eN4324QxdD')
        sage: def func(new_vars_BEgIn423N0eN4324QxdD,new_vars_BEgIn423N1eN4324QxdD):
        ....:     return new_vars_BEgIn423N0eN4324QxdD+new_vars_BEgIn423N1eN4324QxdD^2
        sage: _compile_with_cython_and_then_integrate(func,[new_vars_BEgIn423N0eN4324QxdD,new_vars_BEgIn423N1eN4324QxdD],new_vars_BEgIn423N1eN4324QxdD^3,3,"",[],True,0,0,False)[0] # abs tol 3e-3
        0.2916666666666667

    Test with a constant function:
    
        sage: from sage.libs.cuba.cython_compilation import _compile_with_cython_and_then_integrate
        sage: new_vars_BEgIn423N0eN4324QxdD,new_vars_BEgIn423N1eN4324QxdD=var('new_vars_BEgIn423N0eN4324QxdD,new_vars_BEgIn423N1eN4324QxdD')
        sage: func = 1.0
        sage: _compile_with_cython_and_then_integrate(func,[new_vars_BEgIn423N0eN4324QxdD,new_vars_BEgIn423N1eN4324QxdD],new_vars_BEgIn423N1eN4324QxdD,3,"",[],True,0,0,False)[0] # abs tol 3e-3
        0.5
    
    Test with string function:
    
        sage: from sage.libs.cuba.cython_compilation import _compile_with_cython_and_then_integrate
        sage: new_vars_BEgIn423N0eN4324QxdD,new_vars_BEgIn423N1eN4324QxdD,new_vars_BEgIn423N2eN4324QxdD=var('new_vars_BEgIn423N0eN4324QxdD,new_vars_BEgIn423N1eN4324QxdD,new_vars_BEgIn423N2eN4324QxdD')
        sage: func = ("""
        ....: cpdef the_supplied_integrand(x,y,z):
        ....:     return x*y^2+z
        ....: """)
        sage: _compile_with_cython_and_then_integrate(func,[new_vars_BEgIn423N0eN4324QxdD,new_vars_BEgIn423N1eN4324QxdD,new_vars_BEgIn423N2eN4324QxdD],new_vars_BEgIn423N1eN4324QxdD,3,"",[],True,0,0,False)[0] # abs tol 3e-3
        0.375
    
    
    
    
    AUTHORS:
    
    - Svetlin Tassev (2013-12)
    
    
    '''
    
    #replace_dictionary=[]
    #Dictionary between sage and C math.h names, so that the code is fast.
    replace_dictionary.append(('pi','3.14159265358979323846264338328'))
    replace_dictionary.append(('e','exp(1)'))
    replace_dictionary.append(('abs','fabs'))
    replace_dictionary.append(('arctan','atan'))
    replace_dictionary.append(('arcsin','asin'))
    replace_dictionary.append(('arccos','acos'))
    replace_dictionary.append(('arctanh','atanh'))
    replace_dictionary.append(('arcsinh','asinh'))
    replace_dictionary.append(('arccosh','acosh'))
    
    
    
    if verbose>=2:
        print "Attempting to compile function using Cython ..."
    
    
    tmpv=str(vals)
    tmpj=str(jac)
    
    for r in replace_dictionary:
        tmpj=_replace_word(r[0],r[1],tmpj)
        tmpv=_replace_word(r[0],r[1],tmpv)
        if (preprocess_string and isinstance(func,str)):
            func=_replace_word(r[0],r[1],func)
            
        
    #The following is ugly!
    # Couldn't find a way to import division from __future__ 
    # in the cython code to 
    # fix int's division without rounding. This should be fairly safe as ranges is
    # supposed to contain only expressions. So, no vector indices 
    # for example to mess with it. If it results in a mess, the compilation should hopefully fail.
    tmpv=preparse(tmpv).replace('Integer','float').replace('RealNumber','').replace('(\'','(').replace('\')',')')
    tmpj=preparse(tmpj).replace('Integer','float').replace('RealNumber','').replace('(\'','(').replace('\')',')')
    
    x_will_be_vector="True"
    
    if num_vars>1:
        tmpv=tmpv.replace('new_vars_BEgIn423N','float(X[').replace('eN4324QxdD','])')
        tmpj=tmpj.replace('new_vars_BEgIn423N','float(X[').replace('eN4324QxdD','])')
    else:
        # make sure that we pass a non-vector argument when ndim=1, so that 
        # multidim_integration_unit_cube can use GSL
        tmpv=repl("new_vars_BEgIn423N([0-9]*)eN4324QxdD","X",tmpv)
        tmpj=repl("new_vars_BEgIn423N([0-9]*)eN4324QxdD","X",tmpj)
        x_will_be_vector="False"
        
    if (preprocess_string and isinstance(func,str)):
            func = preparse(func).replace('Integer','float').replace('RealNumber','').replace('(\'','(').replace('\')',')')
    
    
    if (verbose>=3):
        print "Jacobian = " +tmpj
        print "Function arguments expressed with unit hypercube variables: "+tmpv
    
    
    
    #Couldn't find another way to pass the function to the cython code.
    global _THEFUNCTION45624634563457584g5435b6hbv564c
    _THEFUNCTION45624634563457584g5435b6hbv564c=func
    
    if verbose>=2:
        success_string="\nprint \"Success!\""
    else:
        success_string=""
    
#############################################################################################
    if callable(func): # The indentation below is ugly as hell, but otherwise compiling fails...
        source=("""  
cargs -O3
from libc.math cimport *   
"""+additional_imports+"""
from sage.libs.cuba.cuba import multidim_integration_unit_cube
from sage.libs.cuba.cython_compilation import _THEFUNCTION45624634563457584g5435b6hbv564c
cpdef double func1(X):
    jac1="""+tmpj+"""
    vals1="""+tmpv+"""
    return _THEFUNCTION45624634563457584g5435b6hbv564c(*vals1)*jac1
"""+success_string+"""
x=multidim_integration_unit_cube(func1,x_is_vector=""" + x_will_be_vector + 
",ndim="+str(num_vars)+",verbose="+str(cuba_verbose)+",**"+str(kwargs)+")")
            
            
###################################################################################################
    elif not isinstance(func,str):
        source=("""
cargs -O3
from libc.math cimport *
"""+additional_imports+"""
from sage.libs.cuba.cuba import multidim_integration_unit_cube
cpdef double func1(X):
    jac1="""+tmpj+"""
    return """+str(float(func))+"*jac1"+success_string+"""
x=multidim_integration_unit_cube(func1,x_is_vector=""" + x_will_be_vector + 
",ndim="+str(num_vars)+",verbose="+str(cuba_verbose)+",**"+str(kwargs)+")")

###################################################################################################
    elif isinstance(func,str):
        source=("""
cargs -O3
from libc.math cimport *   
"""+additional_imports+"""
from sage.libs.cuba.cuba import multidim_integration_unit_cube

"""+func+"""

cpdef double func1(X):
    jac1="""+tmpj+"""
    vals1="""+tmpv+"""
    return the_supplied_integrand(*vals1)*jac1
"""+success_string+"""
x=multidim_integration_unit_cube(func1,x_is_vector=""" +x_will_be_vector + 
",ndim="+str(num_vars)+",verbose="+str(cuba_verbose)+",**"+str(kwargs)+")")
         
###################################################################################################
    if (verbose>=3) or (only_inspect):
        print "The source code for the compiled function to be integrated over"
        print "the unit hypercube is as follows:"
        print "-------------------------"
        print source.replace("_THEFUNCTION45624634563457584g5435b6hbv564c","the_supplied_integrand")
        print "-------------------------"
            
    if only_inspect:
        print "Only source code inspection requested."
        return ""
            
    source=source.replace("\n","\\n").replace("\"","\\\"")
    source="cython(\""+source+"\"),x"
    
    return eval(source,{'cython':cython})[1]








def _replace_word(input,out,str):
    """
    Replace word ``input`` with ``out`` in string ``str``. See 
    ``multidim_integration``.
    
    TESTS::
        
        sage: from sage.libs.cuba.cython_compilation import _replace_word
        sage: _replace_word('pi','3.14','pi=3.14, pi^2,pi*2. But multipi is not multi3.14.')
        '3.14=3.14, 3.14^2,3.14*2. But multipi is not multi3.14.'
        
    """
    return repl("\\b"+input+"\\b",out,str)
    
    

