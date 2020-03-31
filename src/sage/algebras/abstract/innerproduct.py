# from sage.calculus.var import function
from sage.calculus.var import var
# from sage.modules.free_module import VectorSpace
from sage.modules.free_module_element import vector
from sage.functions.other import sqrt
class Error(Exception):
    """Base class for exceptions in this module."""
    pass

class IncorrectInnerProductError(Error):
    """Raised when an InnerProduct Doesn't satisfy innerproduct Properties

    Attributes:
        message -- explanation of why the specific transition is not allowed
    """

    def __init__(self, message):
        print(message)

class AbstractInnerProduct(object):
    r"""
    An Inner product on X is a mapping of XxX into the scalar field K of X.
    Construct an abstract defination of Inner Product
    Inner product must satisfy the following properties
    IP1: <x+y,z> = <x,z>+<y,z>
    IP2: <ax,y>  = a<x,y>
    IP3: <x,y>   = compliment(<y,x>)
    IP4: <x,x>   >= 0
    IP5: <x,x>   = 0 <-> x=0
	
	for every pair x,y,z in Vector field X and a on scalar filed of X
	 <.,.> returns a scalar value from scalar field of xX
    
    A inner product space is a vector space X with a inner product defined on X. 
    INPUT:

     - ``function``- a function of 2 vectors of any dimension
     - ``Vector Space``  - A Vector space element

    EXAMPLES::
        sage: def ev(x,y): return x*y
        sage: abstractInnerProduct(ev,RR)

 
    """

    def __init__(self,fn,vectorspace):
        """
        Initiates the abstract Inner Product Class.

        EXAMPLES::

        sage: def ev(x,y): return (x).dot_product(y)
        sage: VS = VectorSpace(RR,4)
        sage: g=abstractInnerProduct(ev,VS)

        It would display your inner product


        However, if your inner product space is not valid, it would raise
        error IncorrectInnerProductError with a message.

        sage: def ev(x,y): return (x.appy_map(sqrt)).dot_product(y.appy_map(sqrt))
        sage: VS = VectorSpace(RR,4)
        sage: g=abstractInnerProduct(ev,VS)
        IncorrectInnerProductError: is not <x+y,z> = <x,z>+<y,z>
        
        This implies the that the inner product defined is inccorrect


        In Complex Space, Inner product is defined as <x,y>=Sum(x_i*conjugate(y_i))
        sage: def ev(x,y): return (x).dot_product(y)
        sage: VS = VectorSpace(CC,4)
        sage: g=abstractInnerProduct(ev,VS)
        will throw you error IncorrectInnerProductError: is not <x,y> = congjugate of <y,x>

        However, if
        sage: def ev(x,y): return (x).dot_product(y.conjugate())
        sage: VS = VectorSpace(CC,4)
        sage: g=abstractInnerProduct(ev,VS)
        <x,y>= x0*conjugate(y0) + x1*conjugate(y1) + x2*conjugate(y2) + x3*conjugate(y3)
        Will be approved as Inner Product

        """
        # if not isinstance(fn, function):
            # raise TypeError("Argument must be a function")

        self.dimensions = vectorspace.dimension()
        self.field      = vectorspace.base_field()

        self.eval_function = fn
        self._check_Inner_Product_()

    def initialize_symbolic_variable(self):
        """
        Returns an instance of symbolic variable x belonging to the vector space X assosiated with this innerproduct space
        """
        xv = list(var('x%d'%i,domain=self.field) for i in range(self.dimensions))
        return vector(xv)

    def norm_defination(self):
        """
        Returns a instance of symbolic norm assosiated with inner product space defined here
        """
        xt = self.initialize_symbolic_variable()
        return sqrt(self.eval_function(xt,xt))

    def innerproduct(self,arg1,arg2):
        """
        Computes Inner Product between two Vectors in inner product space defined in the object
        INPUT:

        - ``Vector 1``- a vector on vector space defined exactly as in abstracted defination 
        - ``Vector 2``  -same as vector 1 
        """
        return self.eval_function(arg1,arg2)

    def norm(self,arg1):
        """
        Computes Norm in inner product space defined in the object
        INPUT:

        - ``Vector 1``- a vector on vector space defined exactly as in abstracted defination 
        """
        return sqrt(self.innerproduct(arg1,arg1))


    def _check_Inner_Product_(self):
        """
        A function to check whether the given function satisfies the conditions of Vector Space
        """

        xv = list(var('x%d'%i,domain=self.field) for i in range(self.dimensions))
        yv = list(var('y%d'%i,domain=self.field) for i in range(self.dimensions))
        zv = list(var('z%d'%i,domain=self.field) for i in range(self.dimensions))

        xt = vector(xv)
        yt = vector(yv)
        zt = vector(zv)
        

        try:
            k = self.eval_function(x= xt+zt,y=yt) - (self.eval_function(x=xt,y=yt)+self.eval_function(x=zt,y=yt))
            k = k.simplify_full()
            if(int(k)!=0):
                raise IncorrectInnerProductError("is not <x+y,z> = <x,z>+<y,z>")
        except Exception:
            raise IncorrectInnerProductError("is not <x+y,z> = <x,z>+<y,z>")

        a= var('a',domain = self.field)

        try:
            j = self.eval_function(x=a*xt,y=yt) - (a*self.eval_function(x=xt,y=yt))
            j = j.simplify_full()
            if(int(j)!=0):
                raise IncorrectInnerProductError("is not <a*x,y> = a<x,y>")
        except Exception:
            raise IncorrectInnerProductError("is not <a*x,y> = a<x,y>")

        try:
            l = self.eval_function(x=xt,y=yt) - self.eval_function(x=yt,y=xt).conjugate()
            l = l.simplify_full()
            if(int(l)!=0):
                raise IncorrectInnerProductError("is not <x,y> = congjugate of <y,x>")
        except Exception:
            raise IncorrectInnerProductError("is not <x,y> = congjugate of <y,x>")
        
        mn = min(0,self.eval_function(x=xt,y=xt))
        if(mn!=0):
            raise IncorrectInnerProductError("<x,x> should be greater than equal to 0")

        print("<x,y>=", self.eval_function(xt,yt))

