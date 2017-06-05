#as it turns out, computing the nodes can easily turn out to be more
#expensive than computing the integrals. So it's worth optimizing this.
#making the function into a cython routine helps a little bit. If we really
#want to we can optimize this further, probably to a point where
#we don't have to bother with node computation routines that have a better order
#than this naive approach (which is quadratic)
from sage.libs.mpfr cimport *
import math
from sage.rings.real_mpfr import RealField
from sage.misc.cachefunc import cached_function
from sage.rings.real_mpfr cimport RealNumber, RealField_class

#t3,t2,t1=t2,t1,(R(2*j1-1)*r*t1 - R(j1-1)*t2)/R(j1)

@cached_function
def gl_nodes(degree,prec):
    cdef int j,j1,n
    cdef RealNumber r,t1,t2,t3,t4,a,w
    cdef mpfr_t u,v
    cdef RealField_class R
    R=RealField(int(prec*3/2))
    Rout=RealField(prec)
    mpfr_init2(u,R.__prec)
    mpfr_init2(v,R.__prec)
    ZERO=R.zero()
    ONE=R.one()
    HALF=ONE/2
    TWO=2*ONE
    rnd=R.rnd
    epsilon=R(1)>>(prec+8)
    if degree == 1:
        x=R(3)/5
        w=R(5)/18
        nodes = [((1-x)/2,w),(HALF,R(4)/9),((1+x)/2,w)]
    else:
        nodes=[]
        n=3*2**(degree-1)
        upto = n//2+1
        for j in xrange(1,upto):
            r=R(math.cos(math.pi*(j-0.25)/(n+0.5)))
            while True:
                t1,t2=ONE,ZERO
                for j1 in xrange(1,n+1):
                    mpfr_mul(u,r.value,t1.value,rnd)
                    mpfr_mul_si(u,u,2*j1-1,rnd)
                    mpfr_mul_si(v,t2.value,j1-1,rnd)
                    mpfr_sub(u,u,v,rnd)
                    mpfr_div_si(u,u,j1,rnd)
                    t2=t1
                    t1=R._new()
                    mpfr_set(t1.value,u,rnd)
                t4=R(n)*(r*t1-t2)/(r**2-ONE)
                a=t1/t4
                r=r-a
                if a.abs()<epsilon:
                    break
            x=r
            w = ONE/((ONE-r**2)*t4**2)
            nodes.append(((ONE+x)/TWO,w))
            nodes.append(((ONE-x)/TWO,w))
    nodes=[(Rout(x),Rout(w)) for x,w in nodes]
    nodes.sort()
    mpfr_clear(u)
    mpfr_clear(v)
    return nodes

def estimate_error(results,prec,epsilon):
    if len(results)==2:
        return max((results[0][i]-results[1][i]).abs() for i in xrange(len(results[0])))
    e=[]
    for i in xrange(len(results[0])):
        try:
            if results[-1][i] == results[-2][i] == results[-3][i]:
                e.append(0*epsilon)
            D1=(results[-1][i]-results[-2][i]).abs().log()
            D2=(results[-1][i]-results[-3][i]).abs().log()
        except ValueError:
            e.append(epsilon)
        D4=min(0,max(D1**2/D2,2*D1,-prec))
        e.append(D4.exp())
    return max(e)

def gl_integrate(f,prec,epsilon=None):
    results=[]
    degree=1
    #this is fishy: epsilon is an *absolute* error
    #whereas prec measures relative precision.
    #if the integral has a value that is "about one"
    #(which is essentially true for any non-zero value in multiprecision)
    #then this is not unreasonable:
    #(definitely for our purposes, we want to control absolute error
    #in the integral values, but I don't know how)
    if epsilon is None:
        epsilon = (2.0)**(-prec)
    while True:
        nodes=gl_nodes(degree,prec)
        I=nodes[0][1]*f(nodes[0][0])
        for i in xrange(1,len(nodes)):
            I+=nodes[i][1]*f(nodes[i][0])
        results.append(I)
        if degree > 1:
            err= estimate_error(results,prec,epsilon)
            print degree,err
            if err <= epsilon:
                return I
        degree +=1
########################################################################