r"""
Multiple Zeta Values or Euler-Zagier numbers

It computes multiple zeta values with very high accuracy for the
composition input and also helps to find any possible relations
between the multiple zeta values.

AUTHOR:

- Akhilesh P., IMSc Chennai, INDIA (2016-08-22): initial version

References:

P. Akhilesh, *Double tails of multiple zeta values*, Journal of Number Theory, Volume 170, January 2017, Pages 228-249
http://www.sciencedirect.com/science/article/pii/S0022314X16301718
"""
# ****************************************************************************
#       Copyright (C) 2013 Akhilesh P. <akhileshp.clt@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from __future__ import print_function

from sage.functions.log import log
from mpmath import mpf, pslq, mp
from sage.rings.real_mpfi import RealIntervalField
from sage.rings.real_mpfr import RealField
from sage.functions.other import factorial

#converting composition to binary word
def comptobin(a):
    r"""
    This program converts composition to binary word.

    INPUT:

    array corresponding to the composition

    OUTPUT:

    array corresponding to the binary word

    EXAMPLES::

        sage: from sage.functions.multiple_zeta import *
        sage: comptobin([2,3,1])
        [0, 1, 0, 0, 1, 1]
        sage: comptobin([2,8,7])
        [0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1]
    """
    word = []
    for i in range(len(a)):
        word += [0] * (a[i]-1) + [1]
    return word


# finding dual of the word
def dual(a):
    r"""
    Return the dual of the given word.

    EXAMPLES::

        sage: from sage.functions.multiple_zeta import *
        sage: dual([0,1,1])
        [0, 0, 1]
        sage: dual([0,1,0,1,0,0,1])
        [0, 1, 1, 0, 1, 0, 1]
    """
    b = list(reversed(a))
    for i in range(len(b)):
        b[i] = 1 - b[i]
    return b


# The program compute multiplRIF=RealIntervalField(D)e zeta value.
def multizeta(a, n, D):
    r"""
    Compute multiple zeta values of the given composition a
    and with n steps iteration each steps computed with accuracy D.

    EXAMPLES::

        sage: from sage.functions.multiple_zeta import *
        sage: multizeta([2],170,100)
        1.6449340668482264364724151666460251892189499012067984377355582293700074704032008738336289006197587?
        sage: multizeta([2,3],170,100)
        0.7115661975505724320969738060864026120925612044383392364922224964576860857450582651154252344636008?
        sage: multizeta([2,1],170,100)
        1.2020569031595942853997381615114499907649862923404988817922715553418382057863130901864558736093353?
    """
    D=D*log(10)/log(2)
    RIF=RealIntervalField(D)
    w=comptobin(a)
    wd=dual(w)
    Z=0
    U=RIF('1')/2
    B=[]
    S=[]
    count=-1
    k=len(w)
    for i in range(k):
            B.append(RIF('0'))
            if(w[i]==1 and i<k-1):
                    S.append(RIF('0'))
                    count=count+1
    Bd=[]
    Sd=[]
    countd=-1
    kd=len(wd)
    for i in range(kd):
            Bd.append(RIF('0'))
            if(wd[i]==1 and i<kd-1):
                    Sd.append(RIF('0'))
                    countd=countd+1

    for m in range(0,n+1):
            B[k-1]=RIF('1')/(m+1)
            j=count
            for i in range(k-2,-1,-1):
                    if(w[i]==0):
                            B[i]=B[i+1]/(m+1)
                    elif(w[i]==1):
                            B[i]=S[j]/(m+1)
                            S[j]=S[j]+B[i+1]
                            j=j-1
            Bd[kd-1]=RIF('1')/(m+1)
            jd=countd
            for i in range(kd-2,-1,-1):
                    if(wd[i]==0):
                            Bd[i]=Bd[i+1]/(m+1)
                    elif(wd[i]==1):
                            Bd[i]=Sd[jd]/(m+1)
                            Sd[jd]=Sd[jd]+Bd[i+1]
                            jd=jd-1
            P=0
            m1=m+1
            F=B
            Fd=Bd
            for i in range(1,len(F)):
                    if(w[i-1]==1 and w[i]==0):
                            L=1
                    if((w[i-1]==0 and w[i]==0)or (w[i-1]==1 and w[i]==1)):
                            L=2
                    if(w[i-1]==0 and w[i]==1):
                            L=3
                    P=P+F[i]*Fd[len(w)-i]*L
            Z=Z+P*U
            U=U*(m1+1)/(4*m1+2)
    return Z


def bintocomp(a):
    r"""
    This program converts binary word to composition

    Input- array  correspond to the binery word

    Output- array correspond to the composition:

    EXAMPLES::

        sage: from sage.functions.multiple_zeta import *
        sage: bintocomp([0,1,0,1,0,1])
        [2, 2, 2]
        sage: bintocomp([0,1,0,1,0,1,0,0,0,0,1])
        [2, 2, 2, 5]
    """
    b = []
    count = 1
    for j in range(len(a)):
        if a[j] == 0:
            count += 1
        else:
            b.append(count)
            count = 1
    return b


def numtocomp(n):
    r"""
    Convert a number to a composition.

    EXAMPLES::

        sage: from sage.functions.multiple_zeta import *
        sage: numtocomp(5)
        [3, 1]
        sage: numtocomp(70)
        [5, 1, 2]
    """
    if n < 0:
        raise ValueError
    if n != 0:
        w = [int(bin(n)[j]) for j in range(2, len(bin(n)))]
        w = [0] + w[1:len(w)] + [1]
        return bintocomp(w)
    else:
        return []


def IMF(N):
    r"""
    Return intial, middle and final words of all first N-1 words.

    EXAMPLES::

        sage: IMF(10)
        ([0, 0, 0, 1, 0, 2, 1, 3, 0, 4],
        [0, 0, 0, 0, 0, 1, 0, 0, 0, 2],
         [0, 0, 1, 0, 2, 3, 1, 0, 4, 5])
        sage: IMF(20)
        ([0, 0, 0, 1, 0, 2, 1, 3, 0, 4, 2, 5, 1, 6, 3, 7, 0, 8, 4, 9],
        [0, 0, 0, 0, 0, 1, 0, 0, 0, 2, 1, 3, 0, 1, 0, 0, 0, 4, 2, 5],
        [0, 0, 1, 0, 2, 3, 1, 0, 4, 5, 6, 7, 2, 3, 1, 0, 8, 9, 10, 11])
    """
    Int=[0]
    Fin=[0]
    for n in range(1,N):
            t=2
            while(n%t==0):
                    t=2*t
            t=t/2
            Int=Int+[((n/t)-1)/2]
            t=1
            while(n/t>1):
                    t=t*2
            n1=n
            mark=0
            while((n1/t)==1):
                    n1=n1-t
                    if(t/2>0):
                            t=t/2
                    else:
                            mark=1
            if(mark==1):
                    t=t/2
            Fin=Fin+[n1+t]
    Mid=[Fin[Int[i]] for i in range(N)]
    return(Int,Mid,Fin)


def allmultizeta(N, n, D):
    r"""
    Compute the first N-1 multiple zeta values.

    EXAMPLES::

        sage: from sage.functions.multiple_zeta import *
        sage: allmultizeta(10,170,100)
        [0.50000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000?,
        1.644934066848226436472415166646025189218949901206798437735558229370007470403200873833628900619758706?,
        1.202056903159594285399738161511449990764986292340498881792271555341838205786313090186455873609335258?,
        1.202056903159594285399738161511449990764986292340498881792271555341838205786313090186455873609335258?,
        1.082323233711138191516003696541167902774750951918726907682976215444120616186968846556909635941699917?,
        0.270580808427784547879000924135291975693687737979681726920744053861030154046742211639227408985424980?,
        0.811742425283353643637002772405875927081063213939045180762232161583090462140226634917682226956274938?,
        1.082323233711138191516003696541167902774750951918726907682976215444120616186968846556909635941699918?,
        1.036927755143369926331365486457034168057080919501912811974192677903803589786281484560043106557133337?,
        0.0965511599894437344656455314289427640320103723436914152525630787528921454259587614177018405925170654?]
    """
    L=[0,2,3,3]
    t=4
    t1=2
    count=0
    for i in range(4,N):
            L=L+[t1+2]
            count=count+1
            if(count==t):
                    t=2*t
                    t1=t1+1
                    count=0
    DD = D*log(10)/log(2)
    RIF = RealIntervalField(DD)
    Z = [RIF('0') for i in range(N)]
    Imf = IMF(N)
    nf = factorial(n)
    n2f = factorial(2*n)
    for i in range(n-1,-1,-1):
            Z[0] = RIF((nf*nf))/(n2f)
            for j in range(N-1,0,-1):
                    Z[j]=Z[j]+(RIF('1')/((i+1)**(L[j]-L[Imf[0][j]])))\
                    *Z[Imf[0][j]]+\
                    (RIF('1')/(i+1)**(L[j]-L[Imf[1][j]]))*Z[Imf[1][j]]+\
                    (RIF('1')/(i+1)**(L[j]-L[Imf[2][j]]))*Z[Imf[2][j]]
            if(i>0):
                    nf=nf/(i+1)
                    n2f=n2f/(2*(i+1)*(2*(i+1)-1))

    return Z


def allmultizetaprint(N, n, D):
    r"""
    Print the first N-1 multiple zeta values.

    EXAMPLES::

        sage: from sage.functions.multiple_zeta import *
        sage: allmultizetaprint(10,170,100)
        multizeta( [2] )= 1.644934066848226436472415166646025189218949901206798437735558229370007470403200873833628900619758706?
        multizeta( [3] )= 1.202056903159594285399738161511449990764986292340498881792271555341838205786313090186455873609335258?
        multizeta( [2, 1] )= 1.202056903159594285399738161511449990764986292340498881792271555341838205786313090186455873609335258?
        multizeta( [4] )= 1.082323233711138191516003696541167902774750951918726907682976215444120616186968846556909635941699917?
        multizeta( [3, 1] )= 0.270580808427784547879000924135291975693687737979681726920744053861030154046742211639227408985424980?
        multizeta( [2, 2] )= 0.811742425283353643637002772405875927081063213939045180762232161583090462140226634917682226956274938?
        multizeta( [2, 1, 1] )= 1.082323233711138191516003696541167902774750951918726907682976215444120616186968846556909635941699918?
        multizeta( [5] )= 1.036927755143369926331365486457034168057080919501912811974192677903803589786281484560043106557133337?
        multizeta( [4, 1] )= 0.0965511599894437344656455314289427640320103723436914152525630787528921454259587614177018405925170654?
    """
    A = allmultizeta(N, n, D)
    for i in range(1, N):
        print("multizeta(", numtocomp(i), ")=", A[i])


def Li(word, D):
     r"""
     This intermediate program computes the polylogarithm at 1/2.
     """
     DD=D+int(log(D)/log(10))+4
     RIF=RealIntervalField(DD)
     n=int(DD*log(10)/log(2))+1
     B=[]
     L=[]
     S=[]
     count=-1
     k=len(word)
     for i in range(k):
          B.append(RIF('0'))
          L.append(RIF('0'))
          if(word[i]==1 and i<k-1):
               S.append(RIF('0'))
               count=count+1
     T=RIF('1')
     for m in range(n):
          T=T/2
          B[k-1]=RIF(1)/(m+1)
          j=count
          for i in range(k-2,-1,-1):
               if(word[i]==0):
                    B[i]=B[i+1]/(m+1)
               elif(word[i]==1):
                    B[i]=S[j]/(m+1)
                    S[j]=S[j]+B[i+1]
                    j=j-1
               L[i]=T*B[i]+L[i]
          L[k-1]=T*B[k-1]+L[k-1]
     return(L)


def mzeta(a, D=1000):
    r"""
    Compute multiple zeta values for any given precision.

    INPUT:

    - a -- composition

    - D -- number of digits ?

    EXAMPLES::

        sage: from sage.functions.multiple_zeta import *
        sage: mzeta([2,1])
        1.2020569031595942853997381615114499907649862923404988817922715553418382057863130901864558736093352581461991577952607194184919959986...
        sage: mzeta([2,1],100)
        1.20205690315959428539973816151144999076498629234049888179227155534183820578631309018645587360933526?
    """
    a = comptobin(a)
    D = D*log(10)/log(2)
    DD = D+int(log(D)/log(10))+4
    RIF = RealIntervalField(DD)
    b = dual(a)
    l1 = Li(a, D) + [1]
    l2 = Li(b, D) + [1]
    Z = RIF('0')
    for i in range(len(l1)):
        Z += l1[i] * l2[len(a)-i]
    return Z


def Rmultizeta(a, D=100, M=100, m=10000):
    r"""
    This program allows you to find the linear relationship
    between the multiple zeta values:

    EXAMPLES::

        sage: from sage.functions.multiple_zeta import *
        sage: Rmultizeta([[2,1],[3]])
        [1, -1]
        sage: Rmultizeta([[2,1],[3]],100,100,1000)
        [1, -1]
    """
    R = RealField()
    Z = []
    mp.dps = D
    for i in range(len(a)):
        r = multizeta(a[i],int(R(1.7)*D),int(R(3.14)*D))
        Z.append(mpf((r.lower()+r.upper())/2))
    u=pslq(Z,tol=mpf(10)**-D,maxcoeff=M,maxsteps=m)
    return u
