r"""
 This is the library for Multiple Zeta Values or Euler-Zagier numbers
  
  It computes multiple zeta values with very high accuracy for the composition input and also 
  helps to find any possible relations between the multiple zeta values.
   
AUTHOR:

- Akhilesh P., HRI Allahabad, INDIA (2015-03-19): initial version


EXAMPLES::

"""

#*****************************************************************************
#       Copyright (C) 2013 Akhilesh P. <akhileshp.clt@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
import math
from mpmath import*
import itertools

r'''This program converts binery word to composition
Input- array  correspond to the binery word 
Output- array correspond to the composition
Example::
sage: bintocomp([0,0,0,1,0,1])
[4, 2]
sage: bintocomp([0,0,0,1,0,7])
'error'
'''
def bintocomp(a):
     #checking the given array is binery word
     for i in range(len(a)):
          if(a[i]!=1 and a[i]!=0):
               return("error")
     #converting the binery word to composition
     b=[]
     count=1
     for j in range(len(a)):
          if(a[j]==0):
               count=count+1
          else:
               b.append(count)
               count=1     
     return(b)



r''' This program converts binery word to number
Input- binery word as an array
output- the corresponding number
Example::
sage: bintonum([])
0
sage: bintonum([0,1])
1
sage: bintonum([0,0,1])
2
sage: bintonum([0,1,1])
3
sage: bintonum([0,1,0])
'error'
sage: bintonum([1,1,0])
'error'
sage: bintonum([1,2,0])
'error'
sage: 


'''

def bintonum(a):
     # if the word is empty return 0
     if(a==[]):
          return(0)
     # checking that the given word is binery word
     for i in range(len(a)):
          if(a[i]!=0 and a[i]!=1):
               return('error')
     # checking that given word is admissible
     if(a[0]!=0 or a[len(a)-1]!=1):
          return('error')
     # converting word to respective number
     b=a[1:len(a)-1]
     b=[1]+b
     t=1
     s=0
     for i in range(len(b)-1,-1,-1):
          s=s+t*b[i]
          t=t*2
     return(s)
          
r''' this program convert composition to binery word
Input- composition as an array
Output- binery word as an array
Example::
sage: comptobin([])
[]
sage: comptobin([2,3])
[0, 1, 0, 0, 1]
sage: 


'''
def comptobin(a):
     word=[]
     for i in range(len(a)):
          word=word+[0]*(a[i]-1)+[1]
     return(word)
r''' This program computes dual of given binery word
Input- binery word as an array
Output- dual of binery word as an array
Example::
sage: dual([0,1,1])
[0, 0, 1]
sage: dual([0,1,0,0,1])
[0, 1, 1, 0, 1]
sage: dual([0,1,0,0,4])
'error'
sage: 

'''
def dual(a):
     #checking the given array is binery word
     for i in range(len(a)):
          if(a[i]!=1 and a[i]!=0):
               return("error")
     #computing the dual
     b=list()
     b=a
     b=b[::-1]
     for i in range(len(b)):
          b[i]=1-b[i]               
     return(b)


def Li(word,D):
     DD=D+int(math.log(D)/math.log(10))+4
     mp.dps=DD
     n=int(DD*math.log(10)/math.log(2))+1
     B=[]
     L=[]
     S=[]
     count=-1
     k=len(word)
     for i in range(k):
          B.append(mpf('0'))
          L.append(mpf('0'))
          if(word[i]==1 and i<k-1):
               S.append(mpf(0))
               count=count+1
     T=mpf(1)
     for m in range(n):
          T=T/2
          B[k-1]=mpf(1)/(m+1)
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
r''' This program computes multiple zeta values for any given precision 
Example::
sage: mzeta([2,1])
mpf('1.2020569031595942853997381615114499907649862923404988817922715553418382057863130901864558736093352581461991577952607
19418491995998673283213776396837207900161453941782949360066719191575522242494243961563909664103291159095780965514651279918
40510571525598801543710978110203982753256678760352233698494166181105701471577863949973752378527793703095602570185318279000
30765471075630488433208697115737423807934450316076253177145354444118311781822497185263570918244899879620350833575617202260
339378587032813126780799005417734869115253706562370574409662217129026273207323614922429130405285553723410330775777980642420
243048828152100091460265382206962715520208227433500101529480119869011762595167636699817183557523488070371955574234729408359
5208861666202572853755813079282586487282173705566196898952662018776810629200817792338135876828426412432431480282173674506720
69350762689530434593937503296636377575062473323992348288310773390527680200757984356793711505090050273660471140085335034364672
24856531518117766181092227918431')
sage: mzeta([2,5],100)
mpf('0.6587533875711093581412522186346254271044356998380703541143384794612078112362164544354656566174209515505697938')
sage: mzeta([2],100)
mpf('1.644934066848226436472415166646025189218949901206798437735558229370007470403200873833628900619758705304004264')
sage: zeta([3],100)-zeta([2,1],100)
mpf('0.0')
sage: 


'''

def mzeta(a,D=1000):
     a=comptobin(a)
     DD=D+int(math.log(D)/math.log(10))+4
     mp.dps=DD
     b=dual(a)
     l1=Li(a,D)+[1]
     l2=Li(b,D)+[1]
     Z=mpf('0')
     for i in range(len(l1)):
          Z=Z+l1[i]*l2[len(a)-i]
     return(Z)

r''' This program allows you to find the linear relationship
between the multiple zeta values
Example ::
sage: Rmzeta([[2,1],[3]],100,100,1000)
[1, -1]
sage: Rmzeta([[4],[2,2]],100,100,1000)
[-3, 4]
sage: 

 '''
def Rmzeta(a,D=1000,M=100,m=10000):
     Z=[]
     for i in range(len(a)):
          Z.append(mzeta(a[i],D))
     u=pslq(Z,tol=mpf(10)**-D,maxcoeff=M,maxsteps=m)
     return(u)




