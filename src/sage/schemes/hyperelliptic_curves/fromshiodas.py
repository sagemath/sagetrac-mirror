from sage.rings.all import ZZ
from sage.rings.all import PolynomialRing


#Case C2 x S4 y^2 = x^8 + 14*x^4 + 1 in  char > 7 see [MaShShVo2002].
def G3ModelsInCharFF_G48_48(J_i):
      
    FF=J_i[0].parent()
    R, (x,) = PolynomialRing(FF, 1,'x').objgens()
    f = x**8 + 14*x**4 + 1
    return [f]


#Case V8  y^2 = x^8 - 1 in  char > 7 see [MaShShVo2002].
def G3ModelsInCharFF_G32_9(J_i):
      
    FF=J_i[0].parent()
    R, (x,) = PolynomialRing(FF, 1,'x').objgens()
    f=x**8 - 1
    return [f]


# Case U6 y^2 = x (x^6 - 1) in  char > 7, see [MaShShVo2002].
def G3ModelsInCharFF_G24_5(J_i):
      
    FF=J_i[0].parent()
    R, (x,) = PolynomialRing(FF, 1,'x').objgens()
    f=x*(x**6 - 1)
    return [f]

#Case C14 y^2 = x^7 - 1 in  char > 7, see [MaShShVo2002].
def G3ModelsInCharFF_C14(J_i):
      
    FF=J_i[0].parent()
    R, (x,) = PolynomialRing(FF, 1,'x').objgens()
    f=x**7-1
    return [f]
 
#Case C2xD8 y^2 = x^8 + a4*x^4 + a0 in  char > 7
def G3ModelsInCharFF_G16_11(J_i):
      
    FF=J_i[0].parent()
    R, (x,) = PolynomialRing(FF, 1,'x').objgens()
    J2, J3, J4, J5, J6, J7, J8= J_i

    if -30*J3**2+J2**3==0:
	a4=35/2*(J5*J2+6*J4*J3)/(-30*J3**2+J2**3)
    elif J4!=0:
	a4=35/3*J5/J4
    elif J5!=0 :
	a4 = 7/3*(6*J4*J2+30*J3**2-J2**3)/J5
    else:
	print "[G3Twists] G48_48 case trapped in G3ModelsInCharFF_G16_11 by error at JI = "
    a0 = -1/140*a4**2+1/2*J2
    
    f = x**8 + a4*x**4 + a0
    return [f]
    

#Case D12 y^2 = x * (x^6 + a4*x^3 + a1)  in  char > 7 
def G3ModelsInCharFF_D12(J_i):
      
    FF =J_i[0].parent()
    R, (x,) = PolynomialRing(FF, 1,'x').objgens()
    J2, J3, J4, J5, J6, J7, J8 = J_i

    if 30*J3**2-J2**3!=0:
	a4 = 280*(-J5*J2+4*J4*J3)/(30*J3**2-J2**3)
    elif J4!=0:
	a4 = 35/3*J5/J4
    elif J5!=0:
	a4 = 7/108*(96*J4*J2+30*J3**2-J2**3)/J5
    else:
	print "[G3Twists] : G48_48 group trapped in G3ModelsInCharFF_D12()"
    a1 = 2/35*a4**2-4*J2
    f = x * (x**6 + a4*x**3 + a1)
    return [f]


# Case C2xC4 y^2 = a^2*Z^8+2*a^2*Z^6+8*a*Z^2-16; in  char > 7,
def G3ModelsInCharFF_C2xC4(J_i):

    FF =J_i[0].parent()
    R, (x,) = PolynomialRing(FF, 1,'x').objgens()
    J2, J3, J4, J5, J6, J7, J8 = J_i
    
    if 6*J4-J2**2==0:
	a = FF(196/3)
    elif 96*J4-J2**2==0:
	a = FF(-196)
    elif 147*J4-2*J2**2==0:
	if J6-2/3087*J2**3==0:
	     return [x*(x-1)*(x+1)*(x**2+1)**2]
	     
	a = FF(-84)
    else:
	
	a = 98/9*(36288*J4**2-3906*J4*J2**2+14400*J6*J2+43*J2**4)/(96*J4-J2**2)/(147*J4-2*J2**2)
    f = a**2*x**8+2*a**2*x**6+8*a*x**2-16
    return [f]
    

def G3Models(J_i):
  
      J2, J3, J4, J5, J6, J7, J8, J9, J10 = J_i;
      p=J_i[0].parent().characteristic()

      if p in [2,3,5,7]:
          print "[G3Twists] currently, the algorithm is not implemented for finite fields of char. 2,3,5,7"

	  
      # C2 x S4 : y**2 = x**8 + 14*x**4 + 1
      if J2**3-30*J3**2==0 and J4==0 and J5==0 and J6==0 and J7==0 and J8==0 and J9==0 and J10==0:
          print "Automorphism group C2 x S4, curve y^2 = x^8 + 14*x^4 + 1\n"
 	  aut="C2  x S4"   
          twists=G3ModelsInCharFF_G48_48(J_i)
          return twists, aut
      
      #V8 : y**2 = x**8 - 1  	  
      if  J3==0 and J2**2-6*J4==0  and J5==0 and J2**3+36*J6==0 and J7==0 and J2**4+420*J8==0 and J9==0 and 2520*J10 - J2**5==0:
          print "Automorphism group V8, curve y^2 =  x^8 - 1\n"
          aut ="V8"
          twists=G3ModelsInCharFF_G32_9(J_i)
          return twists,aut
   
      #U6 : y^2 = x (x^6 - 1)
      if J3==0 and J2**2-96*J4==0  and J5==0 and J2**3+2304*J6==0 and J7==0 and J2**4-17920*J8==0  and J9==0 and 430080*J10+J2**5==0:
          aut="U6" 
          twists=G3ModelsInCharFF_G24_5(J_i)
          return twists,aut 
	  

      #C14 : y^2 = x^7 - 1 */
      if J2==0 and J3==0 and J4==0 and J5==0 and J6==0 and J8==0 and J9==0 and J10==0:
            print "Automorphism group C14, curve y^2 = x^7 - 1\n"
            aut="C14"
	    twists=G3ModelsInCharFF_C14(J_i)
	    return twists,aut

      #C2 x D8 : y^2 = x^8 + a*x^4 + 1 or a*x^8 + a*x^4 + b
      if J4^3 - 3/2*J4**2*J2**2 - 20*J4*J3**2*J2 + 2/3*J4*J2**4 - 200/3*J3**4 + 40/9*J3**2*J2**3 - 2/27*J2**6==0 and J5*J3 + 3/10*J4**2 - 1/4*J4*J2**2 - J3**2*J2 + 1/30*J2**4==0 and J5*J4 - 2/3*J5*J2**2 + 5*J4*J3*J2 + 20*J3**3 - 2/3*J3*J2**3==0 and J5**2 - 6/5*J4**2*J2 - 6*J4*J3**2 + 1/5*J4*J2**3==0 and J6 -1/2*J4*J2 - 10/3*J3**2 + 1/9*J2**3==0 and J7 - 1/3*J5*J2 + J4*J3==0 and J8 - 1/70*J4**2 - 1/20*J4*J2**2 - 1/3*J3**2*J2 + 1/90*J2**4==0 and 	J9 - 1/9*J5*J2**2 + 5/6*J4*J3*J2 + 10/3*J3**3 - 1/9*J3*J2**3==0 and J10 - 1/42*J4**2*J2 - 1/21*J4*J3**2 + 1/630*J4*J2**3==0:
        print "Automorphism group C2xD8, curve y^2 = x^8 + a*x^4 + 1\n"
	aut ="C2xD8"
	twists = G3ModelsInCharFF_G16_11(J_i)
	return twists, aut



      #D12 : y^2 = x * (x^6+a*x^3+1) 
      if J4**3 - 3/32*J4**2*J2**2 - 5/64*J4*J3**2*J2 + 1/384*J4*J2**4 - 25/1536*J3**4+5/4608*J3**2*J2**3 - 1/55296*J2**6==0 and J5*J3 - 16/5*J4**2 + 1/6*J4*J2**2 + 1/24*J3**2*J2 - 1/720*J2**4==0 and J5*J4 - 1/24*J5*J2**2 - 5/24*J4*J3*J2 - 5/96*J3**3 + 1/576*J3*J2**3==0 and J5**2 - 8/15*J4**2*J2 - 1/6*J4*J3**2 + 1/180*J4*J2**3==0 and J6 - 1/8*J4*J2 - 5/96*J3**2 + 1/576*J2**3==0 and J7 - 1/12*J5*J2 - 1/6*J4*J3==0 and J8 - 26/105*J4**2 + 1/120*J4*J2**2 + 1/288*J3**2*J2 - 1/8640*J2**4==0 and J10 - 5/252*J4**2*J2 - 13/1008*J4*J3**2 + 13/30240*J4*J2**3==0 and J9 - 1/144*J5*J2**2 - 5/144*J4*J3*J2 - 5/576*J3**3 + 1/3456*J3*J2**3==0:
	print "Automorphism group D12, curve y^2 * (x^6 + a*x^3 + 1)\n"
	aut="D12"
	twists=G3ModelsInCharFF_D12(J_i)
	return twists, aut
	
      #C2 x C4 : y^2 = x * (x^2 - 1) * (x^4 + a * x^2 + 1) 
      if J3==0 and J5==0 and J7==0 and J6**2 - 11/20*J6*J4*J2 + 19/2880*J6*J2**3 - 162/125*J4**3 + 929/8000*J4**2*J2**2 -161/72000*J4*J2**4 + 1/81000*J2**6==0 and J8 + 22/135*J6*J2 + 78/175*J4**2 - 47/1350*J4*J2**2 + 2/6075*J2**4==0 and J9==0 and J10 + 173/630*J6*J4 + 7/6480*J6*J2**2 + 173/3600*J4**2*J2 - 89/32400*J4*J2**3 + 1/36450*J2**5==0:
	print "Automorphism group C2xC4, curve y^2 = x * (x^2 - 1) * (x^4 + a * x^2 + 1)\n"
	aut="C2xC4"
	twists=G3ModelsInCharFF_C2xC4(J_i)
	return twists, aut
    

      return []," "


  



def HyperellipticPolynomialsFromShiodaInvariants(J_i):
   #computes the Shiodas from the nine Shioda invariants

   return G3Models(J_i)