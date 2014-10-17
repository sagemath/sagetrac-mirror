from  sage_coding_project.BlockCode import *
from  sage_coding_project.McEliece import *
from  sage_coding_project.Errors import *
import time

class ReedSolomonCode(BlockCodeTools):

    def __init__(self, F, n, k, alpha):
        self.F = F
        self.n = n
        self.k = k
        self.alpha = alpha

    def __eq__(self, other):
        return isinstance(other, ReedSolomonCode) \
                and self.F == other.F \
                and self.n == other.n \
                and self.k == other.k \
                and self.alpha == other.alpha

    def __repr__(self):
        return "[%s, %s, %s] Reed-Solomon Code over %s" % (self.n, self.k, self.minimum_distance(), self.F)

    def __latex__(self):
        return "TODO"

    def generator_matrix(self):
        r"""
        Returns a generator matrix for self.

        OUTPUT:

        - A generator matrix.

        EXAMPLES:

        sage: C = ReedSolomonCode(GF(11), 10, 5, 2)
        sage: C.generator_matrix()
        [ 1  1  1  1  1  1  1  1  1  1]
        [ 1  2  4  8  5 10  9  7  3  6]
        [ 1  4  5  9  3  1  4  5  9  3]
        [ 1  8  9  6  4 10  3  2  5  7]
        [ 1  5  3  4  9  1  5  3  4  9]
        """
        return matrix(self.F, self.k, self.n, lambda i,j : (self.alpha^j)^i)
    
    def minimum_distance(self):
        r"""
        Returns minimum distance for self.

        OUTPUT : 

        - The minimum distance of self.

        EXAMPLES:

        sage: C = ReedSolomonCode(GF(11), 10 ,5, 2)
        sage: C.minimum_distance()
        6
        """


        return self.n-self.k+1

    def encode(self, word):
        r"""
        Returns word as a codeword of self.

        INPUT:

        - ``word`` -- list of numbers of F. Must be
          of size k.

        OUTPUT:

        - A word of self in vector format.

        EXAMPLES:

        sage: C = ReedSolomonCode(GF(11), 10 , 5, 2)
        sage: word =  [1, 2, 3, 4, 5]
        sage: C.encode(word)
        (4, 8, 9, 0, 4, 3, 2, 9, 8, 7)
        """
        return vector(word)*self.generator_matrix()

    def random_codeword(self):
        r"""
        Returns a random word of self.

        OUTPUT:

        - A word of self in vector format.

        EXAMPLES

        sage: C = ReedSolomonCode(GF(11), 10 , 5, 2) 
        sage: C.random_codeword() # random
        """

        return random_vector(self.F, self.k)*self.generator_matrix()

         

    def decode_welch_berlekamp(self,r):
        r"""
        Returns decoded version of r using Welch-Berlekamp decoding algorithm.
        
        INPUT:

        - ``r`` received message.

        OUTPUT:

        - A decoded version of r. If there was no errors in r, the algorithm stops
          before computing anything. If there is too many errors it will fail and
          return an exception. 

        EXAMPLES:

        First we create a Reed-Solomon code:

        sage: C = ReedSolomonCode(GF(11), 10, 5, 2)

        Then we pick a word we encode:

        sage: r = vector(GF(11), [5, 9, 0, 6, 0, 1, 0, 7, 0, 4])
        sage: r
        (5, 9, 0, 6, 0, 1, 0, 7, 0, 4)
        
        We add some errors to the word:

        sage: r[3] = 9
        sage: r[9] = 5
        sage: r
        (5, 9, 0, 9, 0, 1, 0, 7, 0, 5)

        We can now try to decode it:
        
        sage: C.decode_welch_berlekamp(r)
        (5, 9, 0, 6, 0, 1, 0, 7, 0, 4)
        
        If we try to decode a word with more than 2 errors (maximum decoding for this code), it
        returns -1 because it can't decode it:

        sage: r[0]=7
        sage: C.decode_welch_berlekamp(r)
        -1

        And if we try to decode a word without any errors, it returns the word without computing 
        anything else than the syndrome : 

        sage: r = vector(GF(11), [5, 9, 0, 6, 0, 1, 0, 7, 0, 4])
        sage: C.decode_welch_berlekamp(r)
        (5, 9, 0, 6, 0, 1, 0, 7, 0, 4)


        REFERENCES:

        .. Hoholdt, Justesen 
           "A Course In Error Correcting Codes"
           pages 51 to 53
        """
        if self.syndrome(r) == 0:
            return r

        t = floor((self.minimum_distance()-1)/2)
        l0 = self.n-1-t
        l1 = self.n-1-t-(self.k-1)
        S = matrix(self.F,self.n,l0+l1+2,lambda i,j : (self.alpha^i)^j if j<(l0+1)\
                else r[i]*(self.alpha^i)^(j-(l0+1)))
        Sk = S.right_kernel()
        Sk = Sk.basis_matrix().row(0)
        R.<x> = self.F[]

        Q0 = R(Sk.list_from_positions([0..l0]))
        Q1 = R(Sk.list_from_positions([l0+1..l0+l1+1]))
        
        if Q1.divides(Q0) == False:
            return -1
        g = (-Q0)//Q1
        dec = vector( g(self.alpha^i) for i in range(self.n))

        return dec


    def unencode(self, c):
        points=[]
        for i in range (0, self.k):
            points.append((self.alpha^i, c[i]))
        R = PolynomialRing(self.F, 'x')
        P = R.lagrange_polynomial(points)
        m = vector(C.F, P.coeffs())
        return m
    


class DecodingError(Exception):

    def __init__(self, message):
        self.message = message

    def __str__(self):
        return self.message



##################################################################
#                                                                #
#                   RUNS AND TESTS COME HERE                     #
#                                                                #
##################################################################

#For now, I'm just performing "brutal" testing printing some
#testing values directly from the file
#Change values and comment/uncomment lines to customize
#your tests

#p=11
#F.<a> = GF(p, 'alpha')
#GF(p)
#n=10
#k=5
#alpha=2
    
#C=ReedSolomonCode(F,n,k,alpha)
#r=vector(C.F, n, [5,9,0,9,0,1,0,7,0,5])


#r=C.random_codeword()
#print "Clean codeword           ", r
#r=r+random_error_vect(C.n, C.F, random_error_pos(C.n, 2))
#print "Codeword with errors     ", r


#print C.generator_matrix()
#print "Parity"
#print C.parity_check_matrix()
#print "Codeword after decoding  ",C.decode_welch_berlekamp(r)

#total=0
#for j in range (0,100):
#
#    for i in range (0, 100):
#        r=C.random_codeword()
#        err=r+random_error_vect(C.n, C.F, random_error_pos(C.n, 3))
#        res = C.decode_welch_berlekamp(err)
#        if res != -1:
#            if res != r:
#                if C.syndrome(res) != 0:
#                    total+=1
#    
#
#print "On 10,000 tests, %d unexpected results encountered. Mean : %d percent of wrong results per block of 100 tests"%(total, total/1000000*100)


##################################################################
#                                                                #
#                      PLOTTING COMES HERE                       #
#                                                                #
##################################################################

#p=4001
#F=GF(p)
#alpha=F.primitive_element()
#
#list_n_time = []
#n=100
#
#while n<4001:
#    print "Hello"
#    C=ReedSolomonCode(F,n,floor(0.9*n),alpha)
#    t = 0
#    time_list = []
#    for i in range (0,15):
#        r=C.random_codeword()
#        r=r+random_error_vect(C.n, C.F, random_error_pos(C.n, floor((C.n-C.k)/2)))
#
#        before = time.time()
#        C.decode_welch_berlekamp(r)
#        elapsed = time.time() - before
#        time_list.append(elapsed)
#
#    t = median(time_list)
#    list_n_time.append((n,t))
#    n += 400
#
#print list_n_time
#list_plot(list_n_time, plotjoined=True, color='teal', scale='loglog')

C = ReedSolomonCode(GF(11), 10, 5, 2)
(pubkey,privkey) = key_gen(C)

m = random_vector(C.F, C.k)
print "Original message :   ", m
enc = encrypt(m, pubkey)
print "Encrypted message :  ", enc
print "Unencode encrypted : ", C.unencode(enc)
print "Deciphered message : ", decrypt(enc, privkey)


