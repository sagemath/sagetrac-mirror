r"""
Decoder for binary Goppa codes whose defining polynomials are square free which uses Patterson decoding algorithm in order to correct errors in words.
The Patterson decoder makes use of the terminated extended Euclidean algorithm.
For further details see 

N. Patterson: The algebraic decoding of Goppa codes, IEEE Transactions on
Information Theory 21.2 (1975), pp. 203-207.

AUTHORS: 
    - Giuseppe Cesare (2021-04-24)
    - Ferdinando Zullo (2021-04-24)
"""

from sage.coding.decoder import Decoder
from sage.coding.linear_code import AbstractLinearCode

class PattersonDecoder(Decoder):

    def __init__(self, code):

        r"""
        TESTS:

        If ``code`` is not a Binary Goppa code, an error is raised::

            sage: F = GF(64)
            sage: R.<x> = F[]
            sage: t = 5
            sage: a = F.gen()
            sage: f = x^2
            sage: L = [a for a in F.list() if f(a)!= 0]
            sage: n = len(L)
            sage: C = codes.GoppaCode(f,L)
            sage: D = C.decoder("PattersonDecoder")

            Traceback (most recent call last):
            ...
            ValueError: code has to be a binary Goppa and Goppa polynomial has to be squarefree
        """

        if not isinstance(code, sage.coding.goppa_code.GoppaCode) or 2 != code.base_ring().characteristic() or not code._generating_pol.is_squarefree():
            raise ValueError("code has to be a binary Goppa and Goppa polynomial has to be squarefree")
        super(PattersonDecoder, self).__init__(code, code.ambient_space(), "GoppaEncoder")

    def __eq__(self, other):
        r"""
        Tests equality of PattersonDecoder objects.

        EXAMPLES::

            sage: F = GF(64)
            sage: R.<x> = F[]
            sage: t = 5 
            sage: a = F.gen() 
            sage: f = x^5 + a^3*x^3 + (a^3 + a^2 + a)*x^2 + (a^3 + a)*x + a^4 + a^3 + a + 1
            sage: L = [a for a in F.list() if f(a)!= 0]
            sage: n = len(L)
            sage: C = codes.GoppaCode(f,L)
            sage: D1 = C.decoder("PattersonDecoder")
            sage: D2 = C.decoder("PattersonDecoder")
            sage: D1.__eq__(D2)
            True
        """ 
        return isinstance(other, PattersonDecoder) and self.code() == other.code() and self.input_space() == other.input_space()

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: F = GF(64)
            sage: R.<x> = F[]
            sage: t = 5 
            sage: a = F.gen() 
            sage: f = x^5 + a^3*x^3 + (a^3 + a^2 + a)*x^2 + (a^3 + a)*x + a^4 + a^3 + a + 1
            sage: L = [a for a in F.list() if f(a)!= 0]
            sage: n = len(L)
            sage: C = codes.GoppaCode(f,L)
            sage: D = C.decoder("PattersonDecoder")
            sage: D
            Patterson decoder decoder for [64, 34] Goppa code over GF(2)
        """
        return "Patterson decoder for %s" % self.code()

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: F = GF(64)
            sage: R.<x> = F[]
            sage: t = 5 
            sage: a = F.gen() 
            sage: f = x^5 + a^3*x^3 + (a^3 + a^2 + a)*x^2 + (a^3 + a)*x + a^4 + a^3 + a + 1
            sage: L = [a for a in F.list() if f(a)!= 0]
            sage: n = len(L)
            sage: C = codes.GoppaCode(f,L)
            sage: D = C.decoder("PattersonDecoder")
            sage: latex(D)
            \textnormal{Patterson decoder for }[64, 34]\text{ Goppa code over }\Bold{F}_{2}
        """
        return "\\textnormal{Patterson decoder for }%s" % self.code()._latex_()

    def split(self, p):

        r"""
        When considering a polynomial p(x) in GF(2^n)[x] the algorithm returns p0(x) and p1 (x) in GF(2^n)[x] such that p(x) = p0(x)^2 + x*p1(x)^2.

        INPUT:

        - ``p`` -- polynomial over ``GF(2^n)[x]``

        OUTPUT:

        - a tuple of polynomials

        EXAMPLES::

        sage: F = GF(64)
        sage: R.<x> = F[]
        sage: t = 5
        sage: a = F.gen() 
        sage: f = x^5 + a^3*x^3 + (a^3 + a^2 + a)*x^2 + (a^3 + a)*x + a^4 + a^3 + a + 1
        sage: L = [a for a in F.list() if f(a)!= 0]
        sage: n = len(L)
        sage: C = codes.GoppaCode(f,L)
        sage: D = C.decoder("PattersonDecoder")
        sage: p = R.random_element()
        sage: p
        (z6^3 + 1)*x^4 + z6^5*x^3 + (z6^3 + z6^2 + 1)*x^2 + (z6^4 + z6^3 + z6^2)*x
        sage: D.split(p)
        ((z6^5 + z6^4 + z6^3 + z6^2 + z6)*x^2 + (z6^5 + z6^4 + z6^3 + z6^2)*x,
        (z6^5 + z6^2 + 1)*x + z6^5 + z6^4 + z6^3 + 1)
        """

        Phi=p.parent()
        p0=Phi([sqrt(c) for c in p.list() [0::2]]);
        p1=Phi([sqrt(c) for c in p.list() [1::2]]);
        return (p0,p1);

    def _partial_xgcd_gen(self, a, b, cond, PolRing):

        r"""
        Performs an Euclidean algorithm on ``a`` and ``b`` until a remainder
        has degree less than or equal to `cond`, `cond` being an arbitrary parameter regarding the remainder degree which determines 
        the stop condition and returns `(r, s)` such that in the step just before termination, `r = a\times s + b\times   t`.

        INPUT:

        - ``a, b`` -- polynomials over ``PolRing``
        
        - ``cond`` -- positive integer 

        - ``PolRing`` -- polynomial ring of the output
        

        OUTPUT:

        - a tuple of polynomials

        EXAMPLES::

            sage: F = GF(64)
            sage: R.<x> = F[]
            sage: t = 5
            sage: a = F.gen() 
            sage: f = x^5 + a^3*x^3 + (a^3 + a^2 + a)*x^2 + (a^3 + a)*x + a^4 + a^3 + a + 1
            sage: L = [a for a in F.list() if f(a)!= 0]
            sage: n = len(L)
            sage: C = codes.GoppaCode(f,L)
            sage: D = C.decoder("PattersonDecoder")
            sage: P = PolynomialRing(F,'x')
            sage: x = P.parameter()
            sage: a = x^7 + x^5+ x^4 + x^2 + 1
            sage: b = x^3 + x + 1
            sage: D._partial_xgcd_gen(a, b, cond, PolRing)
            (x^2 + 1, x^4)

        """

        stop = cond 
        s = PolRing.one()
        prev_s = PolRing.zero()
        r = b
        prev_r = a
        while(r.degree() > stop):
            q = prev_r.quo_rem(r)[0]
            (prev_r, r) = (r, prev_r - q * r)
            (prev_s, s) = (s, prev_s - q * s)
        return (r, s)

    def _decode_to_code_and_message(self, r):

        r""" 
        Decodes ``r`` to an element in message space of ``self`` and its
        representation in the ambient space of the code associated to ``self``.

        INPUT:

        - ``r`` -- a vector of the ambient space of ``self.code()``

        OUTPUT:

        - ``(w, m)`` -- ``w`` is the representation of ``m`` decoded in the ambient
          space of the associated code of ``self``, ``m`` its representation in
          the message space of ``self``.

        EXAMPLES::


            sage: F = GF(64)
            sage: R.<x> = F[]
            sage: t = 5
            sage: a = F.gen() 
            sage: f = x^5 + a^3*x^3 + (a^3 + a^2 + a)*x^2 + (a^3 + a)*x + a^4 + a^3 + a + 1
            sage: L = [a for a in F.list() if f(a)!= 0]
            sage: n = len(L)
            sage: C = codes.GoppaCode(f,L)
            sage: D = C.decoder("PattersonDecoder")
            sage: c = C.random_element()
            sage: Chan = channels.StaticErrorRateChannel(C.ambient_space(), t)
            sage: y = Chan(c)
            sage: w, m = D._decode_to_code_and_message(y)
            sage: m == D.connected_encoder().unencode(c)
            True
            sage: w == c
            True
        """

        w = copy(r)
        C = self.code()
        f = C._generating_pol
        L = C._defining_set
        G = C.generator_matrix()
        R = f.parent()
        t = f.degree();
        cond = floor(t // 2)
        #checking if r is a codeword
        if w in C:
            m = G.transpose() \ w
            return w, m
        else:
            #calculation of the polynomial syndrome s
            s = 0
            for i in range(len(L)):
                s = s + w[i]*(x - L[i]).inverse_mod(f)
            T = s.inverse_mod(f)  #calculation of the inverse of s
            if T == x:     #cheching if we already got a solution
                w[0] = w[0] + 1
            else:
                f0_f1 = self.split(f) #calling the split function	
                T0_T1 = self.split(T+x)	
                radx = f0_f1[0] * f0_f1[1].inverse_mod(f) 	
                z = (T0_T1[0] + radx * T0_T1[1]).mod(f)
                #solving a0 = a1*z (mod(f))
                sols = self._partial_xgcd_gen(f,z,cond,R) #using the euclide function defined previously
                #locator polynomial construction
                #it will be = a0^2 + x * a1^2 (unless a multiplicative constant)
                a0 = sols[0]
                a1 = sols[1]
                sigma = a0**2 + x * a1**2
                #getting the roots
                for i in range(len(L)):
                    if sigma(L[i]) == R(0):
                        w[i] = w[i] + 1  #correct the error
                if w not in C:
                    raise DecodingError("Decoding failed because the number of errors exceeded the decoding radius"); 
                else:
                    m = G.transpose() \ w      
                    return w, m

    def decode_to_message(self, r):
        r"""
        Decode ``r`` to an element in message space of ``self``

        INPUT:

        - ``r`` -- a vector of the ambient space of ``self.code()``

        OUTPUT:

        - a vector of ``self`` message space

        EXAMPLES::


            sage: F = GF(64)
            sage: R.<x> = F[]
            sage: t = 5
            sage: a = F.gen() 
            sage: f = x^5 + a^3*x^3 + (a^3 + a^2 + a)*x^2 + (a^3 + a)*x + a^4 + a^3 + a + 1
            sage: L = [a for a in F.list() if f(a)!= 0]
            sage: n = len(L)
            sage: C = codes.GoppaCode(f,L)
            sage: D = C.decoder("PattersonDecoder")
            sage: c = C.random_element()
            sage: Chan = channels.StaticErrorRateChannel(C.ambient_space(), t)
            sage: y = Chan(c)
            sage: D.connected_encoder().unencode(c) == D.decode_to_message(y)
            True
        """            
        return self._decode_to_code_and_message(r)[1]

    def decode_to_code(self, r):
        r"""
        Correct the errors ``r`` and returns a codeword.

        INPUT:

        - ``r`` -- a vector of the ambient space of ``self.code()``

        OUTPUT:

        - a vector of ``self.code()``

        EXAMPLES::

            sage: F = GF(64)
            sage: R.<x> = F[]
            sage: t = 5
            sage: a = F.gen() 
            sage: f = x^5 + a^3*x^3 + (a^3 + a^2 + a)*x^2 + (a^3 + a)*x + a^4 + a^3 + a + 1
            sage: L = [a for a in F.list() if f(a)!= 0]
            sage: n = len(L)
            sage: C = codes.GoppaCode(f,L)
            sage: D = C.decoder("PattersonDecoder")
            sage: c = C.random_element()
            sage: Chan = channels.StaticErrorRateChannel(C.ambient_space(), t)
            sage: y = Chan(c)
            sage: c == D.decode_to_code(y)
            True
        """
        return self._decode_to_code_and_message(r)[0]

    def decoding_radius(self):
        r"""
        Return maximal number of errors that ``self`` can decode

        OUTPUT:

        - the number of errors as an integer

        EXAMPLES::

            sage: F = GF(2^3)
            sage: R.<x> = F[]
            sage: g = x^2 + x + 1
            sage: L = [a for a in F.list() if g(a) != 0]
            sage: C = codes.GoppaCode(g, L)
            sage: D = C.decoder("PattersonDecoder")
            sage: D.decoding_radius()
            2
        """
        return (self.code().minimum_distance()-1)//2