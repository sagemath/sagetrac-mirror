# coding: UTF-8
r"""
Rabin information dispersal

Implements the Rabin information dispersal [Rabin1989]_ protocol
originally intended for secure distributed data storage. From
a cryptographic standpoint it is clear, that it cannot be secure,
however, the idea inspired many new protocols. 
The protocol relies on the idea of Reed-Solomon codes as well as
to Shamir secret sharing [Shamir1979]_.

AUTHORS:

- Thomas Loruenser (2013): initial version

REFERENCES:

.. [Rabin1989] Rabin, M. O. (1989). Efficient dispersal of information for security, 
   load balancing, and fault tolerance. Journal of the ACM, 36(2), 335–348. 
   :doi:`10.1145/62044.62050`
"""
###############################################################################
# Copyright 2013, Thomas Loruenser <thomas.loruenser@ait.ac.at>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
###############################################################################

from shamir_ss import ShamirSS

class RabinIDS(ShamirSS):
    r"""
    Rabin information dispersal.
   
    This class implements the original version of information dispersal as proposed
    by Rabin [Rabin1989]_.  This is a very basic implementation intended for 
    educational purposes only.

    EXAMPLES::

        sage: from sage.crypto.smc.rabin_ids import RabinIDS

    Generate shares::

            sage: n, k, order = 7, 3, 2**8
            sage: data = [i for i in range(15)]
            sage: ids = RabinIDS(n, k, order)
            sage: shares = ids.share(data)

    Reconstruct data::

            sage: data == ids.reconstruct(shares)
            True
    """

    ### begin module private api

    def _latex_(self):
        r"""
        Return Latex representation of self.

        EXAMPLES::

            sage: from sage.crypto.smc.rabin_ids import RabinIDS
            sage: ids = RabinIDS()
            sage: latex(ids)
            `(7,3)`-Rabin information dispersal over the field `\Bold{F}_{2^{8}}`
        """
        from sage.misc.latex import latex
        return "`({},{})`-Rabin information dispersal over the field `{}`".format(
            self._n, self._k, latex(self._F))


    def _repr_(self):
        r"""
        Return String representation of self.

        EXAMPLES::

            sage: from sage.crypto.smc.rabin_ids import RabinIDS
            sage: ids = RabinIDS()
            sage: print(ids)
            (7,3)-Rabin information dispersal over Finite Field in a of size 2^8
        """
        return "({},{})-Rabin information dispersal over {}".format(self._n, self._k, 
                                                              self._F)
    ### begin public api

    def reconstruct(self, shares, decoder='lg'):
        r"""
        Reconstruct shares.

        INPUT:

        - ``shares`` -- a list of shares ((x,y)-tuples of integer) or list of it.
        - ``decoder`` -- (default: ``'lg'``) decoder used to reconstruct secret. Must
            be one of the supported types ``'lg'`` or ``'bw'``.

        OUTPUT:

        The reconstructed data.

        EXAMPLES::

            sage: from sage.crypto.smc.rabin_ids import RabinIDS

            sage: n, k, order = 7, 3, 2**8
            sage: data = [i for i in range(15)]
            sage: ids = RabinIDS(n, k, order)
            sage: shares = ids.share(data)
            sage: data == ids.reconstruct(shares)
            True

        """
        # make shares iterable
        if type(shares[0]) == tuple:
            shares = [shares]

        # set decoder
        if decoder == 'lg':
            decode = self._rec_lagrange
        elif decoder == 'bw':
            decode = self._rec_berlekamp_welsh
        else:
            raise ValueError("unknown decoder.")

        # reconstruct data
        secret = []
        for element in shares:
            # convert to field
            points = [(self._to_GF(x), self._to_GF(y)) for x,y in element]
            # call decoder
            secret.extend([self._to_Int(i) for i in decode(points)])
        return secret


    def share(self, secret):
        r"""
        Generate shares.

        A polynomial of degree `k-1` is generated from input data.
        It is then evaluated at points starting from `1`.

        INPUT:

        - ``secret`` -- the data to be shared as list of integer.

        OUTPUT:

        The shares.

        EXAMPLES::

            sage: from sage.crypto.smc.rabin_ids import RabinIDS

            sage: n, k, order = 7, 3, 2**8
            sage: data = [i for i in range(15)]
            sage: ids = RabinIDS(n, k, order)
            sage: shares = ids.share(data)
            sage: data == ids.reconstruct(shares)
            True
        """
        # check input list size (padding is not supported)
        if len(secret)%self._k:
            raise TypeError("input list must be multiple of k (padding is not supported).")
        
        # generate shares
        shares = []
        for i in range(len(secret)/self._k):
            poly = 0
            for j in range(self._k):
                # gen polynomial form data
                poly += self._to_GF(secret[i*self._k + j]) * self._P.gen()**j
            # evaluate polynomial at different points (shares)
            shares.append([(i, self._to_Int(poly(self._to_GF(i)))) for i in range(1, self._n+1)])
        return shares


# vim: set fileencoding=UTF-8 filetype=python :
