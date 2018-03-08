# ****************************************************************************
#       Copyright (C) 2018 Christian Stump <christian.stump@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
# ****************************************************************************
from operator import itemgetter

from sage.matrix.matrix0 cimport Matrix
from sage.graphs.bliss import canonical_form_from_edge_list
from __builtin__ import frozenset

cpdef matrix_to_edge_list(Matrix M):
    r"""
    Return a bliss digraph of an integer matrix

    INPUT:

    - ``M`` (an integer (n+m) x m matrix)
    """
    cdef Py_ssize_t i, j
    cdef int a,b, x

    cdef list Vout = list()
    cdef list Vin  = list()
    cdef list labels = list()
    cdef dict edge_labels = dict()

    cdef bint pos

    cdef int n = M._ncols
    cdef int m = M._nrows - n

    cdef int Lnr = 1
    cdef int lab = 0

    for i from 0 <= i < n:
        for j from 0 <= j < i:
            a = M.get_unsafe(i,j)
            b = M.get_unsafe(j,i)
            if a == 0:
                pass
            elif a == 1 and b == -1:
                Vout.append(i)
                Vin.append(j)
                labels.append(0)
            elif a == -1 and b == 1:
                Vout.append(j)
                Vin.append(i)
                labels.append(0)
            else:
                if a > 0:
                    pos = True
                else:
                    pos = False
                try:
                    if pos:
                        lab = edge_labels[(a,b)]
                    else:
                        lab = edge_labels[(b,a)]
                except KeyError:
                    lab = Lnr
                    Lnr = Lnr+1
                    if pos:
                        edge_labels[(a,b)] = lab
                    else:
                        edge_labels[(b,a)] = lab
                finally:
                    if pos:
                        Vout.append(i)
                        Vin.append(j)
                        labels.append(lab)
                    else:
                        Vout.append(j)
                        Vin.append(i)
                        labels.append(lab)

    for i from n <= i < n+m:
        for j from 0 <= j < n:
            a = M.get_unsafe(i,j)
            if a == 0:
                pass
            elif a == 1:
                Vout.append(i)
                Vin.append(j)
                labels.append(0)
            elif a == -1:
                Vout.append(j)
                Vin.append(i)
                labels.append(0)
            else:
                if a > 0:
                    pos = True
                else:
                    pos = False
                b = -a
                try:
                    if pos:
                        lab = edge_labels[(a,b)]
                    else:
                        lab = edge_labels[(b,a)]
                except KeyError:
                    lab = Lnr
                    if pos:
                        edge_labels[(a,b)] = lab
                    else:
                        edge_labels[(b,a)] = lab
                finally:
                    if pos:
                        Vout.append(i)
                        Vin.append(j)
                        labels.append(lab)
                    else:
                        Vout.append(j)
                        Vin.append(i)
                        labels.append(lab)

    lab_relabels = [0] + [ lab for _,lab in sorted(edge_labels.iteritems(), key=itemgetter(0)) ]
    labels = [lab_relabels[i] for i in labels]

    return Vout, Vin, Lnr, labels

cpdef matrix_canonical_hash_cython(Matrix M):
    cdef int n = M._ncols
    cdef int m = M._nrows-n
    cdef list Vout
    cdef list Vin
    cdef list labels
    cdef int Lnr
    cdef list partition

    Vout, Vin, Lnr, labels = matrix_to_edge_list(M)

    #   def canonical_form_from_edge_list(int Vnr, list Vout, list Vin, int Lnr=1, list labels=[], list partition=[], bint directed=False, bint certificate=False):
    edges = canonical_form_from_edge_list(n+m, Vout, Vin, Lnr, labels, partition=[list(range(n)),list(range(n,n+m))], directed=True, certificate=False)
    return hash(frozenset(edges))
