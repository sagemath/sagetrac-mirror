# ****************************************************************************
#       Copyright (C) 2018 Christian Stump <christian.stump@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
# ****************************************************************************
from sage.matrix.matrix0 cimport Matrix

cpdef matrix_to_edge_list(Matrix M):
    r"""
    Return a bliss digraph of an integer matrix

    INPUT:

    - ``M`` (an integer (n+m) x m matrix)
    """
    cdef Py_ssize_t i, j
    cdef int a,b, x

    cdef list edges = list()
    cdef dict edge_labels = dict()
    cdef list new_partition = list()

    cdef bint pos

    cdef int n = M._ncols
    cdef int m = M._nrows - n

    cdef int Vnew = n+m

    for i from 0 <= i < n:
        for j from 0 <= j < i:
            a = M.get_unsafe(i,j)
            b = M.get_unsafe(j,i)
            if a == 0:
                pass
            elif a == 1 and b == -1:
                edges.append((i,j))
            elif a == -1 and b == 1:
                edges.append((j,i))
            else:
                if a > 0:
                    pos = True
                else:
                    pos = False
                try:
                    if pos:
                        x = edge_labels[(a,b)]
                    else:
                        x = edge_labels[(b,a)]
                except KeyError:
                    x = len(new_partition)
                    if pos:
                        edge_labels[(a,b)] = x
                    else:
                        edge_labels[(b,a)] = x
                    new_partition.append([])
                finally:
                    if pos:
                        edges.extend([(i,Vnew),(Vnew,j)])
                    else:
                        edges.extend([(j,Vnew),(Vnew,i)])
                    new_partition[x].append(Vnew)
                    Vnew += 1

    for i from n <= i < n+m:
        for j from 0 <= j < n:
            a = M.get_unsafe(i,j)
            if a == 0:
                pass
            elif a == 1:
                edges.append((i,j))
            elif a == -1:
                edges.append((j,i))
            else:
                if a > 0:
                    pos = True
                else:
                    pos = False
                b = -a
                try:
                    if pos:
                        x = edge_labels[(a,b)]
                    else:
                        x = edge_labels[(b,a)]
                except KeyError:
                    x = len(new_partition)
                    if pos:
                        edge_labels[(a,b)] = x
                    else:
                        edge_labels[(b,a)] = x
                    new_partition.append([])
                finally:
                    if pos:
                        edges.extend([(i,Vnew),(Vnew,j)])
                    else:
                        edges.extend([(j,Vnew),(Vnew,i)])
                    new_partition[x].append(Vnew)
                    Vnew += 1

    partition = [list(range(n))]
    if m > 0:
        partition.append(list(range(n,n+m)))
    if new_partition:
        partition.extend(new_partition)

    return edges, Vnew, partition
