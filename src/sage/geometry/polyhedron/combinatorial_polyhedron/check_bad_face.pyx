from cython.parallel        cimport prange, threadid
from cysignals.memory               cimport sig_free, sig_calloc
from libc.stdio                     cimport FILE, fopen, fclose, fwrite, fread
from sage.ext.memory_allocator  cimport MemoryAllocator

cdef extern from "Python.h":
    int unlikely(int) nogil  # Defined by Cython

cdef check_from_file2(size_t m, size_t n_coatoms, size_t **PolyIneq, path, size_t n_threads):
    cdef FILE *fp
    fp = fopen(path.encode('utf-8'), "r")
    if (fp==NULL):
        raise IOError("cannot open file {}".format(path))
    cdef size_t len_bad_faces
    cdef size_t n_bad_faces
    cdef MemoryAllocator mem = MemoryAllocator()
    fread(&len_bad_faces, 1, sizeof(size_t), fp)
    fread(&n_bad_faces, 1, sizeof(size_t), fp)

    cdef uint8_t *bad_faces = <uint8_t *> mem.calloc(1, sizeof(uint8_t))
    cdef uint64_t *bad_faces_LHS = <uint64_t *> mem.calloc(1, sizeof(uint64_t))
    cdef uint8_t **bad_faces_pt = <uint8_t **> mem.allocarray(1, sizeof(uint8_t*))
    cdef size_t i,j
    cdef size_t counter
    cdef size_t e
    cdef size_t len_this_face
    cdef size_t l = 0
    while len_bad_faces:
        bad_faces = <uint8_t *> mem.reallocarray(<void *> bad_faces, len_bad_faces, sizeof(uint8_t))
        bad_faces_LHS = <uint64_t *> mem.reallocarray(<void *> bad_faces_LHS, n_bad_faces, sizeof(uint64_t))
        bad_faces_pt = <uint8_t **> mem.reallocarray(<void *> bad_faces_pt, n_bad_faces, sizeof(uint8_t*))
        fread(bad_faces, len_bad_faces, sizeof(uint8_t), fp)
        fread(bad_faces_LHS, n_bad_faces, sizeof(uint64_t), fp)
        counter = 0
        for i in range(n_bad_faces):
            len_this_face = bad_faces[counter]
            bad_faces_pt[i] = &bad_faces[counter]
            counter += 2 + len_this_face


        for i in prange(n_bad_faces, nogil=True, num_threads=n_threads, schedule='dynamic'):
            par_check(m, n_coatoms, PolyIneq, bad_faces_pt[i], bad_faces_LHS[i])

        fread(&len_bad_faces, 1, sizeof(size_t), fp)
        fread(&n_bad_faces, 1, sizeof(size_t), fp)
    fclose(fp)

cdef inline int par_check(size_t m, size_t n_coatoms, size_t **PolyIneq, uint8_t *bad_face, uint64_t LHS) nogil except -1:
    cdef size_t len_this_face, e
    cdef int test
    len_this_face = bad_face[0]
    e = bad_face[1]
    test = check_bad_face_uint8_t(
            PolyIneq, n_coatoms, m,
            LHS,
            &bad_face[2], len_this_face, e)
    if unlikely(test == 1):
        with gil:
            print("Discovered non-empty face")
    elif unlikely(test == 2):
        with gil:
            print("e=", e)
            print("Hrep is {}".format(tuple(bad_face[j] for j in range(2, len_this_face))))
            raise ValueError("Wilf conjection does not hold, counterexample")
