from .list_of_faces             cimport face_list_struct, face_struct

cdef int Vrep_list_to_bit_rep(tuple Vrep_list, face_struct& output) except -1

cdef int incidences_to_bit_rep(tuple incidences, face_struct& output) except -1

cdef size_t bit_rep_to_Vrep_list(face_struct& face, size_t *output) except -1
