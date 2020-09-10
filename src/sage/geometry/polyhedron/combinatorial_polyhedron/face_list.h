#ifndef COMBINATORIAL_FACE_LIST_STRUCTURE
#define COMBINATORIAL_FACE_LIST_STRUCTURE

struct face_list_struct{
    face_struct* faces;
    size_t n_faces;
    size_t total_n_faces;
    size_t n_atoms;
    int polyhedron_is_simple;
};

void initialize_faces(face_list_struct& faces, size_t n_faces, size_t n_atoms);

void faces_copy(face_list_struct& dst, face_list_struct& src);

void faces_shallow_copy(face_list_struct& dst, face_list_struct& src);

void intersection(face_list_struct& dest, face_list_struct& A, face_struct& B);

size_t get_next_level(\
        face_list_struct& faces, \
        face_list_struct& new_faces, \
        face_list_struct& visited_all);

size_t bit_rep_to_coatom_rep(face_struct& face, face_list_struct& coatoms, size_t *output);

int test_alignment(face_list_struct& faces);

#endif
