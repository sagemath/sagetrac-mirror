# This file automatically generated from some bmi header files

cdef extern from "bmi.h":
    ctypedef int bool
    ctypedef long ba0_int_p
    ctypedef long M_INT
    enum bmi_balsa_typeof_object : 
        bmi_balsa_function_object
        bmi_balsa_table_object
        bmi_balsa_list_object
        bmi_balsa_string_object
        bmi_balsa_bool_object
        bmi_balsa_integer_object
        bmi_balsa_error_object 
    struct bmi_balsa_object :
        bmi_balsa_typeof_object type
        void* value
        bool dynamic
        int  nbref
    ctypedef bmi_balsa_object* ALGEB
    struct bmi_balsa_list :
        ba0_int_p alloc
        ba0_int_p size
        ALGEB *tab
    struct bmi_balsa_table :
        ALGEB type
        ALGEB notation
        ALGEB ordering
        ALGEB equations
    ALGEB bmi_balsa_default_notation ()
    ALGEB bmi_balsa_new_ALGEB (
        bmi_balsa_typeof_object type, void* value)
    ALGEB bmi_balsa_new_string (char*)
    ALGEB bmi_balsa_new_error ()
    ALGEB bmi_balsa_new_function (ALGEB, ba0_int_p)
    void bmi_balsa_increment_nbref (ALGEB)
    void bmi_balsa_decrement_nbref (ALGEB)
    void bmi_balsa_clear_ALGEB (void*)
    void bmi_balsa_printf_ALGEB (ALGEB)
    ctypedef void* MKernelVector
    ctypedef bool M_BOOL
    struct struct_RTableSettings :
        int data_type
        int order
        bool read_only
        int num_dimensions
    ctypedef struct_RTableSettings RTableSettings
    void RTableGetDefaults (MKernelVector, RTableSettings*)
    void* RTableCreate (MKernelVector, RTableSettings*, ALGEB, M_INT*) 
    void* RTableDataBlock (MKernelVector, ALGEB)
    void* MapleAlloc (MKernelVector, long)
    void MapleGcAllow (MKernelVector, ALGEB)
    void MapleDispose (MKernelVector, ALGEB)
    void MapleGcProtect (MKernelVector, ALGEB)
    void MapleCheckInterrupt (MKernelVector)
    ctypedef void bmi_balsa_error_proc (char*, void*)
    void MaplePushErrorProc (MKernelVector, bmi_balsa_error_proc*, void*)
    void MaplePopErrorProc (MKernelVector)
    ALGEB ToMapleName (MKernelVector, char*, M_BOOL)
    ALGEB ToMapleInteger (MKernelVector, ba0_int_p)
    long MapleToInteger32 (MKernelVector, ALGEB)
    ALGEB EvalMapleProc (MKernelVector, ALGEB, int, ...)
    bool IsMapleString (MKernelVector, ALGEB)
    bool IsMapleTable (MKernelVector, ALGEB)
    bool IsMapleName (MKernelVector, ALGEB)
    bool MapleToM_BOOL (MKernelVector, ALGEB)
    char* MapleToString (MKernelVector, ALGEB)
    ALGEB MapleTableSelect (MKernelVector, ALGEB, ALGEB)
    void MapleRaiseError (MKernelVector, char*)
    long MapleNumArgs (MKernelVector, ALGEB)
    ALGEB MapleListAlloc (MKernelVector, M_INT)
    void MapleListAssign (MKernelVector, ALGEB, M_INT, ALGEB)
    ALGEB ToMapleBoolean (MKernelVector, long)
    ALGEB bmi_balsa_new_differential_ring (ALGEB)
    ALGEB bmi_balsa_new_regchain (ALGEB)
    struct bmi_balsa_object_string :
        bmi_balsa_typeof_object type
        char* value
        bool dynamic
        int  nbref
    ctypedef bmi_balsa_object_string* ALGEB_string
    struct bmi_balsa_object_table :
        bmi_balsa_typeof_object type
        bmi_balsa_table* value
        bool dynamic
        int  nbref
    ctypedef bmi_balsa_object_table* ALGEB_table
    struct bmi_balsa_object_list :
        bmi_balsa_typeof_object type
        bmi_balsa_list* value
        bool dynamic
        int  nbref
    ctypedef bmi_balsa_object_list* ALGEB_list
    struct bmi_balsa_listof_table :
        ba0_int_p alloc
        ba0_int_p size
        ALGEB_table *tab
    struct bmi_balsa_object_listof_table :
        bmi_balsa_typeof_object type
        bmi_balsa_listof_table* value
        bool dynamic
        int  nbref
    ctypedef bmi_balsa_object_listof_table* ALGEB_listof_table
    struct bmi_balsa_listof_string :
        ba0_int_p alloc
        ba0_int_p size
        ALGEB_string *tab
    struct bmi_balsa_object_listof_string :
        bmi_balsa_typeof_object type
        bmi_balsa_listof_string* value
        bool dynamic
        int  nbref
    ctypedef bmi_balsa_object_listof_string* ALGEB_listof_string
    bool bmi_sage_is_error (void*)
    char* bmi_sage_mesgerr (void*)
    ALGEB bmi_sage_differential_ring (char*, char*, char*, char*)
    ALGEB bmi_sage_pretend_regular_differential_chain (
                char*, ALGEB, char*, bool, 
                char*, char*, ba0_int_p, ba0_int_p)
    ALGEB_string bmi_sage_ranking (ALGEB)
    ALGEB_string bmi_sage_attributes (ALGEB)
    ALGEB_string bmi_sage_indets (
                char*, ALGEB, char*, char*, 
                char*, char*, ba0_int_p, ba0_int_p)
    ALGEB_string bmi_sage_differential_prem (
                char*, char*, ALGEB, char*, char*, ba0_int_p, ba0_int_p)
    ALGEB_string bmi_sage_differential_prem2 (
                char*, char*, char*, ALGEB, 
                char*, char*, ba0_int_p, ba0_int_p)
    ALGEB_string bmi_sage_coeffs (
                char*, char*, char*, char*, ALGEB,
                char*, char*, ba0_int_p, ba0_int_p)
    ALGEB_string bmi_sage_equations (
                ALGEB, bool, bool, char*, char*, ba0_int_p, ba0_int_p)
    ALGEB_string bmi_sage_equations_with_criterion_DR (
                char*, ALGEB, bool, bool, char*, 
                char*, char*, ba0_int_p, ba0_int_p)
    ALGEB_string bmi_sage_equations_with_criterion_RDC (
                ALGEB, bool, bool, char*, char*, char*, ba0_int_p, ba0_int_p)
    ALGEB_string bmi_sage_factor_derivative (
                char*, ALGEB, char*, char*, ba0_int_p, ba0_int_p)
    ALGEB_string bmi_sage_is_constant (
                char*, char*, ALGEB, char*, char*, ba0_int_p, ba0_int_p)
    ALGEB_string bmi_sage_leading_derivative (
                char*, ALGEB, bool, char*, char*, ba0_int_p, ba0_int_p)
    ALGEB_string bmi_sage_leading_rank (
                char*, ALGEB, bool, bool, char*, char*, ba0_int_p, ba0_int_p)
    ALGEB_string bmi_sage_leading_coefficient (
                char*, ALGEB, bool, char*, char*, char*, ba0_int_p, ba0_int_p)
    ALGEB_string bmi_sage_tail (
                char*, ALGEB, bool, char*, char*, char*, ba0_int_p, ba0_int_p)
    ALGEB_string bmi_sage_separant (
                char*, ALGEB, bool, char*, char*, char*, ba0_int_p, ba0_int_p)
    ALGEB_string bmi_sage_integrate (
                ALGEB, char*, char*, bool, char*, char*, ba0_int_p, ba0_int_p)
    ALGEB_string bmi_sage_sort (
                char*, char*, ALGEB, char*, char*, ba0_int_p, ba0_int_p)
    ALGEB_string bmi_sage_number_of_equations (ALGEB)
    ALGEB_string bmi_sage_preparation_equation (
                char*, ALGEB, char*, char*, ba0_int_p, char*,
                char*, char*, ba0_int_p, ba0_int_p)
    ALGEB_list bmi_sage_RosenfeldGroebner (
                char*, char*, char*, char*, 
                char*, ALGEB, char*, char*, bool, char*, 
                char*, ba0_int_p, ba0_int_p)
    ALGEB bmi_sage_pardi (
                ALGEB, char*, bool,
                char*, char*, ba0_int_p, ba0_int_p)
    ALGEB_string bmi_sage_normal_form (
                ALGEB, char*, char*, char*, ba0_int_p, ba0_int_p)
    ALGEB_string bmi_sage_differentiate (
                ALGEB, char*, bool, char*,
                char*, char*, ba0_int_p, ba0_int_p)
    ALGEB_listof_string bmi_sage_process_equations (
                char*, ALGEB, char*, char*, ba0_int_p, ba0_int_p)
    ALGEB_string bmi_sage_base_field_generators (
                char*, char*, ALGEB, char*, char*, ba0_int_p, ba0_int_p)
