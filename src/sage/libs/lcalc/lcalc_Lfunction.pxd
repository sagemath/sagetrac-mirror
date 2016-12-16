cdef extern from "sage/libs/lcalc/lcalc_sage.h":

    ctypedef struct doublevec "std::vector<double>":
        int (*size)()
        double ind "operator[]"(int i)
        void (* clear)()

    doublevec doublevec_factory "std::vector<double>"(int len)

    cdef void initialize_globals()

    ctypedef struct c_Complex "Complex":
        double real()
        double imag()

    #######################
    #L function with (I)nteger Coefficients
    ######################

    ctypedef struct c_Lfunction_I "L_function<int>":
        c_Complex (* value) (c_Complex s, int derivative, char *whattype)
        int (* compute_rank) ()
        double (* N) (double T)
        void  (* find_zeros_v)(double T1, double T2, double stepsize, doublevec result )
        void (*find_zeros_via_N_v)(long count,int do_negative,double max_refine, int rank, int test_explicit_formula, doublevec result)
        void (*print_data_L)()

        #Constructor and destructor
    c_Lfunction_I *new_c_Lfunction_I "new L_function<int>"(char *NAME, int what_type, int N, int *coeff, long long Period, double q,  c_Complex w, int A, double *g, c_Complex *l, int n_poles, c_Complex *p, c_Complex *r)
    cdef void del_c_Lfunction_I "delete"(c_Lfunction_I *L)

    ################################
    #L function with (D)ouble Coefficients
    ################################

    ctypedef struct c_Lfunction_D "L_function<double>":
        c_Complex (* value) (c_Complex s, int derivative, char *whattype)
        int (* compute_rank) ()
        double (* N) (double T)
        double *dirichlet_coefficient
        void  (* find_zeros_v)(double T1, double T2, double stepsize, doublevec result )
        void (*find_zeros_via_N_v)(long count,int do_negative,double max_refine, int rank, int test_explicit_formula, doublevec result)
        void (*print_data_L)()

        #Constructor and destructor
    c_Lfunction_D *new_c_Lfunction_D "new L_function<double>"(char *NAME, int what_type, int N, double *coeff, long long Period, double q,  c_Complex w, int A, double *g, c_Complex *l, int n_poles, c_Complex *p, c_Complex *r)
    cdef void del_c_Lfunction_D "delete"(c_Lfunction_D *L)

    #######################
    #L function with (C)omplex Coefficients
    ######################

    ctypedef struct c_Lfunction_C "L_function<Complex>":
        c_Complex (* value) (c_Complex s, int derivative, char *whattype)
        int (* compute_rank) ()
        double (* N) (double T)
        void  (* find_zeros_v)(double T1, double T2, double stepsize, doublevec result )
        void (*find_zeros_via_N_v)(long count,int do_negative,double max_refine, int rank, int test_explicit_formula, doublevec result)
        void (*print_data_L)()

        #Constructor and destructor
    c_Lfunction_C *new_c_Lfunction_C "new L_function<Complex>"(char *NAME, int what_type, int N, c_Complex *coeff, long long Period, double q,  c_Complex w, int A, double *g, c_Complex *l, int n_poles, c_Complex *p, c_Complex *r)
    cdef void del_c_Lfunction_C "delete"(c_Lfunction_C *L)

    #######################
    #Zeta function
    ######################

    ctypedef struct c_Lfunction_Zeta "L_function<int>":
        c_Complex (* value) (c_Complex s, int derivative, char *whattype)
        int (* compute_rank) ()
        double (* N) (double T)
        void  (* find_zeros_v)(double T1, double T2, double stepsize, doublevec result )
        void (*find_zeros_via_N_v)(long count,int do_negative,double max_refine, int rank, int test_explicit_formula, doublevec result)#puts result in vector<double> result
        void (*find_zeros_via_N)(long count,int do_negative,double max_refine, int rank, int test_explicit_formula, char *filename) #puts result in filename

        #Constructor and destructor
    c_Lfunction_Zeta *new_c_Lfunction_Zeta "new L_function<int>"()
    cdef void del_c_Lfunction_Zeta "delete"(c_Lfunction_Zeta *L)

    #######################
    # Below are helper functions
    ######################

    cdef int *new_ints(int l)
    cdef void del_ints(int *)
    cdef double *new_doubles(int l)
    cdef void del_doubles(double *)
    cdef c_Complex *new_Complexes(int l)
    cdef void del_Complexes(c_Complex *)
    cdef c_Complex new_Complex(double r, double i)
    void delete "delete "(void *ptr)
    cdef void testL(c_Lfunction_C  *ptr)

################
#
#Below are definition of Lfunction classes with
# (I)nteger  (D)ouble or (C)omplex coefficients
#
################

# strange bug, I can't compile without this trick ???
# it's only used in __typedN
ctypedef double Double

cdef class Lfunction:
    cdef void *thisptr
    cdef void __init_fun(self, char *NAME, int what_type, dirichlet_coeff, long long Period, double q,  c_Complex w, int A, double *g, c_Complex *l, int n_poles, c_Complex *p, c_Complex *r)
    cdef c_Complex __value(self,c_Complex s,int derivative)
    cdef c_Complex __hardy_z_function(self,c_Complex s)
    cdef int __compute_rank(self)
    #strange bug, replacing Double with double gives me a compile error
    cdef Double __typedN(self, double T)
    cdef void __find_zeros_v(self, double T1, double T2, double stepsize,doublevec *result)
    cdef void __find_zeros_via_N_v(self, long count,int do_negative,double max_refine, int rank, int test_explicit_formula,doublevec *result)

    cdef str _repr

cdef class Lfunction_I(Lfunction):
    pass

cdef class Lfunction_D(Lfunction):
    pass

cdef class Lfunction_C(Lfunction):
    pass

cdef class Lfunction_Zeta(Lfunction):
    pass
