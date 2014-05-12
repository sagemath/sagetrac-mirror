cdef extern from "cuba.h":
    ctypedef int (*integrand_t)(const int *ndim, const double x[],
        const int *ncomp, double f[], void *userdata)
    
    ctypedef void (*peakfinder_t)(const int *ndim, const double b[],
    int *n, double x[])
    
    void Vegas(const int ndim, const int ncomp,
    integrand_t integrand, void *userdata,
    const double epsrel, const double epsabs,
    const int flags, const int seed,
    const int mineval, const int maxeval,
    const int nstart, const int nincrease, const int nbatch,
    const int gridno, const char *statefile,
    int *neval, int *fail,
    double integral[], double error[], double prob[])
    
    void Suave(const int ndim, const int ncomp,
    integrand_t integrand, void *userdata,
    const double epsrel, const double epsabs,
    const int flags, const int seed,
    const int mineval, const int maxeval,
    const int nnew, const double flatness,
    const char *statefile,
    int *nregions, int *neval, int *fail,
    double integral[], double error[], double prob[])
    
    
    void Divonne(const int ndim, const int ncomp,
    integrand_t integrand, void *userdata,
    const double epsrel, const double epsabs,
    const int flags, const int seed,
    const int mineval, const int maxeval,
    const int key1, const int key2, const int key3, const int maxpass,
    const double border, const double maxchisq, const double mindeviation,
    const int ngiven, const int ldxgiven, double xgiven[],
    const int nextra, peakfinder_t peakfinder,
    const char *statefile,
    int *nregions, int *neval, int *fail,
    double integral[], double error[], double prob[])
    
    
    void Cuhre(const int ndim, const int ncomp,
    integrand_t integrand, void *userdata,
    const double epsrel, const double epsabs,
    const int flags, const int mineval, const int maxeval,
    const int key,
    const char *statefile,
    int *nregions, int *neval, int *fail,
    double integral[], double error[], double prob[])
    
    
    void cubasetinit(void (*)(), void *)
    void cubasetexit(void (*)(), void *)
    void cubaruninit()
    void cubaruninit()


