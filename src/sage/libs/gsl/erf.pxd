# distutils: libraries = GSL_LIBRARIES
from .types cimport *

cdef extern from "gsl/gsl_sf_erf.h":

  double  gsl_sf_erf(double x)

  int  gsl_sf_erf_e(double x, gsl_sf_result * result)

  double  gsl_sf_erfc(double x)

  int  gsl_sf_erfc_e(double x, gsl_sf_result * result)

  double  gsl_sf_log_erfc(double x)

  int  gsl_sf_log_erfc_e(double x, gsl_sf_result * result)

  double  gsl_sf_erf_Z(double x)

  int  gsl_sf_erf_Z_e(double x, gsl_sf_result * result)

  double  gsl_sf_erf_Q(double x)

  int  gsl_sf_erf_Q_e(double x, gsl_sf_result * result)

  double  gsl_sf_hazard(double x)

  int  gsl_sf_hazard_e(double x, gsl_sf_result * result)

