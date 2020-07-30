#include <functional>
#include "getCorrections.h"
#include <gsl/gsl_interp.h>

double getExperimental(double E, double E_th, double sigma, double* visibleCS, gsl_interp *interp, gsl_interp_accel * accel, double *En );
