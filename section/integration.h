#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <math.h>
#include <functional>

double func(double x, void * p);
double sIntegral(std::function<double(double)>* f,
                 double a, double b, double* abserr);
double Integral(std::function<double(double)>* f,
                double a, double b, double* abserr);
