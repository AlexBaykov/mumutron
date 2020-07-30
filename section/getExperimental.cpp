#include "getExperimental.h"

double getExperimental(double E, double E_th, double sigma, double* visibleCS, gsl_interp *interp, gsl_interp_accel * accel, double *En ){
    
	double result = 0.;
    	double abserr = 0.;
    
    	std::function<double(double)> f1 = [&En, &visibleCS, interp, accel, E, sigma] (double x){
      		/**
       	  	(x - E)*(x - E)/sigma/sigma --> 0.5 * (x - E)*(x - E)/sigma/sigma
       	  	according to gaussian distribution
       		**/
      		return gsl_interp_eval(interp, En, visibleCS, x, accel) *
      		exp(-0.5 * (x - E) * (x - E) / sigma / sigma) /
      		sqrt(2 * M_PI) / sigma;
    	};
    	result = Integral(&f1, E - 6*sigma, E + 6*sigma, &abserr);
    
    
    	return result;
}
