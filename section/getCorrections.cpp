#include "getCorrections.h"
#include "integration.h"
#include <iostream>
#define me 0.5109989461 
#define mu 105.6583745 

double getCorrections(double E, double E_th) {
  	if (E < E_th) {
	    return 0.;
  	}
	double xm = 1.e-8;
	double xmax = 1 - 4. * me * me / E /E;
	std::function<double(double)> fcn = [E, E_th] (double x) {
		double tx = born(E * E * (1-x), E_th*E_th) * radiator(x, E * E);
      		return tx;
  	};
	double abserr = 0;
	double result;
	if (xm < xmax) {
		result = sIntegral(&fcn, 0, xm, &abserr);
		result += Integral(&fcn, xm, xmax, &abserr);
	} else {
		result = sIntegral(&fcn, 0, xmax, &abserr);
	}
	return result;
}
