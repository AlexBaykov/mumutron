#include "section.h"
#include "getCorrections.h"
#include "getExperimental.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <gsl/gsl_interp.h>

int main(){
	std::ofstream file;
	file.open("sigma.dat");
	double mu = 105.6583745;
	double E[30000];
	double visibleCS[30000];
	double bornCS[30000];
	double expCS[15500];
	double sigma = 0.396;

	double E_th = 2*mu;
	for(int i = 0; i < 30000; i++){
		E[i] = E_th - 10*sigma +  1e-3*i*sigma; 
		visibleCS[i] = getCorrections(E[i], E_th)*0.3894e12; 
		bornCS[i] = born(E[i]*E[i], E_th*E_th)*0.3894e12; 
	}

	gsl_interp * interp = gsl_interp_alloc(gsl_interp_linear, 30000);
	gsl_interp_init(interp, E, visibleCS, 30000);
	gsl_interp_accel *accel = gsl_interp_accel_alloc();

	for(int i = 0; i < 15500; i++){
		std::cout << E[i] << std::endl;
		expCS[i] = getExperimental(E[i + 6000], E_th, sigma, visibleCS, interp, accel, E);
		file << std::setprecision(9) <<  E[i] - E_th  << " " << expCS[i] << std::endl;
	}
	file.close();
	return 0;
}

