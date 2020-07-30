#include "section.h"
#include <iostream>
#define mu 105.6583745
#define a 7.2973525698e-3
double beta1(double s, double s_th){
    	return sqrt(1 - s_th/s);
}

double eta(double s, double s_th){
	return M_PI*a/beta1(s, s_th);
}

double born(double s, double s_th){
  	if (s > s_th) {
    		return 2*M_PI*M_PI*pow(a,3)*(1 - beta1(s, s_th)*beta1(s, s_th)/3)/(1 - exp(-eta(s, s_th)))/s;
  	}
  	if(s == s_th){
	  	return 2*M_PI*M_PI*pow(a, 3)/s;
  	}
  	return 0;
}
