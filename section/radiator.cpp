#include "radiator.h"
#define me 0.5109989461 
#define a 7.2973525698e-3

// double L(double E){
//     return log(2*E/me); 
// }

// double delta2(double E){

//     double res1 = (9.0/8 - 2*k2)*L(E)*L(E) - (45.0/16 - 11.0/2*k2 - 3*k3)*L(E);
//     double res2 = 6.0/5*k2*k2 - 9.0/2*k3 - 6*k2*log(2) + 3.0/8*k2 + 57.0/12;
//     return res1 + res2;
// }

// double Delta(double E){
//     return 1 + a/M_PI*(3.0/2*L(E) + 1.0/3*M_PI*M_PI - 2) + (a/M_PI)*(a/M_PI)*delta2(E);
// }

// double beta(double E){
//     return 2*a/M_PI*(L(E) - 1);
// }

// // E --> E * E !!!
// double radiator(double x, double E){
//     //double res1 = Delta(E)*beta(E)*pow(x, beta(E) - 1) - beta(E)/2*(2-x);
//     //double res2 = beta(E)*beta(E)/8*((2-x)*(3*log(1-x) - 4*log(x)) - 4*log(1-x)/x - 6 + x);
//     return Delta(E)*beta(E)*pow(x, beta(E) - 1) - beta(E)*(2-x)/2 + beta(E)*beta(E)*((2-x)*(3*log(1-x) - 4*log(x)) - 4*log(1-x)/x - 6 + x)/8;
// }

double L (double s) {
  return log(s / me / me);
}

double beta (double s) {
  return (2. * a) / M_PI * (L (s) - 1);
}


double radiator (double x, double s) {
  double lnX = log(x);
  double mX = 1 - x;
  double lnMX = log(mX);
  double sM = s / me / me;
  double logSM = log(sM);
  double E = 0.5 * sqrt(s);
  double res = 0;
    
  double part1 = beta (s) * pow(x, beta(s) - 1) *
    (1 + a / M_PI * (M_PI * M_PI / 3 - 0.5) +
     0.75 * beta (s) - beta (s) * beta (s) / 24 *
     (L (s) / 3 + 2 * M_PI * M_PI - 9.25));

  double part2 = -beta (s) * (1 - 0.5 * x);

  double part3 = 0.125 * beta (s) * beta (s) *
    (-4 * (2 - x) * lnX - (1 + 3 * mX * mX) / x * lnMX -
     6 + x);

  res = part1 + part2 + part3;

  if (x > 2 * me / E) {
    
    double subsubpart1 = 2 * lnX + logSM - 5. / 3;
  
    double subpart1 = pow(x - 2 * me / E , beta (s) ) / 6 / x *
      pow(subsubpart1, 2) * (2 - 2 * x + x * x + beta (s) / 3 * subsubpart1);
    
    double subpart2 = 0.5 * L (s) * L (s) * (2. / 3 * (1 - mX * mX *mX) / mX +
						 (2 - x) * lnMX + 0.5 * x);
    res += (a / M_PI) * (a / M_PI) * (subpart1 + subpart2);
  }

  return res;
}
