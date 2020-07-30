#include <TF1.h>
#include <TH1F.h>
#include <TH1D.h>

/********************
 par[0] - sigma_ecm
 par[1] - norm
********************/

Double_t histfit(Double_t *x, Double_t *par){
	Double_t val;
	Double_t arg;
	arg = x[0];
	int bin;
	bin = arg/0.2 + 1;  
	if(bin > 25){
		std::cout << bin << std::endl;
		std::cout << arg << std::endl;
		bin = 25;
	}	

	for(int i = 0; i < 15; i ++){
		if(par[0] >= sigma[i] && par[0] < sigma[i+1]){ 
			double y1, y2;
			y1 = hists[i]->GetBinContent(bin);
			y2 = hists[i+1]->GetBinContent(bin);
			val = par[1]*(y1*(sigma[i+1] - par[0])/(sigma[i+1] - sigma[i]) + y2*(par[0] - sigma[i])/(sigma[i+1] - sigma[i]));
		}
	}
  	return val;
}
