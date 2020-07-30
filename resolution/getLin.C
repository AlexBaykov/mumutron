TH1D * getLin(double sigma_int, int norm){
	TH1D *lin = new TH1D("", "", 100, 0, 5);
	for(int bin = 1; bin < 101; bin++){
		for(int i = 0; i < 60; i ++){
			if(sigma_int >= sigma[i] && sigma_int <= sigma[i+1]){ 
				double y1, y2;
				y1 = hists[i]->GetBinContent(bin);
				y2 = hists[i+1]->GetBinContent(bin);
				lin->SetBinContent(bin, norm*(y1*(sigma[i+1] - sigma_int)/(sigma[i+1] - sigma[i]) + y2*(sigma_int - sigma[i])/(sigma[i+1] - sigma[i])));
			}
		}
	}
	return lin;
}

