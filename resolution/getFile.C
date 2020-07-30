void getFile(double sigma1, double sigma2, double step, double t, int events_number){
	TFile f(Form("s%.3f-%.3f_n%d.root", sigma1, sigma2, events_number/1000), "recreate");
	TTree * data_tree = new TTree();
	TH1D *hist = new TH1D("name", "title", 25, 0, 5);
	double sigma = sigma1;
	int n = (sigma2 - sigma1)/step;

	for(int i = 0; i < n + 1; i++){
		std::cout << i << std::endl;
		data_tree = getTree(sigma, events_number, t);
		hist->SetName(Form("%.3f", sigma));
		hist->SetTitle(Form("hist for sigma_ecm = %.3f", sigma));
		data_tree->Draw(Form("sqrt(y2*y2 + z2*z2)>>%.3f", sigma));
		hist->Write();
		data_tree->Write();
		sigma += step;
	}
	f.Close();
}
