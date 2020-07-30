{
TFile *f = TFile::Open(Form("s0.300-0.600_n100.root"));
TFile *f1 = TFile::Open(Form("s0.396-0.396_n10.root"));

TH1D * hists[16];
TH1D *h = (TH1D *)f1->FindObjectAny("0.396");
double sigma[16];

for(int i = 0; i < 16; i++){
	hists[i] = (TH1D *)f->FindObjectAny(Form("%.3f", 0.3 + i*0.02));
	sigma[i] = atof(hists[i]->GetName());
	hists[i]->Scale(1./hists[i]->Integral());
}
}
