{
TF1 *func = new TF1("histfit", histfit, 0, 5, 2);
func->SetParameter(0, 0.36);
func->SetParameter(1, 1e4);

h->Fit(func, "LVER");

h->SetMarkerStyle(kFullCircle);
h->SetLineColor(1);
h->Draw("E1 P");
}
