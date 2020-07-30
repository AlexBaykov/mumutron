double sigma_ang(double sigma_ecm){
	double a = (e*e - me*me)*sin(2*ang)/sqrt(2*me*me + 2*(e*e + (e*e - me*me)*cos(2*ang)));
	double e1 = e*(1 + cos(2*ang))/sqrt(2*me*me + 2*(e*e + (e*e - me*me)*cos(2*ang)));
	return sqrt((sigma_ecm*sigma_ecm - 2*e1*e1*sigma_eb*sigma_eb)/(2*a*a));
}


TTree * getTree(double sigma_ecm, int event_number, double t){
	TTree *event_tree = new TTree(Form("s%.3f", sigma_ecm), Form("event tree for sigma_ecm = %.3f", sigma_ecm));
	event_tree->Branch("1 muon y", &cy1, "cy1/D");
	event_tree->Branch("1 muon z", &z1, "z1/D");
	event_tree->Branch("2 muon y", &y2, "y2/D");
	event_tree->Branch("2 muon z", &z2, "z2/D");
	event_tree->Branch("1 muon px", &p1x, "p1x/D");
	event_tree->Branch("1 muon py", &p1y, "p1y/D");
	event_tree->Branch("1 muon pz", &p1z, "p1z/D");
	event_tree->Branch("2 muon px", &p2x, "p2x/D");
	event_tree->Branch("2 muon py", &p2y, "p2y/D");
	event_tree->Branch("2 muon pz", &p2z, "p2z/D");
	event_tree->Branch("E_cm", &ecm, "ecm/D");
	event_tree->Branch("cs_theta", &cs_theta, "cs_theta/D");

	for(int i = 0; i < event_number; i++){
		getEvent(sigma_ang(sigma_ecm), event_tree, t); 
	}
	return event_tree;
}
