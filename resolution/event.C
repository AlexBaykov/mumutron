double me = 0.5109989461;
double alpha = 1/137.035999;
double m_mu = 105.6583745;
double sigma_a = 6.8e-4;
double ang = 75*TMath::Pi()/180;
double e = sqrt((2*m_mu*m_mu - me*me*(1 - cos(2*ang)))/(1 + cos(2*ang)));
double sigma_eb = 7.8e-4*e;
double m_sec = 600;
double thickness = 2;// Window thickness in mm
TLorentzVector p1;
TLorentzVector p2;
TVector3 vcm;
double cy1, z1, y2, z2, p1x, p1y, p1z, p2x, p2y, p2z, ecm, cs_theta, phi;

void getEnergy(){
	double E1 = gRandom->Gaus(e, sigma_eb);
	double E2 = gRandom->Gaus(e, sigma_eb);

	double a1 = gRandom->Gaus(ang, sigma_a);
	double a2 = gRandom->Gaus(ang, sigma_a);

	double alpha = a1 + a2;

	double p1 = sqrt(E1*E1 - me*me);
	double p2 = sqrt(E2*E2 - me*me);

	ecm = sqrt(2*me*me + 2*(E1*E2 + p1*p2*cos(alpha))); 

}

double section(double ecm, double beta, double cs_theta){
	double coef = 0;
	double eta = 0;
	if(beta < 1./1000000){
		coef = 1;
	}else{
		eta = TMath::Pi()*alpha/beta;
		coef = eta/(1 - exp(-eta)); 
	}
	double x = getCS(ecm, m_mu)*(1 + cs_theta*cs_theta + (1 - beta*beta)*(1 - cs_theta*cs_theta));
	return x; 
}

void getMultScat(TLorentzVector p){

	double beta = p.Rho()/p.E();
    	double phi = gRandom->Uniform(0, 1)*2*TMath::Pi();
    	double theta = gRandom->Gaus(0, 0.053*log(128.51)*thickness/(p.Rho()*p.Rho()*beta*beta));    
    	TVector3 v = p.Vect();
    	TVector3 v1 = v.Orthogonal();

    	v1.Rotate(phi, v);
    	p.Rotate(theta, v1);
}

void getEvent(double sigma_ang, TTree * event_tree, double t){	
	thickness = t;
	sigma_a = sigma_ang;
	ecm = 0;	
	while(ecm < 2*m_mu){
		getEnergy();
	}

	double emu = ecm/2;
	double pmu = sqrt(emu*emu - m_mu*m_mu);
	double beta = sqrt(1 - 4*m_mu*m_mu/ecm/ecm);
	double rand1 = gRandom->Uniform(0, 1);
	cs_theta = 2*gRandom->Uniform(0, 1) - 1;

	if(section(ecm, beta, cs_theta) < rand1*m_sec){
		getEvent(sigma_ang, event_tree, t);
	}else{
		phi = 2*TMath::Pi()*gRandom->Uniform(0, 1);
		double s_theta = sqrt(1 - cs_theta*cs_theta);
		vcm.SetXYZ(sin(75*TMath::Pi()/180), 0, 0);

		p1.SetXYZM(pmu*s_theta*cos(phi), pmu*s_theta*sin(phi), pmu*cs_theta, m_mu);
		p2.SetXYZM(-pmu*s_theta*cos(phi), -pmu*s_theta*sin(phi), -pmu*cs_theta, m_mu);
		p1.Boost(vcm);
		p2.Boost(vcm);

		cy1 = 10*p1.Y()/p1.X();
		z1 = 10*p1.Z()/p1.X();
		y2 = 10*p2.Y()/p2.X();
		z2 = 10*p2.Z()/p2.X();

		getMultScat(p1);
		getMultScat(p2);

		cy1 += 140*p1.Y()/p1.X() + gRandom->Gaus(0, 0.01);
		z1 += 140*p1.Z()/p1.X() + gRandom->Gaus(0, 0.01);
		y2 += 140*p2.Y()/p2.X() + gRandom->Gaus(0, 0.01);
		z2 += 140*p2.Z()/p2.X() + gRandom->Gaus(0, 0.01);

		p1x = p1.X();
		p1y = p1.Y();
		p1z = p1.Z();
		p2x = p2.X();
		p2y = p2.Y();
		p2z = p2.Z();

		event_tree->Fill();

	}
}
