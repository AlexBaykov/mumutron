#include <boost/algorithm/string.hpp>
#include <fstream>
#include <vector>

std::vector<std::string> read(){
	std::ifstream sigma("sigma_pol.dat");
	std::vector<std::string> results;
	for(int i = 0; i < 1000; i++){
		results.push_back("0"); 
		getline(sigma, results[i]);
	}
	return results;
}

double getCS(double E, double mu){
	std::vector<std::string> lines;
	std::vector<std::string> words;
	double energy[1000];
	double visible[1000];
	
	lines = read();
	for(int i = 0; i < 1000; i++){
		boost::split(words, lines[i], [](char c){return c == ' ';});
		energy[i] = 2*mu + std::stod(words[0]);
		visible[i] = std::stod(words[2]); 
	}

	for(int i = 0; i < 999; i++){
		if(E < energy[i] && E < energy[i+1]){
			return visible[i-1]*(energy[i] - E)/(energy[i] - energy[i-1]) + visible[i]*(E - energy[i-1])/(energy[i] - energy[i-1]); 
		}
	}
	return 1;
}

