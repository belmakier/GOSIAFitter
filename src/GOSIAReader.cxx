#include "Gosia.h"
#include "GOSIAReader.h"

GOSIAReader::GOSIAReader(Nucleus* nucl, const char *datafilename){
	gosiaData.clear();
	fNucleus = *nucl;
	ReadGOSIAFile(datafilename);
}

GOSIAReader::GOSIAReader(Nucleus* nucl, out_yields &gosia_yields){
	gosiaData.clear();
	fNucleus = *nucl;
	ReadGOSIAFile(gosia_yields);
}

void	GOSIAReader::ReadGOSIAFile(const char* datafilename){

	std::ifstream infile(datafilename);

	int nExpt = 0;

	int initial_state;
	int final_state;
	double iJ;
	double fJ;
	double counts;
	double uncertainty;	// Not required for GOSIA

	ExperimentData tmpExpt;

	bool	flag = false;
	int	counter = 0;
	int	exptCounter = 0;

	std::string line;
	while(std::getline(infile,line)){
		if(flag){
			if(counter > 0){	// First GOSIA line after NORMALIZED YIELD is whitespace
				std::istringstream ss(line);
				ss	>> initial_state;
				ss	>> final_state;
				ss	>> iJ;
				ss	>> fJ;
				ss	>> counts;
				ss	>> uncertainty;
				//std::cout.precision(60); 
				//std::cout	<< line << "\n"
				//		<< counts << std::endl;
				initial_state--;	// GOSIA to Cygnus numbering
				final_state--;		// GOSIA to Cygnus numbering
	
				if(initial_state < fNucleus.GetNstates() && final_state < fNucleus.GetNstates() && initial_state >= 0 && final_state >= 0){
					tmpExpt.AddData(initial_state,final_state,counts,uncertainty);
				}
				if(initial_state == 1 && final_state == 0){
					flag = false;
					nExpt++;
					gosiaData.push_back(tmpExpt);
					continue;
				}
			}
			counter++;
		} 
		std::size_t found = line.find("NORMALIZED YIELD"); // Comment
		if(found != std::string::npos){
			flag = true;
			counter = 0;
			tmpExpt.ClearData();
		}
    std::size_t found_ruth = line.find("INTEGRATED RUTHERFORD"); // Comment
		if(found_ruth != std::string::npos){
      std::string junk;

      std::istringstream ss(line);
      ss >> junk;
      ss >> junk;
      ss >> junk;
      ss >> junk;
      rutherfords.push_back(std::atof(junk.substr(8).c_str()));      
		}

    std::size_t found_angle = line.find("ENERGY RANGE"); // Comment
		if(found_angle != std::string::npos){
      energy_low.push_back(std::atof(line.substr(19,7).c_str()));
      energy_high.push_back(std::atof(line.substr(30,7).c_str()));
      double angle_lo = std::atof(line.substr(67,7).c_str());
      double angle_hi = std::atof(line.substr(77,7).c_str());
      angle_low.push_back(angle_lo);
      angle_high.push_back(angle_hi);
      angle_av.push_back((angle_lo+angle_hi)/2.0);
		}
	}	
}

void	GOSIAReader::ReadGOSIAFile(out_yields &gosia_yields){
  //std::cout << "Reading GOSIA output" << std::endl;
  for (int i=0; i<gosia_yields.nexp; ++i) {
    int nyields = gosia_yields.experiment[i].nyields;
    //std::cout << "Exp: " << i << "  " << nyields << " yields" << std::endl;
    ExperimentData tmpExpt;
    for (int j=0; j<nyields; ++j) {
      int ni = gosia_yields.experiment[i].ni[j];
      int nf = gosia_yields.experiment[i].nf[j];
      double yld = gosia_yields.experiment[i].yield[j];
      ni--;	// GOSIA to Cygnus numbering
      nf--;		// GOSIA to Cygnus numbering
      if(ni < fNucleus.GetNstates() && nf < fNucleus.GetNstates() && ni >= 0 && nf >= 0){
        tmpExpt.AddData(ni,nf,yld,0);        
      }
    }
    gosiaData.push_back(tmpExpt);

    energy_low.push_back(gosia_yields.experiment[i].en_low);
    energy_high.push_back(gosia_yields.experiment[i].en_high);
    angle_low.push_back(gosia_yields.experiment[i].theta_low);
    angle_high.push_back(gosia_yields.experiment[i].theta_high);
    angle_av.push_back((gosia_yields.experiment[i].theta_high+
                        gosia_yields.experiment[i].theta_low)/2.);
    rutherfords.push_back(gosia_yields.experiment[i].ruth);      
  }
}

void	GOSIAReader::Print() const{

	for(size_t i=0;i<gosiaData.size();i++)
		gosiaData.at(i).Print();

}
