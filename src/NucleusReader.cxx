#include "NucleusReader.h"

NucleusReader::NucleusReader() : fNucleus(NULL)
{
}

NucleusReader::NucleusReader(const char *filename, bool GOSIA) : fNucleus(NULL)
{
	if(GOSIA)
		ReadGOSIANucleus(filename);
	else
		ReadNucleusFile(filename);	
}

NucleusReader::NucleusReader(const NucleusReader &n) : fNucleus(n.fNucleus)
{
}

NucleusReader& NucleusReader::operator = (const NucleusReader &n)
{

	fNucleus	= n.fNucleus;

	return *this;

}

void NucleusReader::ReadNucleusFile(const char *filename){
	
	std::ifstream infile(filename);
	
	int counter = 0;
	int tmpA, tmpZ, tmpS, tmpL;
	int tmpI, tmpP; 
	double tmpE, tmpJ;

	int tmpI1, tmpI2, tmp_Lambda;
	std::string tmp_Lambda_s;
	double tmp_ME;

	std::string line;
	while(std::getline(infile,line)){
    if (line.size() == 0) { continue; }  //skip empty lines
		std::size_t found = line.find("!");
		if(found == std::string::npos){
			std::istringstream ss(line);

			if(counter == 0){
				ss >> tmpA >> tmpZ >> tmpS >> tmpL;
				fNucleus = new Nucleus(tmpA,tmpZ,tmpS,tmpL);
        present.resize(tmpL);
        for (int l=0; l<tmpL; ++l) {          
          present.at(l).ResizeTo(tmpS, tmpS);
          for (int i=0; i<tmpS; ++i) {
            for (int f=0; f<tmpS; ++f) {
              present.at(l)[i][f] = 0.0;
            }
          }
        }
			}
			else if(counter <= tmpS){
				ss >> tmpI >> tmpE >> tmpJ >> tmpP;
				fNucleus->SetState(tmpI,tmpE,tmpJ,tmpP);	
			}
			else{
				std::size_t foundE = line.find("E");
				std::size_t foundM = line.find("M");
				if(foundE != std::string::npos || foundM != std::string::npos){ // Found a transition specified E or M
					ss >> tmpI1 >> tmpI2 >> tmp_ME >> tmp_Lambda_s;
					if(foundE != std::string::npos && foundM !=std::string::npos){
						std::cout << "Both E and M found, not possible. Skipping line!" << std::endl;
						continue;
					}
					tmp_Lambda = tmp_Lambda_s[1] - '0'; // Since we know the second digit is the multipolarity
					if(foundM != std::string::npos)
						tmp_Lambda += 6;	
				}
				else{
					ss >> tmpI1 >> tmpI2 >> tmp_ME >> tmp_Lambda;	
				}
        if (tmp_Lambda > 0) {
          fNucleus->SetMatrixElement(tmp_Lambda-1,tmpI1,tmpI2,tmp_ME);
          present.at(tmp_Lambda-1)[tmpI1][tmpI2] = 1.0;
          present.at(tmp_Lambda-1)[tmpI2][tmpI1] = 1.0;
        }
        else {
          fNucleus->AddRhoSquared(tmpI1,tmpI2,tmp_ME*tmp_ME*1000.0);
        }
			}

			counter++;
		}

	}

	//fNucleus->PrintNucleus();

}

void NucleusReader::ReadGOSIANucleus(const char* filename){

	std::ifstream infile(filename);

	std::string line;

	std::vector<int>	tI;
	std::vector<int>	tP;
	std::vector<double>	tJ;
	std::vector<double>	tE;

	int	nExpt;
	int	tmpZ;
	int	tmpA;

	bool	exptFlag = false;
	while(std::getline(infile,line)){
		if(exptFlag){
			std::istringstream ss(line);
			ss	>> nExpt >> tmpZ >> tmpA;	
			break;	// That's all we need from here
		}
		if(line == "EXPT")
			exptFlag = true;
	}

	infile.clear();
	infile.seekg(0);

	bool	levelFlag = false;
	while(std::getline(infile,line)){

		if(levelFlag && (line == "0,0,0,0" || line == "0\t0\t0\t0")){
			levelFlag = false;
		}

		if(levelFlag){
			std::istringstream ss(line);
			int	tmpI,tmpP;
			double	tmpJ,tmpE;
			ss	>> tmpI >> tmpP >> tmpJ >> tmpE;
			tI.push_back(tmpI-1);
			tP.push_back(tmpP);
			tJ.push_back(tmpJ);
			tE.push_back(tmpE);		
		}

		if(line == "LEVE"){
			levelFlag = true;
		}
	}

	fNucleus = new Nucleus(tmpZ,tmpA,tI.size());
	for(size_t i=0;i<tI.size();i++)
		fNucleus->SetState(tI.at(i),tE.at(i),tJ.at(i),tP.at(i));

	infile.clear();
	infile.seekg(0);

	bst_i.clear();
	bst_f.clear();
	bst_l.clear();

	int	tmpLambda = -1;
	bool	MEFlag = false;
	while(std::getline(infile,line)){
		if(MEFlag && (line == "0,0,0,0,0" || line == "0\t0\t0\t0\t0"))
			MEFlag = false;

		if(MEFlag && (line == "1,0,0,0,0" || line == "1\t0\t0\t0\t0")){
			tmpLambda = 0;
		}
		else if(MEFlag && (line == "2,0,0,0,0" || line == "2\t0\t0\t0\t0")){
			tmpLambda = 1;
		}
		else if(MEFlag && (line == "3,0,0,0,0" || line == "3\t0\t0\t0\t0")){
			tmpLambda = 2;
		}
		else if(MEFlag && (line == "4,0,0,0,0" || line == "4\t0\t0\t0\t0")){
			tmpLambda = 3;
		}
		else if(MEFlag && (line == "5,0,0,0,0" || line == "5\t0\t0\t0\t0")){
			tmpLambda = 4;
		}
		else if(MEFlag && (line == "6,0,0,0,0" || line == "6\t0\t0\t0\t0")){
			tmpLambda = 5;
		}
		else if(MEFlag && (line == "7,0,0,0,0" || line == "7\t0\t0\t0\t0")){
			tmpLambda = 6;
		}
		else if(MEFlag && (line == "8,0,0,0,0" || line == "8\t0\t0\t0\t0")){
			tmpLambda = 7;
		}
		else if(MEFlag){
			std::istringstream ss(line);
			int	tmpI,tmpF;
			float	tmpME,tmpLL,tmpUL;
			ss	>> tmpI >> tmpF >> tmpME >> tmpLL >> tmpUL;
            
			fNucleus->SetMatrixElement(tmpLambda,tmpI-1,std::abs(tmpF)-1,tmpME);
			bst_i.push_back(tmpI-1);
			bst_f.push_back(std::abs(tmpF)-1);
			bst_l.push_back(tmpLambda);
		}
			
		if(line == "ME")
			MEFlag = true;
	}
}

void NucleusReader::WriteNucleusFile(const char *nucleusfilename) {
  //prints data file in Cygnus format
  std::ofstream nucleusfile(nucleusfilename);
  nucleusfile << "! Nucleus file created by GOSIAFitter" << std::endl;
  nucleusfile << "! Nucleus: Z, A, NStates, NLambda (7 == Max)" << std::endl;
  nucleusfile << fNucleus->GetZ() << "    " << fNucleus->GetA() << "   " << fNucleus->GetNstates() << "    " << fNucleus->GetMaxLambda() << std::endl;

  nucleusfile << "! States: Index, Energy (MeV), J, Parity" << std::endl;
  for (int i=0; i<fNucleus->GetNstates(); ++i) {
    char str[80];
    sprintf(str, "%i    %7.6f   %3.1f   %i\n", i, fNucleus->GetLevelEnergies()[i], fNucleus->GetLevelJ()[i], fNucleus->GetLevelP()[i]);
    nucleusfile << str;
  }

  nucleusfile << "! Matrix elements: initial index, final index, matrix element, lambda (1-6 = E1-6, 7 = M1)" << std::endl;
  std::vector<TMatrixD> MatrixElements = fNucleus->GetMatrixElements();
  std::string mult[8] = {"E1","E2","E3","E4","E5","E6","M1","M2"};
  for (int l=0; l<fNucleus->GetMaxLambda(); ++l) {
    TMatrixD mes = MatrixElements.at(l);
    for (int i=0; i<fNucleus->GetNstates(); ++i) {
      for (int j=i; j<fNucleus->GetNstates(); ++j) {
        if (mes[i][j] != 0) {
          char str[80];
          sprintf(str, "%i    %i    %7.6f    %s\n", i, j, mes[i][j], mult[l].c_str());
          nucleusfile << str;
        }
      }
    }          
  }

  nucleusfile.close();
}

void NucleusReader::WriteGOSIANucleus(const char *gosiafilename) {
  //writes EXPT and ME blocks for GOSIA input
  //this will write matrix elements for matrix elements present in 
  //also assumes all are fixed (upper bound = lower bound = 1)
  std::ofstream gosiafile(gosiafilename);
  gosiafile << "LEVE" << std::endl;
  for (int i=0; i<fNucleus->GetNstates(); ++i) {
    char str[80];

    sprintf(str, "%i    %i    %3.1f   %7.6f\n", i+1, fNucleus->GetLevelP()[i], fNucleus->GetLevelJ()[i], fNucleus->GetLevelEnergies()[i]);
    gosiafile << str;
  }
  gosiafile << "0,0,0,0" << std::endl;
  gosiafile << "ME" << std::endl;
  for (int l=0; l<fNucleus->GetMaxLambda(); ++l) {
    TMatrixD mes = fNucleus->GetMatrixElements().at(l);
    std::stringstream ss;
    int valid = 0;
    ss << l+1 << ",0,0,0,0" << std::endl;
    for (int i=0; i<fNucleus->GetNstates(); ++i) {
      for (int j=i; j<fNucleus->GetNstates(); ++j) {
        if (present.at(l)[i][j] == 0.0) { continue; }
        //check parity
        int lambda = l+1;
        int abs_lambda = l+1;
        if(l >= 6) {// Magnetic - opposite parity conventions to electric
          lambda -= 5;
          abs_lambda -= 6;
        }
        int	dP_Lambda = TMath::Power(-1,lambda);
        int	dP = fNucleus->GetLevelP().at(i)/fNucleus->GetLevelP().at(j);
        if (dP == dP_Lambda) {
          //triangle inequality
          if (std::abs(fNucleus->GetLevelJ().at(i) - fNucleus->GetLevelJ().at(j)) <= abs_lambda &&
              abs_lambda <= std::abs(fNucleus->GetLevelJ().at(i) + fNucleus->GetLevelJ().at(j))) {
            char str[80];
            sprintf(str, "%i    %i    %7.6f    %i    %i\n", i+1, j+1, mes[i][j], 1, 1);
            ss << str;
            ++valid;
          }
        }
      }
    }
    if (valid > 0) { gosiafile << ss.rdbuf(); }
  }

  gosiafile << "0,0,0,0,0" << std::endl;
}
 
void NucleusReader::WriteBST(const char *bstfilename) {
  //writes (starting) matrix elements to file in "GOSIA" order, the same as
  //the block written in NucleusReader::WriteGOSIANucleus
  //this also sets bst_i, bst_f, bst_l
  bst_i.clear();
  bst_f.clear();
  bst_l.clear();
  std::ofstream bstfile(bstfilename);
  for (int l=0; l<fNucleus->GetMaxLambda(); ++l) {
    TMatrixD mes = fNucleus->GetMatrixElements().at(l);
    for (int i=0; i<fNucleus->GetNstates(); ++i) {
      for (int j=i; j<fNucleus->GetNstates(); ++j) {
        if (present.at(l)[i][j] == 0.0) { continue; }
        //check parity
        int lambda = l+1;
        int abs_lambda = l+1;
        if(l >= 6) {// Magnetic - opposite parity conventions to electric
          lambda -= 5;
          abs_lambda -= 6;
        }
        int	dP_Lambda = TMath::Power(-1,lambda);
        int	dP = fNucleus->GetLevelP().at(i)/fNucleus->GetLevelP().at(j);
        if (dP == dP_Lambda) {
          //triangle inequality
          if (std::abs(fNucleus->GetLevelJ().at(i) - fNucleus->GetLevelJ().at(j)) <= abs_lambda &&
              abs_lambda <= std::abs(fNucleus->GetLevelJ().at(i) + fNucleus->GetLevelJ().at(j))) {
            char str[80];
            sprintf(str, "%7.6f\n", mes[i][j]);
            bstfile << str;
            bst_i.push_back(i);
            bst_f.push_back(j);
            bst_l.push_back(l);
          }
        }
      }
    }
  }
}

void NucleusReader::SetMapping() {
  bst_i.clear();
  bst_f.clear();
  bst_l.clear();
  for (int l=0; l<fNucleus->GetMaxLambda(); ++l) {
    TMatrixD mes = fNucleus->GetMatrixElements().at(l);
    for (int i=0; i<fNucleus->GetNstates(); ++i) {
      for (int j=i; j<fNucleus->GetNstates(); ++j) {
        if (present.at(l)[i][j] == 0.0) { continue; }
        //check parity
        int lambda = l+1;
        int abs_lambda = l+1;
        if(l >= 6) {// Magnetic - opposite parity conventions to electric
          lambda -= 5;
          abs_lambda -= 6;
        }
        int	dP_Lambda = TMath::Power(-1,lambda);
        int	dP = fNucleus->GetLevelP().at(i)/fNucleus->GetLevelP().at(j);
        if (dP == dP_Lambda) {
          //triangle inequality
          if (std::abs(fNucleus->GetLevelJ().at(i) - fNucleus->GetLevelJ().at(j)) <= abs_lambda &&
              abs_lambda <= std::abs(fNucleus->GetLevelJ().at(i) + fNucleus->GetLevelJ().at(j))) {
            bst_i.push_back(i);
            bst_f.push_back(j);
            bst_l.push_back(l);
          }
        }
      }
    }
  }
}
  

