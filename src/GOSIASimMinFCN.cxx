#include "GOSIASimMinFCN.h"
#include "Gosia.h"

void GOSIASimMinFCN::SetupCalculation(){

	exptIndex.resize(std::max(exptData_Beam.size(), exptData_Target.size()));
	for(unsigned int i=0;i<scalingParameters.size();i++) {
		for(unsigned int s=0;s<scalingParameters.at(i).GetExperimentNumbers().size();s++) {
			exptIndex[scalingParameters.at(i).GetExperimentNumbers().at(s)] = (int)i;
    }
  }

  if (verbosity>1) {
    std::cout	<< std::setw(13) << std::left << "Experiment: "
              << std::setw(14) << std::left << "Scaling index: "
              << std::endl;
    for(unsigned int i=0;i<exptIndex.size();i++){
      std::cout	<< std::setw(13) << std::left << i+1
                << std::setw(14) << std::left << exptIndex.at(i)
                << std::endl;
    }
  }

}

double GOSIASimMinFCN::CompareLifetimes(std::vector<LitLifetime> &litLifetimes, TransitionRates &rates, double &chisq_lifetime, int &NDF_lit, int &NDF ) {

  //rates_b.Print();

  double chisq = 0;
  int printct = 0;

  for(unsigned int i=0;i<litLifetimes.size();i++){
    double 	tmp = 0;
    int	index		= litLifetimes.at(i).GetIndex();
    double	lifetime	= litLifetimes.at(i).GetLifetime();
    double	calcLifetime	= rates.GetLifetimes()[index];
    if(fLikelihood){
      double	sigma		= litLifetimes.at(i).GetUpUnc() * litLifetimes.at(i).GetDnUnc();
      double	sigma_prime	= (litLifetimes.at(i).GetUpUnc() - litLifetimes.at(i).GetDnUnc());
      chisq 			+= 0.5 * TMath::Power((calcLifetime - lifetime),2)/(sigma + sigma_prime * (calcLifetime - lifetime));
    }
    else{
      if(calcLifetime > lifetime)
        tmp = (calcLifetime - lifetime) / litLifetimes.at(i).GetUpUnc();
      else
        tmp = (calcLifetime - lifetime) / litLifetimes.at(i).GetDnUnc();
      chisq += tmp * tmp;
      chisq_lifetime += tmp*tmp;
      if (verbosity>1) {
        std::cout << std::setw(5) << std::left << index 
                  << std::setw(12) << std::left << calcLifetime
                  << std::setw(12) << std::left << lifetime
                  << std::setw(12) << std::left << litLifetimes.at(i).GetUpUnc()
                  << std::setw(12) << std::left << tmp*tmp;
        if (printct%2 == 1) {
          std::cout << std::endl;
        }
        else {
          std::cout << std::setw(7) << " "
                    << std::setw(1) << "|"
                    << std::setw(7) << " ";
        }
        printct += 1;
      }
    }
    NDF++;
    NDF_lit++;
  }
  return chisq;
}

double GOSIASimMinFCN::CompareBranchingRatios(std::vector<LitBranchingRatio> &litBranchingRatios, TransitionRates &rates, double &chisq_branch, int &NDF_lit, int &NDF) {

  int printct = 0;
  double chisq = 0;
  for(unsigned int i=0;i<litBranchingRatios.size();i++){
    double 	tmp 		= 0;
    int	index_init 	= litBranchingRatios.at(i).GetInitialIndex();
    int	index_final1	= litBranchingRatios.at(i).GetFinalIndex_1();
    int	index_final2	= litBranchingRatios.at(i).GetFinalIndex_2();
    double	BR		= litBranchingRatios.at(i).GetBranchingRatio();
    double  calcBR		= rates.GetBranchingRatios()[index_final1][index_init] / rates.GetBranchingRatios()[index_final2][index_init];
    if(fLikelihood){
      double	sigma		= litBranchingRatios.at(i).GetUpUnc() * litBranchingRatios.at(i).GetDnUnc();
      double	sigma_prime	= (litBranchingRatios.at(i).GetUpUnc() - litBranchingRatios.at(i).GetDnUnc());
      chisq 			+= 0.5 * TMath::Power((calcBR - BR),2)/(sigma + sigma_prime * (calcBR - BR));
    }
    else{
      if(calcBR > BR) {
        tmp = (BR - calcBR) / litBranchingRatios.at(i).GetUpUnc();
      }
      else {
        tmp = (BR - calcBR) / litBranchingRatios.at(i).GetDnUnc();
      }
      chisq += tmp * tmp;
      chisq_branch += tmp*tmp;      
      if (verbosity>1) {
        std::cout	<< std::setw(5) << std::left << index_init 
                  << std::setw(5) << std::left << index_final1
                  << std::setw(5) << std::left << index_final2
                  << std::setw(13) << std::left << calcBR
                  << std::setw(10) << std::left << BR
                  << std::setw(10) << std::left << litBranchingRatios.at(i).GetUpUnc()
                  << std::setw(13) << std::left << tmp*tmp ;
        if (printct%2 == 1) {
          std::cout << std::endl;
        }
        else {
          std::cout << std::setw(7) << " "
                    << std::setw(1) << "|"
                    << std::setw(7) << " ";
        }
        printct += 1;
      }
    }
    NDF++;
    NDF_lit++;
  }
  return chisq;
}

double GOSIASimMinFCN::CompareMixingRatios(std::vector<LitMixingRatio> &litMixingRatios, TransitionRates &rates, double &chisq_mixingratio, int &NDF_lit, int &NDF) {
  int printct = 0;
  double chisq = 0;

  for(unsigned int i=0;i<litMixingRatios.size();i++){
    double 	tmp = 0;
    int 	index_init	= litMixingRatios.at(i).GetInitialIndex();
    int	index_final	= litMixingRatios.at(i).GetFinalIndex();
    double	delta		= litMixingRatios.at(i).GetMixingRatio();
    double	calcDelta	= rates.GetMixingRatios()[index_final][index_init];
    if(fLikelihood){
      double	sigma		= litMixingRatios.at(i).GetUpUnc() * litMixingRatios.at(i).GetDnUnc();
      double	sigma_prime	= (litMixingRatios.at(i).GetUpUnc() - litMixingRatios.at(i).GetDnUnc());
      chisq 			+= 0.5 * TMath::Power((calcDelta - delta),2)/(sigma + sigma_prime * (calcDelta - delta));
    }
    else{
      if(calcDelta > delta)
        tmp = (delta - calcDelta) / litMixingRatios.at(i).GetUpUnc();
      else
        tmp = (delta - calcDelta) / litMixingRatios.at(i).GetDnUnc();
      chisq += tmp * tmp;		
      chisq_mixingratio += tmp*tmp;
      if (verbosity>1) {
        std::cout	<< std::setw(5) << std::left << index_init 
                  << std::setw(5) << std::left << index_final
                  << std::setw(12) << std::left << calcDelta
                  << std::setw(10) << std::left << delta
                  << std::setw(10) << std::left << litMixingRatios.at(i).GetUpUnc()
                  << std::setw(12) << std::left << tmp*tmp;
        if (printct%2 == 1) {
          std::cout << std::endl;
        }
        else {
          std::cout << std::setw(7) << " "
                    << std::setw(1) << "|"
                    << std::setw(7) << " ";
        }
        printct += 1;

      }
    }
    NDF++;
    NDF_lit++;
  }
  return chisq;
}

double GOSIASimMinFCN::CompareMatrixElements(std::vector<LitMatrixElement> &litMatrixElements, Nucleus &nucl, double &chisq_matel, int &NDF_lit, int &NDF) {
  double chisq = 0;
  int printct = 0;
  for(unsigned int i=0;i<litMatrixElements.size();i++){
    double tmp = 0;
    int	mult		= litMatrixElements.at(i).GetMultipolarity();
    int 	index_init	= litMatrixElements.at(i).GetInitialIndex();
    int	index_final	= litMatrixElements.at(i).GetFinalIndex();
    double	ME		= litMatrixElements.at(i).GetMatrixElement();
    double	calcME		= nucl.GetMatrixElements().at(mult)[index_init][index_final];
    if(fLikelihood){
      double	sigma		= litMatrixElements.at(i).GetUpUnc() * litMatrixElements.at(i).GetDnUnc();
      double	sigma_prime	= (litMatrixElements.at(i).GetUpUnc() - litMatrixElements.at(i).GetDnUnc());
      chisq 			+= 0.5 * TMath::Power((calcME - ME),2)/(sigma + sigma_prime * (calcME - ME));
    }
    else{
      if (std::abs(litMatrixElements.at(i).GetSign()) > 0) {
        if(calcME > ME) {
          tmp = (ME - calcME) / litMatrixElements.at(i).GetUpUnc();
        }
        else {
          tmp = (ME - calcME) / litMatrixElements.at(i).GetDnUnc();
        }
      }
      else { //don't know sign of literature matrix element, compare absolute values
        if(std::abs(calcME) > std::abs(ME)) {
          tmp = (std::abs(ME) - std::abs(calcME)) / litMatrixElements.at(i).GetUpUnc();
        }
        else {
          tmp = (std::abs(ME) - std::abs(calcME)) / litMatrixElements.at(i).GetDnUnc();
        }
      }
        
      chisq += tmp * tmp;		
      chisq_matel += tmp*tmp;
      if (verbosity>1) {
        std::cout	<< std::setw(5) << std::left << index_init 
                  << std::setw(5) << std::left << index_final
                  << std::setw(5) << std::left << mult
                  << std::setw(14) << std::left << calcME
                  << std::setw(10) << std::left << ME
                  << std::setw(10) << std::left << litMatrixElements.at(i).GetUpUnc()
                  << std::setw(14) << std::left << tmp*tmp;
        if (printct%2 == 1) {
          std::cout << std::endl;
        }
        else {
          std::cout << std::setw(7) << " "
                    << std::setw(1) << "|"
                    << std::setw(7) << " ";
        }
        printct += 1;

      }
    }
    NDF++;
    NDF_lit++;
  }
  return chisq;
}

void GOSIASimMinFCN::PrintYieldHeader(int species) {
  if (species == 0) { //beam
    if (exptData_Beam.size() == 0) { return; }
    std::cout 	<< std::setw( 7) << std::left << "Beam:"
                << std::endl;
  }
  else if (species == 1) { //target
    if (exptData_Target.size() == 0) { return; }
    std::cout 	<< std::setw( 7) << std::left << "Target:"
                << std::endl;
  }
  std::cout 	<< std::setw( 4) << std::left << " ";
  std::cout	<< std::setw(10) << std::left << "Scaling:";

  if (!fLikelihood) {
    std::cout 	<< std::setw( 6) << std::left << "Init:"
                << std::setw( 6) << std::left << "Finl:"
                << std::setw(14) << std::left << "Calc:"
                << std::setw(14) << std::left << "Expt:"
                << std::setw(14) << std::left << "Err:"
                << std::setw(14) << std::left << "C/E:"
                << std::setw(14) << std::left << "Chisq:"
                << std::setw(20) << std::left << " "
                << std::setw( 6) << std::left << "Init:"
                << std::setw( 6) << std::left << "Finl:"
                << std::setw(14) << std::left << "Calc:"
                << std::setw(14) << std::left << "Expt:"
                << std::setw(14) << std::left << "Err:"
                << std::setw(14) << std::left << "C/E:"
                << std::setw(14) << std::left << "Chisq:";
  }
  else {
       std::cout 	<< std::setw( 6) << std::left << "Init:"
                << std::setw( 6) << std::left << "Finl:"
                << std::setw(14) << std::left << "Calc:"
                << std::setw(14) << std::left << "Expt:"
                << std::setw(14) << std::left << "C/E:"
                << std::setw(14) << std::left << "-Ln(L) cont.:"
                << std::setw(20) << std::left << " "
                << std::setw( 6) << std::left << "Init:"
                << std::setw( 6) << std::left << "Finl:"
                << std::setw(14) << std::left << "Calc:"
                << std::setw(14) << std::left << "Expt:"
                << std::setw(14) << std::left << "C/E:"
                << std::setw(14) << std::left << "-Ln(L) cont.:";
  } 
  std::cout 	<< std::endl;
}
  
double GOSIASimMinFCN::CompareYields(std::vector<ExperimentData> &exptData,
                                     std::vector<TMatrixD> &EffectiveCrossSection,
                                     std::vector<double> &scaling,
                                     std::vector<double> &exptchisq,
                                     double &chisqspecies,
                                     int &NDFspecies,
                                     int &NDF) {
  double chisq = 0;
  exptchisq.resize(exptData.size());
  for(unsigned int i=0;i<exptData.size();i++){
    exptchisq[i] = 0;
    if(expt_weights.at(i) == 0) 
      continue;
    if(verbosity>1){
      std::cout	<< std::setw(4) << std::left << i+1;
      std::cout	<< std::setw(10) << std::left << scaling.at(i);
    }
    int print_ct = 0;
    for(unsigned int t=0;t<exptData.at(i).GetData().size();++t){
      double 	tmp 		= 0;
      int	index_init 	= exptData.at(i).GetData().at(t).GetInitialIndex();
      int	index_final 	= exptData.at(i).GetData().at(t).GetFinalIndex();
      double 	calcCounts 	= scaling.at(i) * EffectiveCrossSection.at(i)[index_final][index_init];
      double 	exptCounts 	= exptData.at(i).GetData().at(t).GetCounts();
      double	sigma		= exptData.at(i).GetData().at(t).GetUpUnc() * exptData.at(i).GetData().at(t).GetDnUnc();
      double	sigma_prime	= (exptData.at(i).GetData().at(t).GetUpUnc() - exptData.at(i).GetData().at(t).GetDnUnc());
      sigma			/= expt_weights.at(i);
      sigma_prime		/= expt_weights.at(i);

      if (exptCounts == 0) { //treat as limit
        double limit = exptData.at(i).GetData().at(t).GetUpUnc();
        if (calcCounts > limit) {
          exptCounts = 0.01;
        }
      }

      if(exptCounts > 0 && sigma > 0){
        if (verbosity>1) {
          if(fLikelihood){
            std::cout 	<< std::setw( 6) << std::left << index_init 
                        << std::setw( 6) << std::left << index_final 
                        << std::setw(14) << std::left << calcCounts 
                        << std::setw(14) << std::left << exptCounts 
                        << std::setw(14) << std::left << calcCounts/exptCounts
                        << std::setw(14) << std::left << 0.5 * TMath::Power((exptCounts - calcCounts),2)/(sigma + sigma_prime * (exptCounts - calcCounts)) << std::endl << "              ";
          }
          else{
            std::cout 	<< std::setw( 6) << std::left << index_init 
                        << std::setw( 6) << std::left << index_final 
                        << std::setw(14) << std::left << calcCounts 
                        << std::setw(14) << std::left << exptCounts 
                        << std::setw(14) << std::left << exptData.at(i).GetData().at(t).GetUpUnc()
                        << std::setw(14) << std::left << calcCounts/exptCounts
                        << std::setw(14) << std::left << TMath::Power((calcCounts - exptCounts)/exptData.at(i).GetData().at(t).GetUpUnc(),2);              
            if (print_ct%2 == 0) {
              std::cout << std::setw(20) << std::left << " ";
              if (t==(exptData.at(i).GetData().size()+
                      exptData.at(i).GetDoublet().size()-1)) { std::cout << std::endl; }
            }
            else { std::cout << std::endl << std::setw(14) << std::left << " "; }
          }
          print_ct += 1;
        }
        if(fLikelihood){
          chisq		+= 0.5 * TMath::Power((exptCounts - calcCounts),2)/(sigma + sigma_prime * (exptCounts - calcCounts));
          chisqspecies += 0.5 * TMath::Power((exptCounts - calcCounts),2)/(sigma + sigma_prime * (exptCounts - calcCounts));
        }
        else{
          if(calcCounts > exptCounts)
            tmp 		= (calcCounts - exptCounts) / exptData.at(i).GetData().at(t).GetUpUnc();
          else
            tmp 		= (calcCounts - exptCounts) / exptData.at(i).GetData().at(t).GetDnUnc();
          chisq		+= tmp * tmp;
          chisqspecies += tmp * tmp;
          exptchisq[i]	+= tmp * tmp;
        }
        NDF++;
        NDFspecies++;
      }
    }
    for(unsigned int t=0;t<exptData.at(i).GetDoublet().size();++t){
      double 	tmp 		= 0;
      int	index_init1 	= exptData.at(i).GetDoublet().at(t).GetInitialIndex1();
      int	index_final1 	= exptData.at(i).GetDoublet().at(t).GetFinalIndex1();
      int	index_init2 	= exptData.at(i).GetDoublet().at(t).GetInitialIndex2();
      int	index_final2 	= exptData.at(i).GetDoublet().at(t).GetFinalIndex2();
      double 	calcCounts 	= scaling.at(i) * (EffectiveCrossSection.at(i)[index_final1][index_init1] + EffectiveCrossSection.at(i)[index_final2][index_init2]);
      double 	exptCounts 	= exptData.at(i).GetDoublet().at(t).GetCounts();
      double	sigma		= exptData.at(i).GetDoublet().at(t).GetUpUnc() * exptData.at(i).GetDoublet().at(t).GetDnUnc();
      double	sigma_prime	= (exptData.at(i).GetDoublet().at(t).GetUpUnc() - exptData.at(i).GetDoublet().at(t).GetDnUnc());
      sigma			/= expt_weights.at(i);
      sigma_prime		/= expt_weights.at(i);
      if (exptCounts == 0) { //treat as limit
        double limit = exptData.at(i).GetDoublet().at(t).GetUpUnc();
        if (calcCounts > limit) {
          exptCounts = 0.01;
        }
      }
      else if(calcCounts > 0 && sigma > 0){
          if (verbosity>1) {
            if(fLikelihood){
              std::cout 	<< std::setw( 6) << std::left << index_init1*100 + index_init2
                          << std::setw( 6) << std::left << index_final1*100 + index_final2
                          << std::setw(14) << std::left << calcCounts 
                          << std::setw(14) << std::left << exptCounts 
                          << std::setw(14) << std::left << calcCounts/exptCounts
                          << std::setw(14) << std::left << 0.5 * TMath::Power((exptCounts - calcCounts),2)/(sigma + sigma_prime * (exptCounts - calcCounts));
            }
            else{
              std::cout 	<< std::setw( 6) << std::left << index_init1*100 + index_init2 
                          << std::setw( 6) << std::left << index_final1*100 + index_final2
                          << std::setw(14) << std::left << calcCounts 
                          << std::setw(14) << std::left << exptCounts 
                          << std::setw(14) << std::left << exptData.at(i).GetDoublet().at(t).GetUpUnc()
                          << std::setw(14) << std::left << calcCounts/exptCounts
                          << std::setw(14) << std::left << TMath::Power((calcCounts - exptCounts)/exptData.at(i).GetDoublet().at(t).GetUpUnc(),2);
              if (print_ct%2 == 0) {
                std::cout << std::setw(20) << std::left << " ";
                if (t==exptData.at(i).GetDoublet().size()-1) { std::cout << std::endl; }
              }
              else { std::cout << std::endl << std::setw(14) << std::left << " "; }
            }
            print_ct += 1;
          }
        if(fLikelihood){
          chisq		+= 0.5 * TMath::Power((exptCounts - calcCounts),2)/(sigma + sigma_prime * (exptCounts - calcCounts));
          chisqspecies	+= 0.5 * TMath::Power((exptCounts - calcCounts),2)/(sigma + sigma_prime * (exptCounts - calcCounts));
        }
        else{
          if(calcCounts > exptCounts)
            tmp 		= (calcCounts - exptCounts) / exptData.at(i).GetDoublet().at(t).GetUpUnc();
          else
            tmp 		= (calcCounts - exptCounts) / exptData.at(i).GetDoublet().at(t).GetDnUnc();
          chisq		+= tmp * tmp;
          chisqspecies	+= tmp * tmp;
          exptchisq[i]	+= tmp * tmp;
        }
        NDF++;
        NDFspecies++;
      }
    }
    if(verbosity>1)
      std::cout << std::endl;
  }
  return chisq;
}

double GOSIASimMinFCN::operator()(const double* par){
	std::cout 	<< std::setprecision(6);

  if (simanmin) {
    if (simanmin->GetRecorr()) {
      CalcBeamCorrectionFactors();
      CalcTargetCorrectionFactors();
    }
  }
  if (integralAlways) {
    CalcBeamCorrectionFactors();
    CalcTargetCorrectionFactors();
  }
  
	double chisq = 0;
	int 	NDF = 0;

	int	NDF_lit = 0;
	int	NDF_beam = 0;
	int 	NDF_targ = 0;

	typedef std::chrono::high_resolution_clock Clock;
	typedef std::chrono::milliseconds milliseconds;
	Clock::time_point t0 = Clock::now();

	Nucleus &nucl_b = fNucleus_Beam;	
	Nucleus &nucl_t = fNucleus_Target;

    
	parameters.clear();
  int parct = 0;
  int nRelBeam = 0;
  int nRelTarget = 0;
  int nBeamME = 0;
  int nTargetME = 0;

  for (unsigned int i=0; i<FE_Beam.size(); ++i) {
    FE_Beam.at(i)->Propagate(nucl_b,par,parct,0);  //this propagates the numbers from par into nucl_b correctly, as implemented in the FittingElement derived classes
  }

  for (unsigned int i=0; i<FE_Target.size(); ++i) {
    FE_Target.at(i)->Propagate(nucl_t,par,parct,0);  //this propagates the numbers from par into nucl_t correctly, as implemented in the FittingElement derived classes
  }
  
  if (verbosity>1) {
    std::cout << std::endl;
    std::cout << std::setw(4) << parct << " Parameters: ";
    int linect = 0;
    for (int i=0; i<FE_Beam.size(); ++i) {
      if (FE_Beam[i]->GetFixed()) { continue; }
      for (int j=0; j<FE_Beam[i]->GetNPars(); ++j) {
        std::cout << std::setw(10) << FE_Beam[i]->GetName()
                  << std::setw(1) << " "
                  << std::setw(7) << FE_Beam[i]->GetType()
                  << std::setw(2) << " "
                  << std::setw(12) << par[FE_Beam[i]->GetIndex()+j];
        if (linect < 4) {
          std::cout << "  |  ";
          linect += 1;
        }
        else {
          std::cout << std::endl << std::setw(17) << " ";
          linect = 0;
        }
      }      
    }
    for (int i=0; i<FE_Target.size(); ++i) {
      if (FE_Target[i]->GetFixed()) { continue; }
      for (int j=0; j<FE_Target[i]->GetNPars(); ++j) {
        std::cout << std::setw(10) << FE_Target[i]->GetName()
                  << std::setw(1) << " "
                  << std::setw(7) << FE_Target[i]->GetType()
                  << std::setw(2) << " "
                  << std::setw(12) << par[FE_Target[i]->GetIndex()+j];
        if (linect < 4) {
          std::cout << "  |  ";
          linect += 1;
        }
        else {
          std::cout << std::endl << std::setw(17) << " ";
          linect = 0;
        }
      }      
    }
    std::cout << std::endl;
  }

  // 	COMPARE WITH LITERATURE CONSTRAINTS:
  TransitionRates rates_b(&nucl_b);
  TransitionRates rates_t(&nucl_t);

  // 	First, compare with the literature for the beam:
  double lifetime_chisq = 0;
  double br_chisq = 0;
  double mr_chisq = 0;
  double me_chisq = 0;
  
  if (verbosity > 1 &&litLifetimes_Beam.size() > 0) {
    std::cout << "Lifetimes (Beam):" << std::endl;
  }
  chisq += CompareLifetimes(litLifetimes_Beam, rates_b, lifetime_chisq, NDF_lit, NDF);
  
  if(litBranchingRatios_Beam.size()>0 && verbosity>1) {
    std::cout	<< std::endl << "BR (Beam):" 
              << std::endl;
  }
  chisq += CompareBranchingRatios(litBranchingRatios_Beam, rates_b, br_chisq, NDF_lit, NDF);
  
  if(litMixingRatios_Beam.size()>0 && verbosity>1) {
    std::cout	<< std::endl << "Delta (Beam):" 
              << std::endl;
  }
  chisq += CompareMixingRatios(litMixingRatios_Beam, rates_b, mr_chisq, NDF_lit, NDF);

  if(litMatrixElements_Beam.size()>0 && verbosity>1) {
    std::cout 	<< std::endl << "Matrix Elements (Beam)"
                << std::endl;
  }
  chisq += CompareMatrixElements(litMatrixElements_Beam, nucl_b, me_chisq, NDF_lit, NDF);

  // 	Now, compare with the literature for the target:
  if(litLifetimes_Target.size()>0 && verbosity>1) {
    std::cout	<< std::endl << "Lifetimes (Target):" 
              << std::endl;
  }
  chisq += CompareLifetimes(litLifetimes_Target, rates_t, lifetime_chisq, NDF_lit, NDF);

  if(litBranchingRatios_Target.size()>0 && verbosity>1) {
    std::cout	<< std::endl << "BR (Target):" 
              << std::endl;
  }
  chisq += CompareBranchingRatios(litBranchingRatios_Target, rates_t, br_chisq, NDF_lit, NDF);

  if(litMixingRatios_Target.size()>0 && verbosity>1 ) {
    std::cout	<< std::endl << "Delta (Target):" 
              << std::endl;
  }
  chisq += CompareMixingRatios(litMixingRatios_Target, rates_t, mr_chisq, NDF_lit, NDF);
  
  if(litMatrixElements_Target.size()>0 && verbosity > 1) {
    std::cout 	<< std::endl << "Matrix Elements (Target)"
                << std::endl;
  }
  chisq += CompareMatrixElements(litMatrixElements_Target, nucl_t, me_chisq, NDF_lit, NDF);
  if (verbosity > 1) {
    std::cout << std::flush << std::endl;
  }

  double litchisq = chisq;
  //	COULEX AND STUFF:
  int verb = 0;
  Gosia(0,verb, 0);
  Gosia(1,verb, 0);

  //RunGosia(beam_inputfile, workingDir, all_detectors, beam_me, beam_out, verbosity);
  //RunGosia(target_inputfile, workingDir, all_detectors, target_me, target_out, verbosity);

  GOSIAReader beam_gosiaReader(&nucl_b,beam_yields);	//	Grab the GOSIA yields
  GOSIAReader target_gosiaReader(&nucl_t,target_yields);	//	Grab the GOSIA yields

  std::vector<ExperimentData>	beamCalc	= beam_gosiaReader.GetGOSIAData();
  std::vector<ExperimentData>	targetCalc	= target_gosiaReader.GetGOSIAData();
  EffectiveCrossSection_Beam.clear();	
  EffectiveCrossSection_Target.clear();	

  int dim_b = rates_b.GetBranchingRatios().GetNrows();
  int dim_t = rates_t.GetBranchingRatios().GetNrows();

  //do the correction from point to integral yields
  for(size_t i=0; i<beamCalc.size(); i++){
    EffectiveCrossSection_Beam.push_back(beamCalc.at(i).GetEffectiveCrossSection(correctionFactors_Beam.at(i), dim_b));
  };

  for(size_t i=0; i<targetCalc.size(); i++){
    EffectiveCrossSection_Target.push_back(targetCalc.at(i).GetEffectiveCrossSection(correctionFactors_Target.at(i), dim_t));
  }

  if(verbosity>1)
    std::cout << std::endl;

  //fit optimal scaling factors
  std::vector<double>	scaling;
  scaling.resize(std::max(exptData_Beam.size(), exptData_Target.size()));
  for(size_t s=0;s<scalingParameters.size();s++){
    scalingParameters.at(s).Fit(exptData_Beam, EffectiveCrossSection_Beam,
                                exptData_Target, EffectiveCrossSection_Target,
                                expt_weights, scaling);
  }
  
  double			beamchisq = 0;
  std::vector<double>	beamexptchisq;
  double			targchisq = 0;
  std::vector<double>	targexptchisq;

  if (verbosity > 1) { PrintYieldHeader(0); }
  chisq += CompareYields(exptData_Beam, EffectiveCrossSection_Beam, scaling, beamexptchisq, beamchisq, NDF_beam, NDF);
  
  if (verbosity > 1) { PrintYieldHeader(1); }
  chisq += CompareYields(exptData_Target, EffectiveCrossSection_Target, scaling, targexptchisq, targchisq, NDF_targ, NDF);

  if(verbosity>1){
    std::cout	<< std::setw(16) << std::left << "Beam expt.:";
    for(size_t i=0;i<beamexptchisq.size();i++)
      std::cout << std::setw(10) << std::left << i;
    std::cout	<< std::endl;
    std::cout	<< std::setw(16) << std::left << "Chi-squared:";
    for(size_t i=0;i<beamexptchisq.size();i++)
      std::cout << std::setw(10) << std::left << std::setprecision(3) << beamexptchisq.at(i);
    std::cout	<< std::endl;
    std::cout	<< std::setw(16) << std::left << "Target expt.:";
    for(size_t i=0;i<targexptchisq.size();i++)
      std::cout << std::setw(10) << std::left << i;
    std::cout	<< std::endl;
    std::cout	<< std::setw(16) << std::left << "Chi-squared:";
    for(size_t i=0;i<targexptchisq.size();i++)
      std::cout << std::setw(10) << std::left << std::setprecision(3) << targexptchisq.at(i);
    std::cout	<< std::endl;
  }

  if(verbosity>1){
    std::cout 	<< std::setw(26) << std::left << "Literature chi-squared: " 
                << std::setw(12) << std::left << std::setprecision(6) << litchisq
                << std::setw(26) << std::left << "Beam chi-squared: " 
                << std::setw(12) << std::left << std::setprecision(6) << beamchisq 
                << std::setw(26) << std::left << "Target chi-squared: " 
                << std::setw(12) << std::left << std::setprecision(6) << targchisq 
                << std::setw(26) << std::left << "Total chi-squared: " 
                << std::setw(12) << std::left << std::setprecision(6) << chisq 
                << std::setw(16) << std::left << "NDF: " 
                << std::setw( 8) << std::left << std::setprecision(6) << NDF
                << std::setw(18) << std::left << "Chisq / NDF: " 
                << std::setw(12) << std::left << std::setprecision(6) << chisq / (float)NDF;
    if(true){
      std::cout	<< std::endl;
      std::cout	<< std::setw(16) << std::left << "NDF lit: " 
                << std::setw( 8) << std::left << std::setprecision(6) << NDF_lit
                << std::setw(16) << std::left << "NDF beam: " 
                << std::setw( 8) << std::left << std::setprecision(6) << NDF_beam
                << std::setw(16) << std::left << "NDF target: " 
                << std::setw( 8) << std::left << std::setprecision(6) << NDF_targ;
    }
    std::cout	<< std::endl;
  }

  iter++;

  Clock::time_point t1 = Clock::now();
  milliseconds ms = std::chrono::duration_cast<milliseconds>(t1-t0);

  if (verbosity == 1) {
    if((iter % 100) == 0 )
      std::cout 	<< std::setw(12) << std::left << iter 
                  << std::setw(13) << std::left << chisq 
                  << std::setw(7)  << std::left << NDF
                  << std::setw(13) << std::left << chisq/(double)NDF 
                  << std::setw(12) << std::left << litchisq 
                  << std::setw(24) << std::left << ms.count() 
                  << "\r" << std::flush;
  }

  return chisq;	

}

void GOSIASimMinFCN::ClearAll() {

  parameters.clear();			
	
  //ME_Beam.clear();
  FE_Beam.clear();
  exptData_Beam.clear();			
  litLifetimes_Beam.clear();			
  litBranchingRatios_Beam.clear();		
  litMixingRatios_Beam.clear();		
  EffectiveCrossSection_Beam.clear();

  //ME_Target.clear();
  FE_Target.clear();
  exptData_Target.clear();			
  litLifetimes_Target.clear();			
  litBranchingRatios_Target.clear();		
  litMixingRatios_Target.clear();		
  EffectiveCrossSection_Target.clear();

}

 void GOSIASimMinFCN::Gosia(int species, int verbosity, int integral) {
  int dirlen = (int)workingDir.size();
  if (species == 0) { //beam
    Nucleus &nucl_b = fNucleus_Beam;
    for (int i=0; i<beamMapping_i.size(); ++i) {
      bst_me[i] = nucl_b.GetMatrixElements().at(beamMapping_l.at(i))[beamMapping_f.at(i)][beamMapping_i.at(i)];
    }
    if (!integral) {
      gosia_(&beam_inputfile,
             workingDir.c_str(), dirlen,
             &all_detectors,
             &bst_me[0],
             &beam_yields,
             verbosity);
    }
    else {
      gosia_(&beam_int_inputfile,
             workingDir.c_str(), dirlen,
             &all_detectors,
             &bst_me[0],
             &beam_yields,
             verbosity);
    }
  }
  else if (species == 1) { //target
    Nucleus &nucl_t = fNucleus_Target;	
    for (int i=0; i<targetMapping_i.size(); ++i) {
      bst_me[i] = nucl_t.GetMatrixElements().at(targetMapping_l.at(i))[targetMapping_f.at(i)][targetMapping_i.at(i)];
    }

    if (!integral) {
      gosia_(&target_inputfile,
             workingDir.c_str(), dirlen,
             &all_detectors,
             &bst_me[0],
             &target_yields,
             verbosity);
    }
    else {
      gosia_(&target_int_inputfile,
             workingDir.c_str(), dirlen,
             &all_detectors,
             &bst_me[0],
             &target_yields,
             verbosity);
    }
  }
}

void GOSIASimMinFCN::AddBeamCorrectionFactor(TMatrixD corrFac){
	correctionFactors_Beam.push_back(corrFac);
}

void GOSIASimMinFCN::CalcBeamCorrectionFactors() {
    out_yields ptyields;
    out_yields intyields;
    int dirlen = workingDir.size();
    for (int i=0; i<beamMapping_i.size(); ++i) {
      bst_me[i] = fNucleus_Beam.GetMatrixElements().at(beamMapping_l.at(i))[beamMapping_f.at(i)][beamMapping_i.at(i)];
    }

    gosia_(&beam_inputfile,
           workingDir.c_str(), dirlen,
           &all_detectors,
           &bst_me[0],
           &ptyields,
           verbosity);

    gosia_(&beam_int_inputfile,
           workingDir.c_str(), dirlen,
           &all_detectors,
           &bst_me[0],
           &intyields,
           verbosity);

    GOSIAReader gosiaReader_point(&fNucleus_Beam, ptyields);
    GOSIAReader gosiaReader_inti(&fNucleus_Beam, intyields);

  correctionFactors_Beam.clear();
  for(size_t e=0; e < gosiaReader_point.GetGOSIAData().size(); e++){	// Loop over experiments
    size_t 		len = fNucleus_Beam.GetNstates();
		TMatrixD	tmpMat;
		tmpMat.ResizeTo(len,len);
    for(size_t ee=0; ee<gosiaReader_point.GetGOSIAData().at(e).GetData().size();ee++){
      int	init = gosiaReader_inti.GetGOSIAData().at(e).GetDataPoint(ee).GetInitialIndex();
      int	fina = gosiaReader_inti.GetGOSIAData().at(e).GetDataPoint(ee).GetFinalIndex();
      double	point = gosiaReader_point.GetGOSIAData().at(e).GetDataPoint(ee).GetCounts();
      double	inti = gosiaReader_inti.GetGOSIAData().at(e).GetDataPoint(ee).GetCounts();
      //std::cout << point << "   " << inti << std::endl;
      if(point > 1e-8){
        tmpMat[init][fina] = inti/point;
        tmpMat[fina][init] = inti/point;
      }
      else{
        tmpMat[init][fina] = 0;
        tmpMat[fina][init] = 0;
      }
    }
    AddBeamCorrectionFactor(tmpMat);
  }
}
    
void GOSIASimMinFCN::AddTargetCorrectionFactor(TMatrixD corrFac){
	correctionFactors_Target.push_back(corrFac);
}

void GOSIASimMinFCN::CalcTargetCorrectionFactors() {
    out_yields ptyields;
    out_yields intyields;
    int dirlen = workingDir.size();
    for (int i=0; i<targetMapping_i.size(); ++i) {
      bst_me[i] = fNucleus_Target.GetMatrixElements().at(targetMapping_l.at(i))[targetMapping_f.at(i)][targetMapping_i.at(i)];
    }

    gosia_(&target_inputfile,
           workingDir.c_str(), dirlen,
           &all_detectors,
           &bst_me[0],
           &ptyields,
           verbosity);

    gosia_(&target_int_inputfile,
           workingDir.c_str(), dirlen,
           &all_detectors,
           &bst_me[0],
           &intyields,
           verbosity);

    GOSIAReader gosiaReader_point(&fNucleus_Target, ptyields);
    GOSIAReader gosiaReader_inti(&fNucleus_Target, intyields);

  correctionFactors_Target.clear();
  for(size_t e=0; e < gosiaReader_point.GetGOSIAData().size(); e++){	// Loop over experiments
    size_t 		len = fNucleus_Target.GetNstates();
		TMatrixD	tmpMat;
		tmpMat.ResizeTo(len,len);
    for(size_t ee=0; ee<gosiaReader_point.GetGOSIAData().at(e).GetData().size();ee++){
      int	init = gosiaReader_inti.GetGOSIAData().at(e).GetDataPoint(ee).GetInitialIndex();
      int	fina = gosiaReader_inti.GetGOSIAData().at(e).GetDataPoint(ee).GetFinalIndex();
      double	point = gosiaReader_point.GetGOSIAData().at(e).GetDataPoint(ee).GetCounts();
      double	inti = gosiaReader_inti.GetGOSIAData().at(e).GetDataPoint(ee).GetCounts();
      //std::cout << point << "   " << inti << std::endl;
      if(point > 1e-8){
        tmpMat[init][fina] = inti/point;
        tmpMat[fina][init] = inti/point;
      }
      else{
        tmpMat[init][fina] = 0;
        tmpMat[fina][init] = 0;
      }
    }
    AddTargetCorrectionFactor(tmpMat);
  }
}
    
