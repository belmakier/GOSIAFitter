#include "GOSIASimFitter.h"
#include "Gosia.h"

#include "TString.h"
#include "TRandom3.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "Math/GSLMinimizer.h"
#include "Math/GSLSimAnMinimizer.h"

GOSIASimFitter::GOSIASimFitter()
{

	ClearAll();

	first		= true;

	fLikelihood	= false;

	fDoFullUnc	= false;

	maxIter		= 500;
	maxCalls	= 500;
	fitTolerance	= 0.001;
	nThreads	= 1;

	verbosity		= 1; //default

	chisq		= -1;

  workingDir = "./";

}

GOSIASimFitter::GOSIASimFitter(const GOSIASimFitter& g) {
  workingDir = g.workingDir;

	index				= g.index;

	chisq				= g.chisq;
	parameters			= g.parameters;			
	par_LL				= g.par_LL;				
	par_UL				= g.par_UL;				

  fittingElements_Beam = g.fittingElements_Beam;
  fittingElements_Target = g.fittingElements_Target;
	scalingParameters		= g.scalingParameters;		

	correctionFactors_Beam		= g.correctionFactors_Beam;		
	correctionFactors_Target	= g.correctionFactors_Target;	

	exptData_Beam			= g.exptData_Beam;			
	exptData_Target			= g.exptData_Target;		

	litLifetimes_Beam		= g.litLifetimes_Beam;		
	litLifetimes_Target		= g.litLifetimes_Target;		

	litBranchingRatios_Beam		= g.litBranchingRatios_Beam;	
	litBranchingRatios_Target	= g.litBranchingRatios_Target;	

	litMixingRatios_Beam		= g.litMixingRatios_Beam;		
	litMixingRatios_Target		= g.litMixingRatios_Target;		

	litMatrixElements_Beam		= g.litMatrixElements_Beam;		
	litMatrixElements_Target	= g.litMatrixElements_Target;	

	theErrorDef			= g.theErrorDef;

	EffectiveCrossSection_Beam	= g.EffectiveCrossSection_Beam;
	EffectiveCrossSection_Target	= g.EffectiveCrossSection_Target;

	fNucleus_Beam			= g.fNucleus_Beam;
	fNucleus_Target			= g.fNucleus_Target;
	fNucleus_Beam_Base		= g.fNucleus_Beam_Base;
	fNucleus_Target_Base		= g.fNucleus_Target_Base;

	maxIter				= g.maxIter;
	maxCalls			= g.maxCalls;
	fitTolerance			= g.fitTolerance;

	nThreads			= g.nThreads;

	first				= g.first;
	verbosity				= g.verbosity;

	covMat				= g.covMat;
	corMat				= g.corMat;

	fDoFullUnc			= g.fDoFullUnc;

	fLikelihood			= g.fLikelihood;

	beamGOSIAFile_inp		= g.beamGOSIAFile_inp;
	targetGOSIAFile_inp		= g.targetGOSIAFile_inp;

	beamGOSIAFile_out		= g.beamGOSIAFile_out;
	targetGOSIAFile_out		= g.targetGOSIAFile_out;

	beamBSTFile			= g.beamBSTFile;
	targetBSTFile			= g.targetBSTFile;

	beamMapping_i			= g.beamMapping_i;
	beamMapping_f			= g.beamMapping_f;
	beamMapping_l			= g.beamMapping_l;

	targetMapping_i			= g.targetMapping_i;
	targetMapping_f			= g.targetMapping_f;
	targetMapping_l			= g.targetMapping_l;

	expt_weights			= g.expt_weights;

  all_detectors = g.all_detectors;

}
GOSIASimFitter& GOSIASimFitter::operator = (const GOSIASimFitter& g){

  workingDir = g.workingDir;
	index				= g.index;

	chisq				= g.chisq;
	parameters			= g.parameters;			
	par_LL				= g.par_LL;				
	par_UL				= g.par_UL;				

  fittingElements_Beam		= g.fittingElements_Beam;
  fittingElements_Target		= g.fittingElements_Target;		
	scalingParameters		= g.scalingParameters;		

	correctionFactors_Beam		= g.correctionFactors_Beam;		
	correctionFactors_Target	= g.correctionFactors_Target;	

	exptData_Beam			= g.exptData_Beam;			
	exptData_Target			= g.exptData_Target;		

	litLifetimes_Beam		= g.litLifetimes_Beam;		
	litLifetimes_Target		= g.litLifetimes_Target;		

	litBranchingRatios_Beam		= g.litBranchingRatios_Beam;	
	litBranchingRatios_Target	= g.litBranchingRatios_Target;	

	litMixingRatios_Beam		= g.litMixingRatios_Beam;		
	litMixingRatios_Target		= g.litMixingRatios_Target;		

	litMatrixElements_Beam		= g.litMatrixElements_Beam;		
	litMatrixElements_Target	= g.litMatrixElements_Target;	

	theErrorDef			= g.theErrorDef;

	EffectiveCrossSection_Beam	= g.EffectiveCrossSection_Beam;
	EffectiveCrossSection_Target	= g.EffectiveCrossSection_Target;

	fNucleus_Beam			= g.fNucleus_Beam;
	fNucleus_Target			= g.fNucleus_Target;
	fNucleus_Beam_Base		= g.fNucleus_Beam_Base;
	fNucleus_Target_Base		= g.fNucleus_Target_Base;

	maxIter				= g.maxIter;
	maxCalls			= g.maxCalls;
	fitTolerance			= g.fitTolerance;

	nThreads			= g.nThreads;

	first				= g.first;
	verbosity				= g.verbosity;

	covMat				= g.covMat;
	corMat				= g.corMat;

	fDoFullUnc			= g.fDoFullUnc;

	fLikelihood			= g.fLikelihood;

	beamGOSIAFile_inp		= g.beamGOSIAFile_inp;
	targetGOSIAFile_inp		= g.targetGOSIAFile_inp;

	beamGOSIAFile_out		= g.beamGOSIAFile_out;
	targetGOSIAFile_out		= g.targetGOSIAFile_out;

	beamBSTFile			= g.beamBSTFile;
	targetBSTFile			= g.targetBSTFile;

	beamMapping_i			= g.beamMapping_i;
	beamMapping_f			= g.beamMapping_f;
	beamMapping_l			= g.beamMapping_l;

	targetMapping_i			= g.targetMapping_i;
	targetMapping_f			= g.targetMapping_f;
	targetMapping_l			= g.targetMapping_l;

	expt_weights			= g.expt_weights;

  all_detectors = g.all_detectors;

	return *this;

}

void GOSIASimFitter::UpdateMEs() {
  for(unsigned int i=0;i<fittingElements_Beam.size();i++){
    double *par;
    int parct;      
    fittingElements_Beam[i]->Propagate(fNucleus_Beam,
                                   par, parct,
                                   1);
  }

  for(unsigned int i=0;i<fittingElements_Target.size();i++){
    double *par;
    int parct;      
    fittingElements_Target[i]->Propagate(fNucleus_Target,
                                     par, parct,
                                     1);
  }
}

void GOSIASimFitter::WriteBST() {
  std::ofstream	beam_bst(workingDir+"/"+beamBSTFile);
	for(size_t i=0;i<beamMapping_i.size();i++){
		beam_bst << fNucleus_Beam.GetMatrixElements().at(beamMapping_l.at(i))[beamMapping_f.at(i)][beamMapping_i.at(i)] << "\n";
	}
	beam_bst.close();
  std::ofstream	target_bst(workingDir+"/"+targetBSTFile);
	for(size_t i=0;i<targetMapping_i.size();i++){
		target_bst << fNucleus_Target.GetMatrixElements().at(targetMapping_l.at(i))[targetMapping_f.at(i)][targetMapping_i.at(i)] << "\n";
	}
	target_bst.close();
}

std::vector<double> GOSIASimFitter::GetBeamMEs() {
  std::vector<double> retval;
	for(size_t i=0;i<beamMapping_i.size();i++){
	  retval.push_back(fNucleus_Beam.GetMatrixElements().at(beamMapping_l.at(i))[beamMapping_f.at(i)][beamMapping_i.at(i)]);
	}
  return retval;
}

std::vector<double> GOSIASimFitter::GetTargetMEs() {
  std::vector<double> retval;
	for(size_t i=0;i<targetMapping_i.size();i++){
		retval.push_back(fNucleus_Target.GetMatrixElements().at(targetMapping_l.at(i))[targetMapping_f.at(i)][targetMapping_i.at(i)]);
	}
  return retval;
}

void GOSIASimFitter::DoFit(const char* method, const char *algorithm, ROOT::Math::MinimizerOptions *opt){

	GOSIASimMinFCN theFCN(exptData_Beam,exptData_Target);

  theFCN.SetWorkingDir(workingDir);
	theFCN.SetBeamGOSIAInput(GetBeamGOSIAInput());
	theFCN.SetBeamGOSIAOutput(GetBeamGOSIAOutput());
	theFCN.SetTargetGOSIAInput(GetTargetGOSIAInput());
	theFCN.SetTargetGOSIAOutput(GetTargetGOSIAOutput());
	theFCN.SetBeamBST(GetBeamBST());
	theFCN.SetTargetBST(GetTargetBST());
	theFCN.SetBeamMapping(beamMapping_i,beamMapping_f,beamMapping_l);
	theFCN.SetTargetMapping(targetMapping_i,targetMapping_f,targetMapping_l);

  theFCN.SetBeamFittingElements(fittingElements_Beam);
  theFCN.SetTargetFittingElements(fittingElements_Target);
	theFCN.SetScalingParameters(scalingParameters);

	theFCN.SetBeamLitLifetimes(litLifetimes_Beam);
	theFCN.SetBeamLitBranching(litBranchingRatios_Beam);
	theFCN.SetBeamLitMixing(litMixingRatios_Beam);	
	theFCN.SetBeamLitMatrixElements(litMatrixElements_Beam);	
	theFCN.SetTargetLitLifetimes(litLifetimes_Target);
	theFCN.SetTargetLitBranching(litBranchingRatios_Target);
	theFCN.SetTargetLitMixing(litMixingRatios_Target);	
	theFCN.SetTargetLitMatrixElements(litMatrixElements_Target);	

	theFCN.SetBeamNucleus(&fNucleus_Beam);
	theFCN.SetTargetNucleus(&fNucleus_Target);

	theFCN.SetBeamCorrectionFactors(correctionFactors_Beam);
	theFCN.SetTargetCorrectionFactors(correctionFactors_Target);

  theFCN.SetDetectors(all_detectors);
  theFCN.SetInputFiles(beam_inputfile, target_inputfile);
  theFCN.SetIntInputFiles(beam_int_inputfile, target_int_inputfile);
  
	theFCN.SetIter(maxIter);
	theFCN.SetCalls(maxCalls);

	theFCN.SetNthreads(nThreads);

	theFCN.SetVerbosity(verbosity);

	theFCN.SetupCalculation();

	theFCN.SetLikelihoodFit(fLikelihood);

	size_t	Nexpts = 0;
	if(exptData_Beam.size() > exptData_Target.size())
		Nexpts = exptData_Beam.size();
	else
		Nexpts = exptData_Target.size();

	if(expt_weights.size() != Nexpts){
		expt_weights.resize(Nexpts);
		std::fill(expt_weights.begin(),expt_weights.end(),1);
	}

	theFCN.SetWeights(expt_weights);

  std::vector<std::string> names;
	parameters.clear();
	par_LL.clear();
	par_UL.clear();
  for (unsigned int i=0; i<fittingElements_Beam.size(); ++i) {
    fittingElements_Beam[i]->Populate(parameters, par_LL, par_UL);
    if (!fittingElements_Beam[i]->GetFixed()) {
      if (fittingElements_Beam[i]->GetNPars() == 1) {
        names.push_back((std::string)"Beam-"+fittingElements_Beam[i]->GetType()+"-"+std::to_string(i));
      }
      else {
        std::vector<std::string> alphabet = {"a","b","c","d","e","f","g","h","i","j"};
        for (int j=0; j<fittingElements_Beam[i]->GetNPars(); ++j) {
          names.push_back((std::string)"Beam-"+fittingElements_Beam[i]->GetType()+"-"+std::to_string(i)+alphabet[j]);
        }
      }
    }
  }

  for (unsigned int i=0; i<fittingElements_Target.size(); ++i) {
    fittingElements_Target[i]->Populate(parameters, par_LL, par_UL);
    if (!fittingElements_Target[i]->GetFixed()) {
      if (fittingElements_Target[i]->GetNPars() == 1) {
      names.push_back((std::string)"Targ-"+fittingElements_Target[i]->GetType()+"-"+std::to_string(i));
      }
      else {
        std::vector<std::string> alphabet = {"a","b","c","d","e","f","g","h","i","j"};
        for (int j=0; j<fittingElements_Target[i]->GetNPars(); ++j) {
          names.push_back((std::string)"Targ-"+fittingElements_Target[i]->GetType()+"-"+std::to_string(i)+alphabet[j]);
        }
      }
    }
  }

  if (verbosity > 0) {
	for(unsigned int i=0;i<parameters.size();i++) {
		std::cout 	<< std::setw(13) << std::left << parameters.at(i) 
                << std::setw(4) << std::left << "";
	}
	std::cout << std::endl;
  }
  
	theFCN.SetNpar(parameters.size());
  
  theFCN.SetSimAn(NULL);
  ROOT::Math::Minimizer *min = NULL;
  std::string str("GOSIA_SIMAN");
  std::string meth(method);
  if (!meth.compare(str)) {
    min = new GOSIASimAnMinimizer();
    theFCN.SetSimAn(static_cast<GOSIASimAnMinimizer*>(min)); 
  }
  else {
    min = ROOT::Math::Factory::CreateMinimizer(method, algorithm);
  }
	//ROOT::Math::Minimizer *min = new ROOT::Math::GSLMinimizer("conjugatefr");
 // = new ROOT::Math::GSLSimAnMinimizer();
  //    ROOT::Math::Factory::CreateMinimizer(method, algorithm);
	//		ROOT::Math::Factory::CreateMinimizer("Minuit2","Migrad");
	//		ROOT::Math::Factory::CreateMinimizer("Minuit2","Simplex");
	//		ROOT::Math::Factory::CreateMinimizer("Minuit2","Combined");
	//		ROOT::Math::Factory::CreateMinimizer("Minuit2","Scan");
	//		ROOT::Math::Factory::CreateMinimizer("Minuit2","Fumili");
	//		ROOT::Math::Factory::CreateMinimizer("GSLMultiMin", "ConjugateFR");
	//		ROOT::Math::Factory::CreateMinimizer("GSLMultiMin", "ConjugatePR");
	//		ROOT::Math::Factory::CreateMinimizer("GSLMultiMin", "BFGS");
	//		ROOT::Math::Factory::CreateMinimizer("GSLMultiMin", "BFGS2");
	//		ROOT::Math::Factory::CreateMinimizer("GSLMultiMin", "SteepestDescent");

	ROOT::Math::Functor f_init(theFCN,parameters.size());

  if (!min) { std::cerr << "Minimizer not initalized!" << std::endl; exit(1); }
  
  if (opt) 
    min->SetOptions(*opt);

  if (verbosity > 0) {
    min->Options().Print();
  }
  
	min->SetErrorDef(1.);
	//if(fLikelihood)
	if(false)
		min->SetErrorDef(0.5);

  if (verbosity>0) {
    std::cout 	<< "Iterations: " 
                << maxIter << std::endl;
    std::cout 	<< "Calls: " 
                << maxCalls << std::endl;
  }

	min->SetMaxFunctionCalls(maxCalls);
	min->SetMaxIterations(maxIter);
	min->SetTolerance(fitTolerance);
	min->SetFunction(f_init);

  if (verbosity>0){
    std::cout	<< "Tolerance: " << min->Tolerance()
              << std::endl;
  }
  
	for(unsigned int i=0; i<parameters.size(); i++){
    min->SetLimitedVariable(i,names.at(i),parameters.at(i),0.0001,par_LL.at(i),par_UL.at(i));
    if (verbosity > 0) {
      std::cout	<< names.at(i) << std::setw(4) << " " << std::setw(10) << std::left << parameters.at(i) << std::endl;
    }
	}

  if (verbosity>0)
    std::cout << std::endl;


	min->SetPrecision(1e-8);

  if (verbosity>0) {
    std::cout	<< "************************************ INITIAL MINIMIZATION ************************************"
              << std::endl;


    if(!(verbosity>1) && !fLikelihood){
      std::cout 	<< std::setw(12) << std::left << "Iteration:" 
                  << std::setw(13) << std::left << "Chi2 value:" 
                  << std::setw(7)  << std::left << "NDF:"
                  << std::setw(13) << std::left << "Red. Chi2:"
                  << std::setw(12) << std::left << "Lit. Chi2:" 
                  << std::setw(24) << std::left << "Processing time: (ms)" 
                  << std::endl;
    }
    else{
      std::cout 	<< std::setw(12) << std::left << "Iteration:" 
                  << std::setw(13) << std::left << "-Ln(L) value:" 
                  << std::setw(7)  << std::left << "NDF:"
                  << std::setw(13) << std::left << "Red. Ln(L):"
                  << std::setw(12) << std::left << "Lit. -Ln(L):" 
                  << std::setw(24) << std::left << "Processing time: (ms)" 
                  << std::endl;
    }
  }
	
	typedef std::chrono::high_resolution_clock Clock;
	typedef std::chrono::milliseconds milliseconds;
	Clock::time_point t0 = Clock::now();

	min->Minimize();

  if (verbosity > 0) {
    std::cout	<< std::endl;
  }

	int 	counter = 0;
	int	status = min->Status();

	while(counter < 5 && status == 3){
		min->Minimize();
		status	= min->Status();
		counter++;
	}

	std::cout << std::endl;
	Clock::time_point t1 = Clock::now();
	milliseconds ms = std::chrono::duration_cast<milliseconds>(t1-t0);

  if (verbosity > 0) {
    std::cout	<< "**************************************** FIT COMPLETE ****************************************"
              << std::endl;
    std::cout << "   Status = " << status << std::endl;
  }

	chisq	= min->MinValue();

	const double	*res = min->X();
	const double	*unc = min->Errors();

  //Update res -> fittingElements
  int parct = 0;
  for (int i=0; i<fittingElements_Beam.size(); ++i) {
    fittingElements_Beam[i]->Update(res,parct);
  }
  for (int i=0; i<fittingElements_Target.size(); ++i) {
    fittingElements_Target[i]->Update(res,parct);
  }
  
  if (verbosity > 0) {
    std::cout	<< "Beam fitting elements:"
              << std::endl;
  
    for (int i=0; i<fittingElements_Beam.size(); ++i) {

      fittingElements_Beam[i]->Print(std::cout);
    }

    std::cout	<< "Target fitting elements:"
              << std::endl;
  
    for (int i=0; i<fittingElements_Target.size(); ++i) {
      fittingElements_Target[i]->Print(std::cout);
    }
  }
  //Update fittingElements -> fNucleus
  UpdateMEs();

  if (verbosity > 0) {
    std::cout	<< std::endl;
  }

  for(unsigned int i=0;i<parameters.size();i++)
    parameters[i] = res[i];

  covMat.ResizeTo(parameters.size(),parameters.size());
  corMat.ResizeTo(parameters.size(),parameters.size());
  for(unsigned int i=0;i<parameters.size();i++){
    for(unsigned int j=i;j<parameters.size();j++){
      covMat[i][j] = min->CovMatrix(i,j);
      covMat[j][i] = min->CovMatrix(j,i);
      corMat[i][j] = min->Correlation(i,j);
      corMat[j][i] = min->Correlation(j,i);
    }
  }

  if(DoFullUncertainty()){

    //min->SetTolerance(0.01);

    if (verbosity > 0) {
      std::cout	<< "************************************** UNCERTAINTY EVAL. **************************************"
                << std::endl;
      std::cout	<< "MINOS uncertainties (asymmetric):"
                << std::endl;
    }
    std::vector<double> errLowVec, errUpVec;
    for(unsigned int i=0;i<parameters.size();i++){
      double errLow, errUp;
      min->GetMinosError(i,errLow,errUp);
      errLowVec.push_back(errLow);
      errUpVec.push_back(errUp);
    }
    if (verbosity > 0) {
      std::cout	<< "Correlated uncertainty calculation completed"
                << std::endl;
      std::cout	<< std::setw(14) << std::left << "Parameter" 
                << std::setw(14) << std::left << "Value" 
                << std::setw(14) << "+" 
                << std::setw(3) << "/"
                << std::setw(14) << "-" 
                << std::endl;
    
      for(unsigned int i=0;i<parameters.size();i++){      
        std::cout	<< std::setw(14) << std::left << min->VariableName(i)
                  << std::setw(14) << std::left << parameters[i]
                  << std::setw(14) << errUpVec[i] 
                  << std::setw(3) << ""
                  << std::setw(14) << errLowVec[i]
                  << std::endl;
      }
    }

    for(unsigned int i=0;i<parameters.size();i++){
      for(unsigned int j=i;j<parameters.size();j++){
        covMat[i][j] = min->CovMatrix(i,j);
        covMat[j][i] = min->CovMatrix(j,i);
        corMat[i][j] = min->Correlation(i,j);
        corMat[j][i] = min->Correlation(j,i);
      }
    }

  }

  if (verbosity > 0) {
    std::cout	<< "\n"
              << "******** Covariance matrix **********\n"
              << std::endl;

    covMat.Print();

    std::cout	<< "\n"
              << "******** Correlation matrix **********\n"
              << std::endl;

    corMat.Print();


    std::cout	<< "\n"
              << "******** Final Chi-Squared ***********\n"
              << std::endl;

    std::cout	<< chisq
              << std::endl;
  }
  
  delete min;

  if (verbosity > 0) {
    std::cout	<< "********     Complete      ***********\n"
              << std::endl;
  }

}

void GOSIASimFitter::CreateScalingParameter(std::vector<int> expnum){

	ScalingParameter tmpScaling;
	tmpScaling.SetExperimentVector(expnum);

	scalingParameters.push_back(tmpScaling);

}

void GOSIASimFitter::AddBeamFittingMatrixElement(std::string name, int lambda, int init, int fin, double ME, double LL, double UL, bool fx){
	MatrixElement *tmpME = new MatrixElement(name, fittingElements_Beam.size(),lambda,init,fin,ME,LL,UL,fx);
  fittingElements_Beam.push_back(tmpME);
}


void GOSIASimFitter::SetBeamFittingMatrixElement(std::string name, double ME, double LL, double UL){  
    for (int j=0; j<fittingElements_Beam.size(); ++j) {
      if (fittingElements_Beam[j]->GetName() == name) {
        MatrixElement *me = dynamic_cast<MatrixElement*>(fittingElements_Beam[j]);
        me->SetMatrixElement(ME);
        me->SetMatrixElementUpperLimit(UL);
        me->SetMatrixElementLowerLimit(LL);
      }
    }
}


void GOSIASimFitter::AddBeamRelativeMatrixElement(std::string name, int lambda, int init,int fin,int lambda2,int init2,int fin2 ,double ME, double ME_LL, double ME_UL, bool fixed, bool relative, bool arctan) {
  double rel = 1.0;
  double rel_ll = 1.0; 
  double rel_ul = 1.0;
  if (relative == true ) {
    rel = ME;
    rel_ll = ME_LL;
    rel_ul = ME_UL;
  }
  else {
    if (arctan == false ) {
      rel = ME/fNucleus_Beam.GetMatrixElements().at(lambda2 )[init2][fin2];
      rel_ll = ME_LL/fNucleus_Beam.GetMatrixElements().at(lambda2)[init2][fin2];
      rel_ul = ME_UL/fNucleus_Beam.GetMatrixElements().at(lambda2)[init2][fin2];
    }
    else {
      rel = std::atan(ME/fNucleus_Beam.GetMatrixElements().at(lambda2 )[init2][fin2]);
      rel_ll = std::atan(ME_LL/fNucleus_Beam.GetMatrixElements().at(lambda2)[init2][fin2]);
      rel_ul = std::atan(ME_UL/fNucleus_Beam.GetMatrixElements().at(lambda2)[init2][fin2]);
    }
  }    
  RelativeMatrixElement *tmpME = new RelativeMatrixElement(name, fittingElements_Beam.size(),lambda,init,fin,lambda2,init2,fin2,rel,rel_ll,rel_ul,fixed, arctan);
  fittingElements_Beam.push_back(tmpME);
}

void GOSIASimFitter::AddBeamLifetimeMixingElement(std::string name, int init, int fin, int l1, int l2, 
                                                       double wth, double wth_ll, double wth_ul,
                                                       double mix, double mix_ll, double mix_ul, bool fix){
  LifetimeMixingElement *tmpFE = new LifetimeMixingElement(name,init,fin,l1,l2,
                                                           wth,wth_ll, wth_ul,
                                                           mix,mix_ll,mix_ul,fix);
  fittingElements_Beam.push_back(tmpFE);
}

void GOSIASimFitter::AddBeamRelLtMixElement(std::string name, int init, int fin, int l1, int l2, 
                                            double rel_wth, double rel_wth_ll, double rel_wth_ul,
                                            double mix, double mix_ll, double mix_ul,
                                            int init_ref, int final_ref,
                                            bool fix) {
  RelLtMixElement *tmpFE = new RelLtMixElement(name,init,fin,l1,l2,
                                               rel_wth,rel_wth_ll, rel_wth_ul,
                                               mix,mix_ll,mix_ul,
                                               init_ref, final_ref,
                                               fix);
  fittingElements_Beam.push_back(tmpFE);
}

void GOSIASimFitter::AddBeamRelMatWidthElement(std::string name, int l, int init, int fina,
                                               double relmat, double relmat_ll, double relmat_ul,
                                               int init_ref, int final_ref,
                                               bool fx) {
  RelMatWidthElement *tmpFE = new RelMatWidthElement(name, l, init, fina,
                                                     relmat, relmat_ll, relmat_ul,
                                                     init_ref, final_ref,
                                                     fx);
  fittingElements_Beam.push_back(tmpFE);    
}


void GOSIASimFitter::RemoveBeamRelativeMatrixElement(std::string name) {
  for (int i=0; i<fittingElements_Beam.size(); ++i) {
    if (fittingElements_Beam[i]->GetType() != "RelME") { continue; }
    
    if (fittingElements_Beam[i]->GetName() == name) {
      fittingElements_Beam.erase(fittingElements_Beam.begin() + i);
      break;
    }
  }
}

void GOSIASimFitter::AddTargetFittingMatrixElement(std::string name, int lambda, int init, int fin, double ME, double LL, double UL){
	MatrixElement *tmpME = new MatrixElement(name, fittingElements_Target.size(),lambda,init,fin,ME,LL,UL);
  fittingElements_Target.push_back(tmpME);
}

void GOSIASimFitter::AddTargetRelativeMatrixElement(std::string name, int lambda, int init,int fin,int lambda2,int init2,int fin2 ,double ME, double ME_LL, double ME_UL, bool fixed, bool relative, bool arctan) {
  double rel = 1.0;
  double rel_ll = 1.0; 
  double rel_ul = 1.0;
  if (relative == true ) {
    rel = ME;
    rel_ll = ME_LL;
    rel_ul = ME_UL;
  }
  else {
    if (arctan == false) {
      rel = ME/fNucleus_Target.GetMatrixElements().at(lambda2)[init2][fin2];
      rel_ll = ME_LL/fNucleus_Target.GetMatrixElements().at(lambda2)[init2][fin2];
      rel_ul = ME_UL/fNucleus_Target.GetMatrixElements().at(lambda2)[init2][fin2];
    }
    else {
      rel = std::atan(ME/fNucleus_Target.GetMatrixElements().at(lambda2)[init2][fin2]);
      rel_ll = std::atan(ME_LL/fNucleus_Target.GetMatrixElements().at(lambda2)[init2][fin2]);
      rel_ul = std::atan(ME_UL/fNucleus_Target.GetMatrixElements().at(lambda2)[init2][fin2]);
    }
  }    
  RelativeMatrixElement *tmpME = new RelativeMatrixElement(name, fittingElements_Target.size(),lambda,init,fin,lambda2,init2,fin2,rel,rel_ll,rel_ul,fixed, arctan);
  fittingElements_Target.push_back(tmpME);
}

void GOSIASimFitter::AddWeightingFactor(float f){
	expt_weights.push_back(f);
}
void GOSIASimFitter::SetWeightingFactor(int i, float f){
	if(i < (int)expt_weights.size())
		expt_weights.at(i) = f;
	else
		std::cout << "Outside vector range" << std::endl;
}

void GOSIASimFitter::AddBeamCorrectionFactor(TMatrixD corrFac){
	correctionFactors_Beam.push_back(corrFac);
}
void GOSIASimFitter::SetBeamCorrectionFactor(int i, TMatrixD corrFac){
	if((i < (int)correctionFactors_Beam.size()))
		correctionFactors_Beam.at(i) = corrFac;
	else
		std::cout << "Outside vector range" << std::endl;
}
void GOSIASimFitter::AddTargetCorrectionFactor(TMatrixD corrFac){
	correctionFactors_Target.push_back(corrFac);
}
void GOSIASimFitter::SetTargetCorrectionFactor(int i, TMatrixD corrFac){
	if((i < (int)correctionFactors_Target.size()))
		correctionFactors_Target.at(i) = corrFac;
	else
		std::cout << "Outside vector range" << std::endl;
}

void GOSIASimFitter::DefineExperiment(double thetacm){
	ExperimentData tmpExp;
	tmpExp.SetThetaCM(thetacm);
	exptData_Beam.push_back(tmpExp);		
	exptData_Target.push_back(tmpExp);		
}
void GOSIASimFitter::AddBeamData(int nExpt, int init, int fin, double counts, double unc){
	exptData_Beam.at(nExpt).AddData(init,fin,counts,unc);
} 
void GOSIASimFitter::AddTargetData(int nExpt, int init, int fin, double counts, double unc){
	exptData_Target.at(nExpt).AddData(init,fin,counts,unc);
}   

void GOSIASimFitter::AddBeamLifetime(int index, double lifetime, double unc){
	LitLifetime tmpLifetime(index,lifetime,unc);
	litLifetimes_Beam.push_back(tmpLifetime);  
}

void GOSIASimFitter::AddBeamHalfLife(int index, double halflife, double unc){
	LitLifetime tmpLifetime(index,halflife/TMath::Log(2),unc/TMath::Log(2));
	litLifetimes_Beam.push_back(tmpLifetime);  
}
void GOSIASimFitter::AddBeamBranchingRatio(int index_I1, int index_F1, int index_F2, double br, double unc){
	LitBranchingRatio tmpBR(index_I1,index_F1,index_F2,br,unc);
	litBranchingRatios_Beam.push_back(tmpBR);	
}  
void GOSIASimFitter::AddBeamMixingRatio(int index_I, int index_F, double delta, double unc){
	LitMixingRatio tmpMR(index_I,index_F,delta,unc);
	litMixingRatios_Beam.push_back(tmpMR);	
}
void GOSIASimFitter::AddBeamMatrixElement(int mult, int index_I, int index_F, double me, double unc, int sign){
	LitMatrixElement tmpME(mult,index_I,index_F,me,unc,sign);
	litMatrixElements_Beam.push_back(tmpME);	
}
void GOSIASimFitter::AddTargetLifetime(int index, double lifetime, double unc){
	LitLifetime tmpLifetime(index,lifetime,unc);
	litLifetimes_Target.push_back(tmpLifetime);	
}
void GOSIASimFitter::AddTargetHalfLife(int index, double halflife, double unc){
	LitLifetime tmpLifetime(index,halflife/TMath::Log(2),unc/TMath::Log(2));
	litLifetimes_Target.push_back(tmpLifetime);  
}
void GOSIASimFitter::AddTargetBranchingRatio(int index_I1, int index_F1, int index_F2, double br, double unc){
	LitBranchingRatio tmpBR(index_I1,index_F1,index_F2,br,unc);
	litBranchingRatios_Target.push_back(tmpBR);	
}  
void GOSIASimFitter::AddTargetMixingRatio(int index_I, int index_F, double delta, double unc){
	LitMixingRatio tmpMR(index_I,index_F,delta,unc);
	litMixingRatios_Target.push_back(tmpMR);	
}    
void GOSIASimFitter::AddTargetMatrixElement(int mult, int index_I, int index_F, double me, double unc, int sign){
	LitMatrixElement tmpME(mult, index_I,index_F,me,unc,sign);
	litMatrixElements_Target.push_back(tmpME);	
}    

void GOSIASimFitter::ClearAll(){

	index.clear();
	parameters.clear();
			
	fittingElements_Beam.clear();		
	exptData_Beam.clear();			
	litLifetimes_Beam.clear();			
	litBranchingRatios_Beam.clear();		
	litMixingRatios_Beam.clear();		
	litMatrixElements_Beam.clear();		
	EffectiveCrossSection_Beam.clear();
			
  fittingElements_Target.clear();		
	exptData_Target.clear();			
	litLifetimes_Target.clear();			
	litBranchingRatios_Target.clear();		
	litMixingRatios_Target.clear();		
	litMatrixElements_Target.clear();		
	EffectiveCrossSection_Target.clear();

}

void GOSIASimFitter::Print() const{

	if(exptData_Beam.size()>0){
		std::cout 	<< "\n\n"
                << "Experimental data (beam):" << std::endl;
		
		std::cout 	<< exptData_Beam.size() << " experiments" << std::endl;
		
		for(unsigned int i=0;i<exptData_Beam.size();i++){
			std::cout	<< "Experiment " << i+1 << std::endl;
			std::cout	<< "Theta [CM]: " << exptData_Beam.at(i).GetThetaCM() << std::endl;
			std::cout 	<< std::setw(15) << std::left << "Init. index:" 
                  << std::setw(15) << std::left << "Final index:"
                  << std::setw(15) << std::left << "Init. J:" 
                  << std::setw(15) << std::left << "Final J:"
                  << std::setw(10) << std::left << "Counts:"
                  << std::setw(10) << std::left << "Unc:"
                  << std::endl;
			for(unsigned int t=0;t<exptData_Beam.at(i).GetData().size();t++){
				std::cout 	<< std::setw(15) << std::left << exptData_Beam.at(i).GetDataPoint(t).GetInitialIndex()
                    << std::setw(15) << std::left << exptData_Beam.at(i).GetDataPoint(t).GetFinalIndex()
                    << std::setw(15) << std::left << fNucleus_Beam.GetLevelJ().at(exptData_Beam.at(i).GetDataPoint(t).GetInitialIndex())
                    << std::setw(15) << std::left << fNucleus_Beam.GetLevelJ().at(exptData_Beam.at(i).GetDataPoint(t).GetFinalIndex())
                    << std::setw(10) << std::left << exptData_Beam.at(i).GetDataPoint(t).GetCounts()
                    << std::setw(10) << std::left << exptData_Beam.at(i).GetDataPoint(t).GetUpUnc()
                    << std::endl;
			}			
		}
	}
	else	
		std::cout << "No experimental data declared" << std::endl;

	if(litLifetimes_Beam.size()>0){
		std::cout	<< "\n\n"
              << "Literature lifetimes (beam):"
              << std::endl;
		std::cout	<< std::setw(8)  << std::left << "Index"
              << std::setw(6)  << std::left << "J:"
              << std::setw(15) << std::left << "Lifetime (ps)" 
              << std::setw(15) << std::left << "Uncertainty:"
              << std::endl;
		for(unsigned int i=0;i<litLifetimes_Beam.size();i++){
			std::cout 	<< std::setw(8)  << std::left << litLifetimes_Beam.at(i).GetIndex()
                  << std::setw(6)  << std::left << fNucleus_Beam.GetLevelJ().at(litLifetimes_Beam.at(i).GetIndex())
                  << std::setw(15) << std::left << litLifetimes_Beam.at(i).GetLifetime()
                  << std::setw(15) << std::left << litLifetimes_Beam.at(i).GetUpUnc()
                  << std::endl;
		}
	}

	if(litBranchingRatios_Beam.size()>0){
		std::cout	<< "\n\n"
              << "Literature Branching Ratios (beam):"
              << std::endl;
		std::cout	<< std::setw(15) << std::left << "Init. Index"
              << std::setw(15) << std::left << "Final Index 1"
              << std::setw(15) << std::left << "Final Index 2"
              << std::setw(10) << std::left << "J init:"
              << std::setw(12) << std::left << "J final 1:"
              << std::setw(12) << std::left << "J final 2:"
              << std::setw(17) << std::left << "Branching Ratio" 
              << std::setw(15) << std::left << "Uncertainty:"
              << std::endl;
		for(unsigned int i=0;i<litBranchingRatios_Beam.size();i++){
			std::cout 	<< std::setw(15) << std::left << litBranchingRatios_Beam.at(i).GetInitialIndex()
                  << std::setw(15) << std::left << litBranchingRatios_Beam.at(i).GetFinalIndex_1()
                  << std::setw(15) << std::left << litBranchingRatios_Beam.at(i).GetFinalIndex_2()
                  << std::setw(10) << std::left << fNucleus_Beam.GetLevelJ().at(litBranchingRatios_Beam.at(i).GetInitialIndex())
                  << std::setw(12) << std::left << fNucleus_Beam.GetLevelJ().at(litBranchingRatios_Beam.at(i).GetFinalIndex_1())
                  << std::setw(12) << std::left << fNucleus_Beam.GetLevelJ().at(litBranchingRatios_Beam.at(i).GetFinalIndex_2())
                  << std::setw(17) << std::left << litBranchingRatios_Beam.at(i).GetBranchingRatio()
                  << std::setw(15) << std::left << litBranchingRatios_Beam.at(i).GetUpUnc()
                  << std::endl;
		}
	}
	
	if(litMixingRatios_Beam.size()>0){
		std::cout	<< "\n\n"
              << "Literature Mixing Ratios (beam):"
              << std::endl;
		std::cout	<< std::setw(15) << std::left << "Init. Index"
              << std::setw(15) << std::left << "Final Index"
              << std::setw(10) << std::left << "J init:"
              << std::setw(10) << std::left << "J final:"
              << std::setw(14) << std::left << "Mixing Ratio" 
              << std::setw(15) << std::left << "Uncertainty:"
              << std::endl;
		for(unsigned int i=0;i<litMixingRatios_Beam.size();i++){
			std::cout 	<< std::setw(15) << std::left << litMixingRatios_Beam.at(i).GetInitialIndex()
                  << std::setw(15) << std::left << litMixingRatios_Beam.at(i).GetFinalIndex()
                  << std::setw(10) << std::left << fNucleus_Beam.GetLevelJ().at(litMixingRatios_Beam.at(i).GetInitialIndex())
                  << std::setw(10) << std::left << fNucleus_Beam.GetLevelJ().at(litMixingRatios_Beam.at(i).GetFinalIndex())
                  << std::setw(14) << std::left << litMixingRatios_Beam.at(i).GetMixingRatio()
                  << std::setw(15) << std::left << litMixingRatios_Beam.at(i).GetUpUnc()
                  << std::endl;
		}
	}

	std::string	mult[8] = {"E1","E2","E3","E4","E5","E6","M1","M2"};
  
	if(exptData_Target.size()>0){
		std::cout 	<< "\n\n"
                << "Experimental data (target):" << std::endl;
		
		std::cout 	<< exptData_Target.size() << " experiments" << std::endl;
		
		for(unsigned int i=0;i<exptData_Target.size();i++){
			std::cout	<< "Experiment " << i+1 << std::endl;
			std::cout	<< "Theta [CM]: " << exptData_Target.at(i).GetThetaCM() << std::endl;
			std::cout 	<< std::setw(15) << std::left << "Init. index:" 
                  << std::setw(15) << std::left << "Final index:"
                  << std::setw(15) << std::left << "Init. J:" 
                  << std::setw(15) << std::left << "Final J:"
                  << std::setw(10) << std::left << "Counts:"
                  << std::setw(10) << std::left << "Unc:"
                  << std::endl;
			for(unsigned int t=0;t<exptData_Target.at(i).GetData().size();t++){
				std::cout 	<< std::setw(15) << std::left << exptData_Target.at(i).GetDataPoint(t).GetInitialIndex()
                    << std::setw(15) << std::left << exptData_Target.at(i).GetDataPoint(t).GetFinalIndex()
                    << std::setw(15) << std::left << fNucleus_Target.GetLevelJ().at(exptData_Target.at(i).GetDataPoint(t).GetInitialIndex())
                    << std::setw(15) << std::left << fNucleus_Target.GetLevelJ().at(exptData_Target.at(i).GetDataPoint(t).GetFinalIndex())
                    << std::setw(10) << std::left << exptData_Target.at(i).GetDataPoint(t).GetCounts()
                    << std::setw(10) << std::left << exptData_Target.at(i).GetDataPoint(t).GetUpUnc()
                    << std::endl;
			}			
		}
	}
	else	
		std::cout << "No experimental data declared" << std::endl;

	if(litLifetimes_Target.size()>0){
		std::cout	<< "\n\n"
              << "Literature lifetimes (target):"
              << std::endl;
		std::cout	<< std::setw(8)  << std::left << "Index"
              << std::setw(6)  << std::left << "J:"
              << std::setw(15) << std::left << "Lifetime (ps)" 
              << std::setw(15) << std::left << "Uncertainty:"
              << std::endl;
		for(unsigned int i=0;i<litLifetimes_Target.size();i++){
			std::cout 	<< std::setw(8)  << std::left << litLifetimes_Target.at(i).GetIndex()
                  << std::setw(6)  << std::left << fNucleus_Target.GetLevelJ().at(litLifetimes_Target.at(i).GetIndex())
                  << std::setw(15) << std::left << litLifetimes_Target.at(i).GetLifetime()
                  << std::setw(15) << std::left << litLifetimes_Target.at(i).GetUpUnc()
                  << std::endl;
		}
	}

	if(litBranchingRatios_Target.size()>0){
		std::cout	<< "\n\n"
              << "Literature Branching Ratios (target):"
              << std::endl;
		std::cout	<< std::setw(15) << std::left << "Init. Index"
              << std::setw(15) << std::left << "Final Index 1"
              << std::setw(15) << std::left << "Final Index 2"
              << std::setw(10) << std::left << "J init:"
              << std::setw(12) << std::left << "J final 1:"
              << std::setw(12) << std::left << "J final 2:"
              << std::setw(17) << std::left << "Branching Ratio" 
              << std::setw(15) << std::left << "Uncertainty:"
              << std::endl;
		for(unsigned int i=0;i<litBranchingRatios_Target.size();i++){
			std::cout 	<< std::setw(15) << std::left << litBranchingRatios_Target.at(i).GetInitialIndex()
                  << std::setw(15) << std::left << litBranchingRatios_Target.at(i).GetFinalIndex_1()
                  << std::setw(15) << std::left << litBranchingRatios_Target.at(i).GetFinalIndex_2()
                  << std::setw(10) << std::left << fNucleus_Target.GetLevelJ().at(litBranchingRatios_Target.at(i).GetInitialIndex())
                  << std::setw(12) << std::left << fNucleus_Target.GetLevelJ().at(litBranchingRatios_Target.at(i).GetFinalIndex_1())
                  << std::setw(12) << std::left << fNucleus_Target.GetLevelJ().at(litBranchingRatios_Target.at(i).GetFinalIndex_2())
                  << std::setw(17) << std::left << litBranchingRatios_Target.at(i).GetBranchingRatio()
                  << std::setw(15) << std::left << litBranchingRatios_Target.at(i).GetUpUnc()
                  << std::endl;
		}
	}
	
	if(litMixingRatios_Target.size()>0){
		std::cout	<< "\n\n"
              << "Literature Mixing Ratios (target):"
              << std::endl;
		std::cout	<< std::setw(15) << std::left << "Init. Index"
              << std::setw(15) << std::left << "Final Index"
              << std::setw(10) << std::left << "J init:"
              << std::setw(10) << std::left << "J final:"
              << std::setw(14) << std::left << "Mixing Ratio" 
              << std::setw(15) << std::left << "Uncertainty:"
              << std::endl;
		for(unsigned int i=0;i<litMixingRatios_Target.size();i++){
			std::cout 	<< std::setw(15) << std::left << litMixingRatios_Target.at(i).GetInitialIndex()
                  << std::setw(15) << std::left << litMixingRatios_Target.at(i).GetFinalIndex()
                  << std::setw(10) << std::left << fNucleus_Target.GetLevelJ().at(litMixingRatios_Target.at(i).GetInitialIndex())
                  << std::setw(10) << std::left << fNucleus_Target.GetLevelJ().at(litMixingRatios_Target.at(i).GetFinalIndex())
                  << std::setw(14) << std::left << litMixingRatios_Target.at(i).GetMixingRatio()
                  << std::setw(15) << std::left << litMixingRatios_Target.at(i).GetUpUnc()
                  << std::endl;
		}
	}
  
}

void GOSIASimFitter::WriteBeamFittingParameters(std::ostream &outstream) {
  for (int i=0; i<fittingElements_Beam.size(); ++i) {
    FittingElement *fe = fittingElements_Beam[i];
    for (int j=0; j<fe->GetNPars(); ++j) {
      outstream << std::setw(5) << std::left << i 
                << std::setw(15) << std::left << fe->GetName() 
                << std::setw(15) << std::left << fe->GetType()
                << std::setw(5) << std::left << j
                << std::setw(15) << std::left << fe->GetValue(j) << std::endl;
    }
  }
}

void GOSIASimFitter::WriteBeamFittingParameters(std::string filename) {
  std::ofstream fstream(filename);
  WriteBeamFittingParameters(fstream);
  fstream.close();
}

void GOSIASimFitter::WriteTargetFittingParameters(std::ostream &outstream) {
  for (int i=0; i<fittingElements_Target.size(); ++i) {
    FittingElement *fe = fittingElements_Target[i];
    for (int j=0; j<fe->GetNPars(); ++j) {
      outstream << std::setw(5) << std::left << i 
                << std::setw(15) << std::left << fe->GetName() 
                << std::setw(15) << std::left << fe->GetType()
                << std::setw(5) << std::left << j
                << std::setw(15) << std::left << fe->GetValue(j) << std::endl;
    }
  }
}

void GOSIASimFitter::WriteTargetFittingParameters(std::string filename) {
  std::ofstream fstream(filename);
  WriteTargetFittingParameters(fstream);
  fstream.close();
}

void GOSIASimFitter::ReadBeamFittingParameters(std::string filename) {
  std::ifstream fstream(filename);
  int i, j;
  std::string name, type;
  double value;
  while (fstream >> i >> name >> type >> j >> value) {
    if (i >= fittingElements_Beam.size()) { std::cerr << "Error! Attempt to load parameter value for parameter that does not exist" << std::endl; exit(1); }

    FittingElement *fe = fittingElements_Beam[i];
    if (name.compare(fe->GetName())) { std::cerr << "Error! Parameter " << name << " is not the same as " << fe->GetName() << std::endl; exit(1); }
    if (type.compare(fe->GetType())) { std::cerr << "Error! Parameter type " << type << " is not the same as " << fe->GetType() << std::endl; exit(1); }

    fe->SetValue(j, value);    
  }
  fstream.close();
}

void GOSIASimFitter::ReadTargetFittingParameters(std::string filename) {
  std::ifstream fstream(filename);
  int i, j;
  std::string name, type;
  double value;
  while (fstream >> i >> name >> type >> j >> value) {
    if (i >= fittingElements_Target.size()) { std::cerr << "Error! Attempt to load parameter value for parameter that does not exist" << std::endl; exit(1); }

    FittingElement *fe = fittingElements_Target[i];
    if (name.compare(fe->GetName())) { std::cerr << "Error! Parameter " << name << " is not the same as " << fe->GetName() << std::endl; exit(1); }
    if (type.compare(fe->GetType())) { std::cerr << "Error! Parameter type " << type << " is not the same as " << fe->GetType() << std::endl; exit(1); }

    fe->SetValue(j, value);    
  }
  fstream.close();
}

void GOSIASimFitter::WriteYieldGraphs(TFile *file, std::vector<double> angles, std::vector<double> norms, std::vector<double> &scaling, 
                                      std::vector<TMatrixD> &EffectiveCrossSection, std::vector<TMatrixD> &correctionFactors,
                                      std::vector<ExperimentData> &exptData, std::string species) {
  
  std::map<std::pair<int, int>, TGraph* > calcGraphs;
  std::map<std::pair<int, int>, TGraphErrors* > expGraphs;
  std::map<std::pair<int, int>, TMultiGraph* > Graphs;

  std::map<std::pair<int, int>, TGraph* > normCalcGraphs;
  std::map<std::pair<int, int>, TGraphErrors* > normExpGraphs;
  std::map<std::pair<int, int>, TMultiGraph* > normGraphs;
  
  //generate list of init,final pairs we need graphs for
  //includes individual members of doublets and combined doublets

  //go through list, for each make calcGraph, expGraph
  // -> if doublet then no calcGraph
  // -> if no data then no expGraph

  //go through list again to make multigraphs
  //if both calcGraph + expGraph then add to multiGraph
  //if doublet then fetch component graphs + make summed graph + add to multigraph

	for(unsigned int i=0;i<exptData.size();i++){
    
		double	exptchisq	= 0;
		if(expt_weights.at(i) == 0) 
			continue;
    
		for(unsigned int t=0;t<exptData.at(i).GetData().size();++t){
			double 	tmp 		= 0;
			int	index_init 	= exptData.at(i).GetData().at(t).GetInitialIndex();
			int	index_final 	= exptData.at(i).GetData().at(t).GetFinalIndex();

      TGraph *calcGraph;
      TGraphErrors *expGraph;
      TMultiGraph *graph;
      
      TGraph *normCalcGraph;
      TGraphErrors *normExpGraph;
      TMultiGraph *normGraph;
      
      if (calcGraphs.find({index_init, index_final}) == calcGraphs.end()) {
        calcGraph = new TGraph();
        calcGraphs[{index_init, index_final}] = calcGraph;
        TString calcname;
        calcname.Form("calcYields%s_%i_%i", species.c_str(), index_init, index_final);
        calcGraph->SetName(calcname.Data());
        calcGraph->SetLineWidth(2);
        calcGraph->SetLineColor(kRed);
        file->Append(calcGraph);

        normCalcGraph = new TGraph();
        normCalcGraphs[{index_init, index_final}] = normCalcGraph;
        calcname.Form("calcYieldsNrm%s_%i_%i", species.c_str(), index_init, index_final);
        normCalcGraph->SetName(calcname.Data());
        normCalcGraph->SetLineWidth(2);
        normCalcGraph->SetLineColor(kRed);
        file->Append(normCalcGraph);
      }
      else {
        calcGraph =         calcGraphs[{index_init, index_final}];
        normCalcGraph =         normCalcGraphs[{index_init, index_final}];
      }
      if (expGraphs.find({index_init, index_final}) == expGraphs.end()) {
        expGraph = new TGraphErrors();        
        expGraphs[{index_init, index_final}] = expGraph;
        TString expname;
        expname.Form("expYields%s_%i_%i", species.c_str(), index_init, index_final);
        expGraph->SetMarkerStyle(kFullCircle);
        expGraph->SetName(expname.Data());
        file->Append(expGraph);

        normExpGraph = new TGraphErrors();
        normExpGraphs[{index_init, index_final}] = normExpGraph;
        expname.Form("expYieldsNrm%s_%i_%i", species.c_str(), index_init, index_final);
        normExpGraph->SetMarkerStyle(kFullCircle);
        normExpGraph->SetName(expname.Data());
        file->Append(normExpGraph);
      }
      else {
        expGraph =           expGraphs[{index_init, index_final}];
        normExpGraph =           normExpGraphs[{index_init, index_final}];
      }
      if (Graphs.find({index_init, index_final}) == Graphs.end()) {        
        graph = new TMultiGraph();        
        Graphs[{index_init, index_final}] = graph;
        TString name;
        name.Form("yields%s_%i_%i", species.c_str(), index_init, index_final);
        graph->SetName(name.Data());
        graph->Add(calcGraph, "L");
        graph->Add(expGraph, "P");
        file->Append(graph);


        normGraph = new TMultiGraph();
        normGraphs[{index_init, index_final}] = normGraph;
        name.Form("yieldsNrm%s_%i_%i", species.c_str(), index_init, index_final);
        normGraph->SetName(name.Data());
        normGraph->Add(normCalcGraph, "L");
        normGraph->Add(normExpGraph, "P");
        file->Append(normGraph);
      }
      else {
        graph =           Graphs[{index_init, index_final}];
        normGraph =           normGraphs[{index_init, index_final}];       
      }
      
			double 	calcCounts 	= scaling.at(i) * EffectiveCrossSection.at(i)[index_final][index_init];
			double 	exptCounts 	= exptData.at(i).GetData().at(t).GetCounts();

      calcGraph->AddPoint(i+1, calcCounts);
      expGraph->AddPoint(i+1, exptCounts);
      expGraph->SetPointError(expGraph->GetN()-1, 0, exptData.at(i).GetData().at(t).GetUpUnc());

      normCalcGraph->AddPoint(angles[i], calcCounts/(scaling.at(i)*correctionFactors.at(i)[index_final][index_init]));
      normExpGraph->AddPoint(angles[i], exptCounts/(scaling.at(i)*correctionFactors.at(i)[index_final][index_init]));
      normExpGraph->SetPointError(normExpGraph->GetN()-1, 0, exptData.at(i).GetData().at(t).GetUpUnc()/(scaling.at(i)*correctionFactors.at(i)[index_final][index_init]));
		}
    
		for(unsigned int t=0;t<exptData.at(i).GetDoublet().size();++t){
			double 	tmp 		= 0;
			int	index_init1 	= exptData.at(i).GetDoublet().at(t).GetInitialIndex1();
			int	index_final1 	= exptData.at(i).GetDoublet().at(t).GetFinalIndex1();
			int	index_init2 	= exptData.at(i).GetDoublet().at(t).GetInitialIndex2();
			int	index_final2 	= exptData.at(i).GetDoublet().at(t).GetFinalIndex2();

      TGraph *calcGraph;
      TGraphErrors *expGraph;
      TMultiGraph *graph;
      TGraph *calcGraph1;
      TGraph *calcGraph2;
      
      TGraph *normCalcGraph;
      TGraphErrors *normExpGraph;
      TMultiGraph *normGraph;
      TGraph *normCalcGraph1;
      TGraph *normCalcGraph2;

      if (calcGraphs.find({index_init1, index_final1}) != calcGraphs.end()) {
        calcGraph1 = calcGraphs[{index_init1, index_final1}];
        normCalcGraph1 = normCalcGraphs[{index_init1, index_final1}];
      }
      else {
        calcGraph1 = new TGraph();
        calcGraphs[{index_init1, index_final1}] = calcGraph1;
        TString calcname;
        calcname.Form("calcYields%s_%i_%i", species.c_str(), index_init1, index_final1);
        calcGraph1->SetName(calcname.Data());
        calcGraph1->SetLineWidth(1);
        calcGraph1->SetLineColor(kRed);
        file->Append(calcGraph1);

        normCalcGraph1 = new TGraph();
        normCalcGraphs[{index_init1, index_final1}] = normCalcGraph1;
        calcname.Form("calcYieldsNrm%s_%i_%i", species.c_str(), index_init1, index_final1);
        normCalcGraph1->SetName(calcname.Data());
        normCalcGraph1->SetLineWidth(1);
        normCalcGraph1->SetLineColor(kRed);
        file->Append(normCalcGraph1);
      }
        
      if (calcGraphs.find({index_init2, index_final2}) != calcGraphs.end()) {
        calcGraph2 = calcGraphs[{index_init2, index_final2}];
        normCalcGraph2 = normCalcGraphs[{index_init2, index_final2}];
      }
      else {
        calcGraph2 = new TGraph();
        calcGraphs[{index_init2, index_final2}] = calcGraph2;
        TString calcname;
        calcname.Form("calcYields%s_%i_%i", species.c_str(), index_init2, index_final2);
        calcGraph2->SetName(calcname.Data());
        calcGraph2->SetLineWidth(1);
        calcGraph2->SetLineColor(kRed);
        file->Append(calcGraph2);

        normCalcGraph2 = new TGraph();
        normCalcGraphs[{index_init2, index_final2}] = normCalcGraph2;
        calcname.Form("calcYieldsNrm%s_%i_%i", species.c_str(), index_init2, index_final2);
        normCalcGraph2->SetName(calcname.Data());
        normCalcGraph2->SetLineWidth(1);
        normCalcGraph2->SetLineColor(kRed);
        file->Append(normCalcGraph2);
      }
      
      if (calcGraphs.find({index_init1*100+index_init2, index_final1*100+index_final2}) != calcGraphs.end()) {
        calcGraph = calcGraphs[{index_init1*100+index_init2, index_final1*100+index_final2}];
        normCalcGraph = normCalcGraphs[{index_init1*100+index_init2, index_final1*100+index_final2}];
      }
      else {
        calcGraph = new TGraph();
        calcGraphs[{index_init1*100+index_init2, index_final1*100+index_final2}] = calcGraph;
        TString calcname;
        calcname.Form("calcYields%s_%i_%i", species.c_str(), index_init1*100+index_init2, index_final1*100+index_final2);
        calcGraph->SetName(calcname.Data());
        calcGraph->SetLineWidth(2);
        calcGraph->SetLineColor(kRed);
        file->Append(calcGraph);

        normCalcGraph = new TGraph();
        normCalcGraphs[{index_init1*100+index_init2, index_final1*100+index_final2}] = normCalcGraph;
        calcname.Form("calcYieldsNrm%s_%i_%i", species.c_str(), index_init2, index_final2);
        normCalcGraph->SetName(calcname.Data());
        normCalcGraph->SetLineWidth(2);
        normCalcGraph->SetLineColor(kRed);
        file->Append(normCalcGraph);
      }
      
      if (expGraphs.find({index_init1*100+index_init2, index_final1*100+index_final2}) != expGraphs.end()) {
        expGraph = expGraphs[{index_init1*100+index_init2, index_final1*100+index_final2}];
        normExpGraph = normExpGraphs[{index_init1*100+index_init2, index_final1*100+index_final2}];
      }
      else {
        expGraph = new TGraphErrors();
        expGraphs[{index_init1*100+index_init2, index_final1*100+index_final2}] = expGraph;
        TString expname;
        expname.Form("expYields%s_%i_%i", species.c_str(), index_init1*100+index_init2, index_final1*100+index_final2);
        expGraph->SetName(expname.Data());
        expGraph->SetMarkerStyle(kFullCircle);
        file->Append(expGraph);

        normExpGraph = new TGraphErrors();
        normExpGraphs[{index_init1*100+index_init2, index_final1*100+index_final2}] = normExpGraph;
        expname.Form("expYieldsNrm%s_%i_%i", species.c_str(), index_init1*100+index_init2, index_final1*100+index_final2);
        normExpGraph->SetName(expname.Data());
        normExpGraph->SetMarkerStyle(kFullCircle);
        file->Append(normExpGraph);
      }        
      
      if (Graphs.find({index_init1*100+index_init2, index_final1*100+index_final2}) != Graphs.end()) {
        graph = Graphs[{index_init1*100+index_init2, index_final1*100+index_final2}];
        normGraph = normGraphs[{index_init1*100+index_init2, index_final1*100+index_final2}];
      }
      else {
        graph = new TMultiGraph();
        TString name;
        name.Form("yields%s_%i_%i", species.c_str(), index_init1*100+index_init2, index_final1*100+index_final2);
        graph->SetName(name.Data());
        graph->Add(calcGraph1, "L");
        graph->Add(calcGraph2, "L");
        graph->Add(calcGraph, "L");
        graph->Add(expGraph, "P");
        file->Append(graph);

        normGraph = new TMultiGraph();
        name.Form("yieldsNrm%s_%i_%i", species.c_str(), index_init1*100+index_init2, index_final1*100+index_final2);
        normGraph->SetName(name.Data());
        normGraph->Add(normCalcGraph1, "L");
        normGraph->Add(normCalcGraph2, "L");
        normGraph->Add(normCalcGraph, "L");
        normGraph->Add(normExpGraph, "P");
        file->Append(normGraph);
      }

			double 	calcCounts1 	= scaling.at(i) * EffectiveCrossSection.at(i)[index_final1][index_init1];
      double 	calcCounts2 	= scaling.at(i) * EffectiveCrossSection.at(i)[index_final2][index_init2];
			double 	exptCounts 	= exptData.at(i).GetDoublet().at(t).GetCounts();

      calcGraph1->AddPoint(i+1, calcCounts1);
      calcGraph2->AddPoint(i+1, calcCounts2);
      calcGraph->AddPoint(i+1, calcCounts1+calcCounts2);
      expGraph->AddPoint(i+1, exptCounts);
      expGraph->SetPointError(expGraph->GetN()-1, 0, exptData.at(i).GetData().at(t).GetUpUnc());

      normCalcGraph->AddPoint(angles[i], calcCounts1/(scaling.at(i)*correctionFactors.at(i)[index_final1][index_init1]) + calcCounts2/(scaling.at(i)*correctionFactors.at(i)[index_final2][index_init2]));

      double exptNorm = (scaling.at(i)*correctionFactors.at(i)[index_final1][index_init1])*(scaling.at(i)*correctionFactors.at(i)[index_final2][index_init2]) /
        ((scaling.at(i)*correctionFactors.at(i)[index_final1][index_init1]) * calcCounts2 + (scaling.at(i)*correctionFactors.at(i)[index_final2][index_init2]) * calcCounts1 ) * (calcCounts1 + calcCounts2);
      
      normExpGraph->AddPoint(angles[i], exptCounts/exptNorm);
      normExpGraph->SetPointError(normExpGraph->GetN()-1, 0, exptData.at(i).GetData().at(t).GetUpUnc()/exptNorm);

		}
	}
}

void GOSIASimFitter::WriteYieldGraphs(TFile *file, std::vector<double> angles, std::vector<double> norms) {
  file->cd();
  UpdateMEs();
  std::vector<double> beam_mes = GetBeamMEs();
  std::vector<double> target_mes = GetTargetMEs();
  
  TransitionRates rates_b(&fNucleus_Beam);
  TransitionRates rates_t(&fNucleus_Target);

  RunGosia(beam_inputfile,
           workingDir,
           all_detectors,
           beam_mes,
           beam_yields,
           0);

  RunGosia(target_inputfile,
           workingDir,
           all_detectors,
           target_mes,
           target_yields,
           0);
  
	GOSIAReader	beam_gosiaReader(&fNucleus_Beam, beam_yields);	//	Grab the GOSIA yields
  GOSIAReader target_gosiaReader(&fNucleus_Target, target_yields);	//	Grab the GOSIA yields
  
  std::vector<ExperimentData>	beamCalc	= beam_gosiaReader.GetGOSIAData();
  std::vector<ExperimentData>	targetCalc	= target_gosiaReader.GetGOSIAData();
  EffectiveCrossSection_Beam.clear();	
  EffectiveCrossSection_Target.clear();	

  int dim_b = rates_b.GetBranchingRatios().GetNrows();
  int dim_t = rates_t.GetBranchingRatios().GetNrows();
  for(size_t i=0; i<beamCalc.size(); i++){
    EffectiveCrossSection_Beam.push_back(beamCalc.at(i).GetEffectiveCrossSection(correctionFactors_Beam.at(i), dim_b));
  };
  
  for(size_t i=0; i<targetCalc.size(); i++){
    EffectiveCrossSection_Target.push_back(targetCalc.at(i).GetEffectiveCrossSection(correctionFactors_Target.at(i), dim_t));
  }

  std::vector<double>	scaling;
  scaling.resize(exptData_Beam.size());
  for(size_t s=0;s<scalingParameters.size();s++){
     scalingParameters.at(s).Fit(exptData_Beam, EffectiveCrossSection_Beam,
                                exptData_Target, EffectiveCrossSection_Target,
                                expt_weights, scaling);
  }

  WriteYieldGraphs(file, angles, norms, scaling,
                   EffectiveCrossSection_Beam, correctionFactors_Beam,
                   exptData_Beam, "Beam");

  WriteYieldGraphs(file, angles, norms, scaling,
                   EffectiveCrossSection_Target, correctionFactors_Target,
                   exptData_Target, "Target");

  file->Write();  
  
}

FittingElement* GOSIASimFitter::GetBeamFittingElement(std::string name) {
  for (int j=0; j<fittingElements_Beam.size(); ++j) {
    if (!fittingElements_Beam.at(j)->GetName().compare(name)) {
      return fittingElements_Beam.at(j);
    }
  }
  std::cerr << "Warning! Fitting element " << name << " not found!" << std::endl;
  return NULL;
}

FittingElement* GOSIASimFitter::GetTargetFittingElement(std::string name) {
  for (int j=0; j<fittingElements_Target.size(); ++j) {
    if (!fittingElements_Target.at(j)->GetName().compare(name)) {
      return fittingElements_Target.at(j);
    }
  }
  std::cerr << "Warning! Fitting element " << name << " not found!" << std::endl;
  return NULL;
}

void GOSIASimFitter::FixAllBeamFittingElements() {
  for (int j=0; j<fittingElements_Beam.size(); ++j) {
    fittingElements_Beam[j]->SetFixed(true);
  }
}

void GOSIASimFitter::UnFixAllBeamFittingElements() {
  for (int j=0; j<fittingElements_Beam.size(); ++j) {
    fittingElements_Beam[j]->SetFixed(false);
  }
}

void GOSIASimFitter::FixBeamFittingElement(std::string name) {
  for (int j=0; j<fittingElements_Beam.size(); ++j) {
    if (!fittingElements_Beam[j]->GetName().compare(name)) {
      fittingElements_Beam[j]->SetFixed(true);
    }
  }
}

void GOSIASimFitter::UnFixBeamFittingElement(std::string name) {
  for (int j=0; j<fittingElements_Beam.size(); ++j) {
    if (!fittingElements_Beam[j]->GetName().compare(name)) {
      fittingElements_Beam[j]->SetFixed(false);
    }
  }
}

void GOSIASimFitter::FixBeamFittingElements(std::vector<std::string> names) {
  for (int i = 0; i<names.size(); ++i) {
    for (int j=0; j<fittingElements_Beam.size(); ++j) {
      if (fittingElements_Beam[j]->GetName() == names[i]) {
        fittingElements_Beam[j]->SetFixed(true);
      }
    }
  }
}

void GOSIASimFitter::UnFixBeamFittingElements(std::vector<std::string> names) {
  for (int i = 0; i<names.size(); ++i) {
    for (int j=0; j<fittingElements_Beam.size(); ++j) {
      if (fittingElements_Beam[j]->GetName() == names[i]) {
        fittingElements_Beam[j]->SetFixed(false);
      }
    }
  }
}

void GOSIASimFitter::FixTargetFittingElement(std::string name) {
  for (int j=0; j<fittingElements_Target.size(); ++j) {
    if (!fittingElements_Target[j]->GetName().compare(name)) {
      fittingElements_Target[j]->SetFixed(true);
    }
  }
}

void GOSIASimFitter::UnFixTargetFittingElement(std::string name) {
  for (int j=0; j<fittingElements_Target.size(); ++j) {
    if (!fittingElements_Target[j]->GetName().compare(name)) {
      fittingElements_Target[j]->SetFixed(false);
    }
  }
}

void GOSIASimFitter::FixTargetFittingElements(std::vector<std::string> names) {
  for (int i = 0; i<names.size(); ++i) {
    for (int j=0; j<fittingElements_Target.size(); ++j) {
      if (fittingElements_Target[j]->GetName() == names[i]) {
        fittingElements_Target[j]->SetFixed(true);
      }
    }
  }
}

void GOSIASimFitter::UnFixTargetFittingElements(std::vector<std::string> names) {
  for (int i = 0; i<names.size(); ++i) {
    for (int j=0; j<fittingElements_Target.size(); ++j) {
      if (fittingElements_Target[j]->GetName() == names[i]) {
        fittingElements_Target[j]->SetFixed(false);
      }
    }
  }
}


input_file GOSIASimFitter::ReadInputFile(std::string inpfile) {
  std::ifstream inpfile_ifs(inpfile);
  std::string line;
  input_file inputfile;
  inputfile.nlines = 0;
  while (std::getline(inpfile_ifs, line)) {
    if (line.size() > INPUTFILE_LINELEN) {
      std::cerr << "Error! Line number " << inputfile.nlines+1 << " in input file " << inpfile
                << " is too long: " << line.size() << " > " << INPUTFILE_LINELEN << std::endl;
      exit(1);
    }
    for (int i = 0; i<line.size(); ++i) {
      inputfile.lines[inputfile.nlines][i] = line[i];
    }
    inputfile.lines[inputfile.nlines][line.size()] = '\n';
    inputfile.nlines += 1;
  }
  return inputfile;
}

void GOSIASimFitter::ReadBeamInputFile(std::string inpfile, int type=0) {
  if (type==0) { //point calculation
    beam_inputfile = ReadInputFile(inpfile);
  }
  else if (type==1) { //integral calculation
    beam_int_inputfile = ReadInputFile(inpfile);
  }
  else {
    std::cerr << "Invalid type = " << type << ", 0 for point calculation, 1 for integral" << std::endl;
    exit(1);
  }
}

void GOSIASimFitter::ReadTargetInputFile(std::string inpfile, int type=0) {
  if (type==0) { //point calculation
    target_inputfile = ReadInputFile(inpfile);
  }
  else if (type==1) { //integral calculation
    target_int_inputfile = ReadInputFile(inpfile);
  }
  else {
    std::cerr << "Invalid type = " << type << ", 0 for point calculation, 1 for integral" << std::endl;
    exit(1);
  }
}

void GOSIASimFitter::ReadDetectorFile(std::string dfile) {
  std::ifstream dfile_ifs(dfile);
  dfile_ifs >> all_detectors.ndets;
  for (int i=0; i<all_detectors.ndets; ++i) {
    dfile_ifs >> all_detectors.det[i].radius;
    dfile_ifs >> all_detectors.det[i].energy0;
    for (int j=0; j<8; ++j) {
      dfile_ifs >> all_detectors.det[i].C1[j] >> all_detectors.det[i].C2[j] >> all_detectors.det[i].Qk[j];      
    }
  }
  dfile_ifs.close();
}

void GOSIASimFitter::CalcBeamCorrectionFactors() {
  UpdateMEs();
  out_yields ptyields;
  out_yields intyields;
  std::vector<double> beam_mes = GetBeamMEs();
  RunGosia(beam_inputfile,
           workingDir,
           all_detectors,
           beam_mes,
           ptyields,
           0);

  RunGosia(beam_int_inputfile,
           workingDir,
           all_detectors,
           beam_mes,
           intyields,
           0);

  if (ptyields.nexp != intyields.nexp) {
    std::cerr << "Severe error! GOSIA outputs from integral and point yield have different numbers of experiments. There may be an error in one or both of the GOSIA input files" << std::endl;
    exit(1);
  }
  
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

void GOSIASimFitter::CalcTargetCorrectionFactors() {
  UpdateMEs();
  out_yields ptyields;
  out_yields intyields;
  std::vector<double> target_mes = GetTargetMEs();
  RunGosia(target_inputfile,
           workingDir,
           all_detectors,
           target_mes,
           ptyields,
           0);

  RunGosia(target_int_inputfile,
           workingDir,
           all_detectors,
           target_mes,
           intyields,
           0);

  if (ptyields.nexp != intyields.nexp) {
    std::cerr << "Severe error! GOSIA outputs from integral and point yield have different numbers of experiments. There may be an error in one or both of the GOSIA input files" << std::endl;
    exit(1);
  }
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


void GOSIASimFitter::Kick(int seed) {
  TRandom3 rand;
  rand.SetSeed(seed);
  for (int i=0; i<fittingElements_Beam.size(); ++i) {
    fittingElements_Beam[i]->Kick(rand);
  }

  for (int i=0; i<fittingElements_Target.size(); ++i) {
    fittingElements_Target[i]->Kick(rand);
  }

}
