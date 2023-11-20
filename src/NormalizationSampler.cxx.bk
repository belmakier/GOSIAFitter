#include "NormalizationSampler.h"

void NormalizationSampler::SetupCalculation(){

	size_t	size = 0;
	if(GetBeamData().size() > GetTargetData().size())
		size = GetBeamData().size();
	else		
		size = GetTargetData().size();
	exptIndex.resize(size);
	for(unsigned int i=0;i<GetScalingParameters().size();i++)
		for(unsigned int s=0;s<GetScalingParameters().at(i).GetExperimentNumbers().size();s++)
			exptIndex[GetScalingParameters().at(i).GetExperimentNumbers().at(s)] = (int)i;

}

void NormalizationSampler::Run(double *par){

	Nucleus nucl_b = GetBeamNucleus();
	Nucleus nucl_t = GetTargetNucleus();

	for(unsigned int i=0;i<GetBeamMatrixElements().size();i++){
		nucl_b.SetMatrixElement(GetBeamMatrixElements().at(i).GetLambda(),GetBeamMatrixElements().at(i).GetInitialState(),GetBeamMatrixElements().at(i).GetFinalState(),par[i]);
	}
	for(unsigned int i=0;i<GetTargetMatrixElements().size();i++){
		nucl_t.SetMatrixElement(GetTargetMatrixElements().at(i).GetLambda(),GetTargetMatrixElements().at(i).GetInitialState(),GetTargetMatrixElements().at(i).GetFinalState(),par[i + GetBeamMatrixElements().size()]);
	}

	//	COULEX AND STUFF:
	
	std::ofstream		beam_bst(GetBeamBST());
	for(size_t i=0;i<GetBeamMappingInit().size();i++){
		beam_bst << nucl_b.GetMatrixElements().at(GetBeamMappingLambda().at(i))[GetBeamMappingFinal().at(i)][GetBeamMappingInit().at(i)] << "\n";
	}
	beam_bst.close();
	
	std::ofstream		target_bst(GetTargetBST());
	for(size_t i=0;i<GetTargetMappingInit().size();i++){
		target_bst 	<< nucl_t.GetMatrixElements().at(GetTargetMappingLambda().at(i))[GetTargetMappingFinal().at(i)][GetTargetMappingInit().at(i)] << "\n";
	}
	target_bst.close();

	TransitionRates rates_b(&nucl_b);
	TransitionRates rates_t(&nucl_t);

	std::string	str;
	str = "./gosia < "+ GetBeamGOSIAInput()+"> /dev/null";
	const char* 	c_b = str.c_str();
	system(c_b);		// 	Run the beam GOSIA file
	str = "./gosia < "+ GetTargetGOSIAInput()+"> /dev/null";
	const char* 	c_t = str.c_str();
	system(c_t);	// 	Run the target GOSIA file

	const char	*b_out = GetBeamGOSIAOutput().c_str();
	const char	*t_out = GetTargetGOSIAOutput().c_str();
	GOSIAReader	beam_gosiaReader(&nucl_b,b_out);	//	Grab the GOSIA yields
	GOSIAReader	target_gosiaReader(&nucl_t,t_out);	//	Grab the GOSIA yields

	std::vector<ExperimentData>	beamCalc	= beam_gosiaReader.GetGOSIAData();
	std::vector<ExperimentData>	targetCalc	= target_gosiaReader.GetGOSIAData();
	ClearEffectiveCrossSection_Beam();	
	ClearEffectiveCrossSection_Target();	

	for(size_t i=0; i<beamCalc.size(); i++){
		TMatrixD	tmpMat;
		tmpMat.ResizeTo(rates_b.GetBranchingRatios().GetNrows(),rates_b.GetBranchingRatios().GetNcols());
		size_t	nRows = beamCalc.at(i).GetData().size();
		for(size_t j=0; j<nRows; j++){
			int	init		= beamCalc.at(i).GetData().at(j).GetInitialIndex();
			int	fina		= beamCalc.at(i).GetData().at(j).GetFinalIndex();
			double 	counts 		= beamCalc.at(i).GetData().at(j).GetCounts();
			if(false && i==0){
				std::cout	<< std::setw( 6) << std::left << init
						<< std::setw( 6) << std::left << fina
						<< std::setw(12) << std::left << counts	
						<< std::setw(12) << std::left << GetBeamCorrectionFactors().at(i)[init][fina]
						<< std::setw(12) << std::left << GetBeamCorrectionFactors().at(i)[init][fina] * counts
						<< std::endl;
			}
			tmpMat[fina][init]	= counts * GetBeamCorrectionFactors().at(i)[init][fina];
			tmpMat[init][fina]	= counts * GetBeamCorrectionFactors().at(i)[init][fina];
		}
		AddEffectiveCrossSection_Beam(tmpMat);
	}
	for(size_t i=0; i<targetCalc.size(); i++){
		TMatrixD	tmpMat;
		tmpMat.ResizeTo(rates_t.GetBranchingRatios().GetNrows(),rates_t.GetBranchingRatios().GetNcols());
		size_t	nRows = targetCalc.at(i).GetData().size();
		for(size_t j=0; j<nRows; j++){
			int	init		= targetCalc.at(i).GetData().at(j).GetInitialIndex();
			int	fina		= targetCalc.at(i).GetData().at(j).GetFinalIndex();
			double 	counts 		= targetCalc.at(i).GetData().at(j).GetCounts();
			if(false && i==0){
				std::cout	<< std::setw( 6) << std::left << init
						<< std::setw( 6) << std::left << fina
						<< std::setw(12) << std::left << counts	
						<< std::setw(12) << std::left << GetTargetCorrectionFactors().at(i)[init][fina]
						<< std::setw(12) << std::left << GetTargetCorrectionFactors().at(i)[init][fina] * counts
						<< std::endl;
			}
			tmpMat[fina][init]	= counts * GetTargetCorrectionFactors().at(i)[init][fina];
			tmpMat[init][fina]	= counts * GetTargetCorrectionFactors().at(i)[init][fina];
		}
		AddEffectiveCrossSection_Target(tmpMat);
	}


	std::vector<double>	scaling;
	size_t	size = 0;
	if(GetBeamData().size() > GetTargetData().size())
		size = GetBeamData().size();
	else		
		size = GetTargetData().size();
	scaling.resize(size);
	for(size_t s=0;s<GetScalingParameters().size();s++){
		std::vector<double>	sc_expt;
		std::vector<double>	sc_expt_unc;
		std::vector<double>	sc_calc;
		for(size_t ss=0;ss<GetScalingParameters().at(s).GetExperimentNumbers().size();ss++){
			size_t i = GetScalingParameters().at(s).GetExperimentNumbers().at(ss);
			if(i < GetBeamData().size()){

				for(size_t t=0;t<GetBeamData().at(i).GetData().size();++t){
					int	index_init 	= GetBeamData().at(i).GetData().at(t).GetInitialIndex();
					int	index_final 	= GetBeamData().at(i).GetData().at(t).GetFinalIndex();
					double 	calcCounts 	= GetEffectiveCrossSection_Beam().at(i)[index_final][index_init];
					double 	exptCounts 	= GetBeamData().at(i).GetData().at(t).GetCounts();
					double	sigma		= (GetBeamData().at(i).GetData().at(t).GetUpUnc() + GetBeamData().at(i).GetData().at(t).GetDnUnc())/2.;  // Average uncertainty
					if(sigma > 0 && calcCounts > 0 && exptCounts > 0){
						sc_expt.push_back(exptCounts);
						sc_expt_unc.push_back(sigma);
						sc_calc.push_back(calcCounts);
					}				
				}
			}

			if(i < GetTargetData().size()){
				for(size_t t=0;t<GetTargetData().at(i).GetData().size();++t){
					int	index_init 	= GetTargetData().at(i).GetData().at(t).GetInitialIndex();
					int	index_final 	= GetTargetData().at(i).GetData().at(t).GetFinalIndex();
					double 	calcCounts 	= GetEffectiveCrossSection_Target().at(i)[index_final][index_init];
					double 	exptCounts 	= GetTargetData().at(i).GetData().at(t).GetCounts();
					double	sigma		= (GetTargetData().at(i).GetData().at(t).GetUpUnc() + GetTargetData().at(i).GetData().at(t).GetDnUnc())/2.;  // Average uncertainty
					if(sigma > 0 && calcCounts > 0 && exptCounts > 0){
						sc_expt.push_back(exptCounts);
						sc_expt_unc.push_back(sigma);
						sc_calc.push_back(calcCounts);
					}
				}
			}
		}

		ScalingFitFCN theFCN;

		theFCN.SetData(sc_expt,sc_expt_unc,sc_calc);
	
		ROOT::Math::Minimizer *min =
			ROOT::Math::Factory::CreateMinimizer("Minuit2","Combined");
		ROOT::Math::Functor f_init(theFCN,1);
		min->SetFunction(f_init);
		min->SetVariable(0,"Scaling",1,0.0001);
		min->Minimize();

		
		for(size_t ss=0;ss<GetScalingParameters().at(s).GetExperimentNumbers().size();ss++){
			size_t i 	= GetScalingParameters().at(s).GetExperimentNumbers().at(ss);
			scaling[i]	= min->X()[0];
		}

	}

	Normalizations.push_back(scaling);

}
