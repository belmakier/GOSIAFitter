#include "ScalingParameter.h"
#include "ScalingFitFCN.h"

#include "Math/Functor.h"
#include "Math/Factory.h"
#include "Math/Minimizer.h"

void ScalingParameter::SetScalingValue(double v, double vll, double vul){

	scaling 	= v;
	scaling_LL 	= vll;
	scaling_UL 	= vul;	

}

void ScalingParameter::CompareData(const ExperimentData &exptData, const TMatrixD &EffectiveCrossSection, const float &weight,
                                   std::vector<double> &sc_expt, std::vector<double> &sc_expt_unc, std::vector<double> &sc_calc) {
  for(size_t t=0;t<exptData.GetData().size();++t){
    int	index_init 	= exptData.GetData().at(t).GetInitialIndex();
    int	index_final 	= exptData.GetData().at(t).GetFinalIndex();
    double 	calcCounts 	= EffectiveCrossSection[index_final][index_init];
    double 	exptCounts 	= exptData.GetData().at(t).GetCounts();
    double	sigma		= (exptData.GetData().at(t).GetUpUnc() + exptData.GetData().at(t).GetDnUnc())/2.;  // Average uncertainty
    sigma 	/= weight;
    if(sigma > 0 && calcCounts > 0 && exptCounts > 0){
      sc_expt.push_back(exptCounts);
      sc_expt_unc.push_back(sigma);
      sc_calc.push_back(calcCounts);
    }				
  }
  for(size_t t=0;t<exptData.GetDoublet().size();++t){
    int	index_init1 	= exptData.GetDoublet().at(t).GetInitialIndex1();
    int	index_final1 	= exptData.GetDoublet().at(t).GetFinalIndex1();
    int	index_init2 	= exptData.GetDoublet().at(t).GetInitialIndex2();
    int	index_final2 	= exptData.GetDoublet().at(t).GetFinalIndex2();
    double 	calcCounts 	= EffectiveCrossSection[index_final1][index_init1] + EffectiveCrossSection[index_final2][index_init2];
    double 	exptCounts 	= exptData.GetDoublet().at(t).GetCounts();
    double	sigma		= (exptData.GetDoublet().at(t).GetUpUnc() + exptData.GetDoublet().at(t).GetDnUnc())/2.;  // Average uncertainty
    sigma 	/= weight;
    if(sigma > 0 && calcCounts > 0 && exptCounts > 0){
      sc_expt.push_back(exptCounts);
      sc_expt_unc.push_back(sigma);
      sc_calc.push_back(calcCounts);
    }				
  }
  return;
}

void ScalingParameter::Fit(const std::vector<ExperimentData> &exptData_Beam, const std::vector<TMatrixD> &EffectiveCrossSection_Beam,
                           const std::vector<ExperimentData> &exptData_Target, const std::vector<TMatrixD> &EffectiveCrossSection_Target,
                           const std::vector<float> &weights, std::vector<double> &scaling)
{
  std::vector<double>	sc_expt;
  std::vector<double>	sc_expt_unc;
  std::vector<double>	sc_calc;

  for(size_t ss=0;ss<GetExperimentNumbers().size();ss++){
    size_t i = GetExperimentNumbers().at(ss);
    if(weights.at(i) == 0) 
      continue;
      
    if(i < exptData_Beam.size()){
      CompareData(exptData_Beam.at(i), EffectiveCrossSection_Beam.at(i), weights.at(i),
                                          sc_expt, sc_expt_unc, sc_calc);
    }

    if(i < exptData_Target.size()){
      CompareData(exptData_Target.at(i), EffectiveCrossSection_Target.at(i), weights.at(i),
                                          sc_expt, sc_expt_unc, sc_calc);
    }
  }

  if(sc_expt.size() > 0){
    ScalingFitFCN theFCN;

    theFCN.SetData(sc_expt,sc_expt_unc,sc_calc);
		
    ROOT::Math::Minimizer *min =
      ROOT::Math::Factory::CreateMinimizer("Minuit2","Migrad");
    ROOT::Math::Functor f_init(theFCN,1);
    min->SetFunction(f_init);
    min->SetVariable(0,"Scaling",1,0.000001);
    min->SetTolerance(0.001);
    min->Minimize();


    //min->PrintResults();
		
    for(size_t ss=0;ss<GetExperimentNumbers().size();ss++){
      size_t i 	= GetExperimentNumbers().at(ss);
      scaling[i]	= min->X()[0];
    }
      
    delete min;
  }
  else {
    for(size_t ss=0;ss<GetExperimentNumbers().size();ss++){
      size_t i 	= GetExperimentNumbers().at(ss);
      scaling[i]	= 0;
    }
  }
}

