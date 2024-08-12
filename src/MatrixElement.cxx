#include "MatrixElement.h"


MatrixElement& MatrixElement::operator = (const MatrixElement& m){

	index			= m.index;
	lambda			= m.lambda;
	initialstate		= m.initialstate;
	finalstate		= m.finalstate;
	matrixElement		= m.matrixElement;
	matrixElement_ll	= m.matrixElement_ll;
	matrixElement_ul 	= m.matrixElement_ul;
  fixed = m.fixed;

	return *this;

}
MatrixElement::MatrixElement(const MatrixElement& m){

	index			= m.index;
	lambda			= m.lambda;
	initialstate		= m.initialstate;
	finalstate		= m.finalstate;
	matrixElement		= m.matrixElement;
	matrixElement_ll	= m.matrixElement_ll;
	matrixElement_ul 	= m.matrixElement_ul;
  fixed = m.fixed;

}
void MatrixElement::Print(std::ostream &out) const {
  out << type << " " << name << "  :  " << matrixElement << std::endl;
}

void MatrixElement::Kick(TRandom3 &rand) {
  if (GetFixed()) { return; }
  double newme;
  if (GetMatrixElementUpperLimit() * GetMatrixElementLowerLimit() < 0) {
    newme = rand.Gaus(GetMatrixElement(), std::max(std::abs(GetMatrixElementLowerLimit()),
                                                      std::abs(GetMatrixElementUpperLimit()))/4.);
    double newsign = rand.Uniform() - 0.5;
    if (newsign != 0) {
      newme = newme * std::abs(newsign)/newsign;
    }
  }
  else {
    newme = rand.Gaus(GetMatrixElement(), (GetMatrixElementUpperLimit()-GetMatrixElementLowerLimit())/4.);
  }

  if (newme > GetMatrixElementUpperLimit()) { newme = GetMatrixElementUpperLimit(); }
  if (newme < GetMatrixElementLowerLimit()) { newme = GetMatrixElementLowerLimit(); }
  SetMatrixElement(newme);
}

void MatrixElement::Propagate(Nucleus &nucl,
                              const double *par, int &parct,
                              int opt) {  
  if (GetFixed() || opt) {
    nucl.SetMatrixElement(GetLambda(),
                          GetInitialState(),
                          GetFinalState(),
                          GetMatrixElement());
  } 
  else {
    nucl.SetMatrixElement(GetLambda(),
                          GetInitialState(),
                          GetFinalState(),
                          par[parct]);
    indx = parct;
    ++parct;
  }
}

void MatrixElement::Update(const double *res, int &parct) {
  if (GetFixed() ) { return; } 
  SetMatrixElement(res[parct]);
  ++parct;
}

void MatrixElement::Populate(std::vector<double> &parameters,
                             std::vector<double> &par_LL,
                             std::vector<double> &par_UL) {
  if (!GetFixed()) {
    parameters.push_back(GetMatrixElement());
    par_LL.push_back(GetMatrixElementLowerLimit());
    par_UL.push_back(GetMatrixElementUpperLimit());
  }    
}
                             

RelativeMatrixElement& RelativeMatrixElement::operator = (const RelativeMatrixElement& m){

	index			= m.index;
  lambda			= m.lambda;
	initialstate		= m.initialstate;
	finalstate		= m.finalstate;
  lambdaRel			= m.lambdaRel;
  initialstateRel		= m.initialstateRel;
	finalstateRel		= m.finalstateRel;
  rel = m.rel;
  relll = m.relll;
  relul = m.relul;
  fixed = m.fixed;
  arctan = m.arctan;

	return *this;

}
RelativeMatrixElement::RelativeMatrixElement(const RelativeMatrixElement& m){

  index			= m.index;
  lambda			= m.lambda;
	initialstate		= m.initialstate;
	finalstate		= m.finalstate;
  lambdaRel			= m.lambdaRel;
  initialstateRel		= m.initialstateRel;
	finalstateRel		= m.finalstateRel;
  rel = m.rel;
  relll = m.relll;
  relul = m.relul;
  fixed = m.fixed;
  arctan = m.arctan;


}
void RelativeMatrixElement::Print(std::ostream &out) const {
  out << type << " " << name << "  :  " << rel << std::endl;
}

void RelativeMatrixElement::Kick(TRandom3 &rand) {
  if (GetFixed()) { return; }
  double newme;
  if (GetRelativeElementUpperLimit() * GetRelativeElementLowerLimit() < 0) {
    newme = rand.Gaus(GetRelativeElement(), std::max(std::abs(GetRelativeElementLowerLimit()),
                                                      std::abs(GetRelativeElementUpperLimit()))/4.);
    double newsign = rand.Uniform() - 0.5;
    if (newsign != 0) {
    newme = newme * std::abs(newsign)/newsign;
    }
  }
  else {
    newme = rand.Gaus(GetRelativeElement(), (GetRelativeElementUpperLimit()-GetRelativeElementLowerLimit())/4.);
  }

  if (newme > GetRelativeElementUpperLimit()) { newme = GetRelativeElementUpperLimit(); }
  if (newme < GetRelativeElementLowerLimit()) { newme = GetRelativeElementLowerLimit(); }
  SetRelativeElement(newme);
}

void RelativeMatrixElement::Propagate(Nucleus &nucl,
                                      const double *par, int &parct,
                                      int opt) {
  double me = nucl.GetMatrixElements().at(GetLambdaRel())[GetInitialStateRel()][GetFinalStateRel()];
  
  if (GetFixed() || opt) {
    if (GetArcTan()) {
      me = me * std::tan(GetRelativeElement());
    }
    else {
      me = me * GetRelativeElement();
    }
  }
  else {
    if (GetArcTan()) {
      me = me * std::tan(par[parct]);
    }
    else {
      me = me * par[parct];
    }
    indx = parct;
    ++parct;
  }
  nucl.SetMatrixElement(GetLambda(), GetInitialState(), GetFinalState(), me);
}

void RelativeMatrixElement::Update(const double *res, int &parct) {
  if (GetFixed() ) { return; } 
  SetRelativeElement(res[parct]);
  ++parct;
}

void RelativeMatrixElement::Populate(std::vector<double> &parameters,
                               std::vector<double> &par_LL,
                               std::vector<double> &par_UL) {
    if (!GetFixed()) {
      parameters.push_back(GetRelativeElement());
      par_LL.push_back(GetRelativeElementLowerLimit());
      par_UL.push_back(GetRelativeElementUpperLimit());
    }    
  }

void LifetimeMixingElement::Propagate(Nucleus &nucl,
                                      const double *par, int &parct,
                                      int opt) {
  //par[parct] = partial decay width in 1./ps
  //par[parct+1] = mixing ratio of lambda1/lambda2
  // so typically lambda1 = E2, lambda2 = M1

  indx = parct;
  double wth;
  double delt;  
  if (GetFixed() || opt) {
    wth = width;
    delt = mixingRatio;
  }
  else {
    wth = par[parct];
    delt = par[parct+1];
    parct += 2;
  }

  double Egamma = nucl.GetLevelEnergies()[initialState] - nucl.GetLevelEnergies()[finalState];
  double decay = wth*std::pow(10,12);
  double factor = 1./std::sqrt((multfactor[lambda1]/std::pow(100,nbarns[lambda1]))/
                            (multfactor[lambda2]/std::pow(100,nbarns[lambda2])));
  double ratio = delt/(factor*std::pow(Egamma, (int)((power[lambda1] - power[lambda2])/2)));

  int Ji = nucl.GetLevelJ()[initialState];
  int sign = 1;
  if (decay != 0) {
    sign = std::abs(decay)/decay;
  }
  double me_2 = std::sqrt(std::abs(decay)/(std::pow(100,nbarns[lambda1]) * (ratio*ratio/(2*Ji+1)) / (multfactor[lambda1] * std::pow(Egamma, -power[lambda1])) + std::pow(100,nbarns[lambda2]) * (1./(2*Ji+1)) / (multfactor[lambda2] * std::pow(Egamma, -power[lambda2]))));
  double me_1 = me_2*ratio;

  if (((int)(std::abs(me_1)/me_1)) != (int)(std::abs(wth)/wth) ) {
    me_1 = -me_1;
    me_2 = -me_2;
  }

  nucl.SetMatrixElement(lambda1, initialState, finalState, me_1);
  nucl.SetMatrixElement(lambda2, initialState, finalState, me_2);
}

void LifetimeMixingElement::Update(const double *res, int &parct) {
  if (GetFixed()) { return; } 
  SetWidth(res[parct]);
  SetMixingRatio(res[parct+1]);
  parct += 2;
}

void LifetimeMixingElement::Populate(std::vector<double> &parameters,
                               std::vector<double> &par_LL,
                               std::vector<double> &par_UL) {
    if (!GetFixed()) {
      parameters.push_back(GetWidth());
      par_LL.push_back(GetWidthLowerLimit());
      par_UL.push_back(GetWidthUpperLimit());

      parameters.push_back(GetMixingRatio());
      par_LL.push_back(GetMixingRatioLowerLimit());
      par_UL.push_back(GetMixingRatioUpperLimit());
    }
}

void LifetimeMixingElement::Print(std::ostream &out) const {
  out << type << " " << name << "  :  " << width << " / " << mixingRatio << std::endl;
}

void LifetimeMixingElement::Kick(TRandom3 &rand) {
  if (GetFixed()) { return; }
  double newwidth;
  if (GetWidthUpperLimit() * GetWidthLowerLimit() < 0) {    
    newwidth = rand.Gaus(GetWidth(), std::max(std::abs(GetWidthLowerLimit()), std::abs(GetWidthUpperLimit()))/4.);
    double newsign = rand.Uniform() - 0.5;
    if (newsign != 0) {
    newwidth = newwidth * std::abs(newsign)/newsign;
    }
  }
  else {
    newwidth = rand.Gaus(GetWidth(), (GetWidthLowerLimit() - GetWidthUpperLimit())/4.);
  }
  double newmix = rand.Gaus(GetMixingRatio(), (GetMixingRatioLowerLimit() - GetMixingRatioUpperLimit())/4.);
  
  if (newwidth > GetWidthUpperLimit()) { newwidth = GetWidthUpperLimit(); }
  if (newwidth < GetWidthLowerLimit()) { newwidth = GetWidthLowerLimit(); }

  if (newmix > GetMixingRatioUpperLimit()) { newmix = GetMixingRatioUpperLimit(); }
  if (newmix < GetMixingRatioLowerLimit()) { newmix = GetMixingRatioLowerLimit(); }
  
  SetWidth(newwidth);
  SetMixingRatio(newmix);
}

void RelLtMixElement::Propagate(Nucleus &nucl,
                                const double *par, int &parct,
                                int opt) {
  //par[parct] = relative partial decay width
  //par[parct+1] = mixing ratio of lambda1/lambda2
  // so typically lambda1 = E2, lambda2 = M1

  indx = parct;
  double wth_rel;
  double delt;  
  if (GetFixed() || opt) {
    wth_rel = widthRel;
    delt = mixingRatio;
  }
  else {
    wth_rel = par[parct];
    delt = par[parct+1];
    parct += 2;
  }

  double wth_ref = 0;
  double Eg_ref = nucl.GetLevelEnergies()[initRef] - nucl.GetLevelEnergies()[finalRef];
  int J_ref = nucl.GetLevelJ()[initRef];
  int sign = 0;
  for (int L=0; L<8; ++L) {
    double me = nucl.GetMatrixElements().at(L)[initRef][finalRef];
    if (me == 0) { continue; }
    if (me != 0 && sign == 0) {
      sign = std::abs(me)/me;
    }
    wth_ref += 1./(multfactor[L] * std::pow(Eg_ref,-power[L]) / (me*me/(2*J_ref+1.)) / std::pow(100,nbarns[L]));
  }
  double decay = wth_rel * wth_ref * sign;

  double Egamma = nucl.GetLevelEnergies()[initialState] - nucl.GetLevelEnergies()[finalState];
  double factor = 1./std::sqrt((multfactor[lambda1]/std::pow(100,nbarns[lambda1]))/
                            (multfactor[lambda2]/std::pow(100,nbarns[lambda2])));
  double ratio = delt/(factor*std::pow(Egamma, (int)((power[lambda1] - power[lambda2])/2)));

  int Ji = nucl.GetLevelJ()[initialState];
  double me_2 = std::sqrt(std::abs(decay)/(std::pow(100,nbarns[lambda1]) * (ratio*ratio/(2*Ji+1)) / (multfactor[lambda1] * std::pow(Egamma, -power[lambda1])) + std::pow(100,nbarns[lambda2]) * (1./(2*Ji+1)) / (multfactor[lambda2] * std::pow(Egamma, -power[lambda2]))));
  double me_1 = me_2*ratio;

  if (((int)(std::abs(me_1)/me_1)) != (int)(std::abs(decay)/decay) ) {
    me_1 = -me_1;
    me_2 = -me_2;
  }

  nucl.SetMatrixElement(lambda1, initialState, finalState, me_1);
  nucl.SetMatrixElement(lambda2, initialState, finalState, me_2);
}

void RelLtMixElement::Update(const double *res, int &parct) {
  if (GetFixed()) { return; } 
  SetWidthRel(res[parct]);
  SetMixingRatio(res[parct+1]);
  parct += 2;
}

void RelLtMixElement::Populate(std::vector<double> &parameters,
                               std::vector<double> &par_LL,
                               std::vector<double> &par_UL) {
    if (!GetFixed()) {
      parameters.push_back(GetWidthRel());
      par_LL.push_back(GetWidthRelLowerLimit());
      par_UL.push_back(GetWidthRelUpperLimit());

      parameters.push_back(GetMixingRatio());
      par_LL.push_back(GetMixingRatioLowerLimit());
      par_UL.push_back(GetMixingRatioUpperLimit());
    }
}

void RelLtMixElement::Print(std::ostream &out) const {
  out << type << " " << name << "  :  " << widthRel << " / " << mixingRatio << std::endl;
}

void RelLtMixElement::Kick(TRandom3 &rand) {
  if (GetFixed()) { return; }
  double newwidth;
  if (GetWidthRelUpperLimit() * GetWidthRelLowerLimit() < 0) {    
    newwidth = rand.Gaus(GetWidthRel(), std::max(std::abs(GetWidthRelLowerLimit()), std::abs(GetWidthRelUpperLimit()))/4.);
    double newsign = rand.Uniform() - 0.5;
    if (newsign != 0) {
    newwidth = newwidth * std::abs(newsign)/newsign;
    }
  }
  else {
    newwidth = rand.Gaus(GetWidthRel(), (GetWidthRelLowerLimit() - GetWidthRelUpperLimit())/4.);
  }
  double newmix = rand.Gaus(GetMixingRatio(), (GetMixingRatioLowerLimit() - GetMixingRatioUpperLimit())/4.);
  
  if (newwidth > GetWidthRelUpperLimit()) { newwidth = GetWidthRelUpperLimit(); }
  if (newwidth < GetWidthRelLowerLimit()) { newwidth = GetWidthRelLowerLimit(); }

  if (newmix > GetMixingRatioUpperLimit()) { newmix = GetMixingRatioUpperLimit(); }
  if (newmix < GetMixingRatioLowerLimit()) { newmix = GetMixingRatioLowerLimit(); }
  
  SetWidthRel(newwidth);
  SetMixingRatio(newmix);
}

void RelMatWidthElement::Propagate(Nucleus &nucl,
                                   const double *par, int &parct, int opt) {
  indx = parct;
  double mat_rel;
  if (GetFixed() || opt) {
    mat_rel = matRel;
  }
  else {
    mat_rel = par[parct];
    parct += 1;
  }

  double wth_ref = 0;
  double Eg_ref = nucl.GetLevelEnergies()[initRef] - nucl.GetLevelEnergies()[finalRef];
  double J_ref = nucl.GetLevelJ()[initRef];
  int sign = 0;
  for (int L=0; L<8; ++L) {
    double me = nucl.GetMatrixElements().at(L)[initRef][finalRef];
    if (me == 0) { continue; }
    if (me != 0 && sign == 0) {
      sign = std::abs(me)/me;
    }
    wth_ref += 1./(multfactor[L] * std::pow(Eg_ref,-power[L]) / (me*me/(2*J_ref+1.)) / std::pow(100,nbarns[L]));
  }

  double wth = std::abs(mat_rel) * wth_ref;
  sign = std::abs(mat_rel)/mat_rel * sign;
  double Eg = std::abs(nucl.GetLevelEnergies()[initialState] - nucl.GetLevelEnergies()[finalState]);
  double J = nucl.GetLevelJ()[initialState];
  double mat = std::sqrt((multfactor[lambda]*std::pow(Eg, -power[lambda]))/((1./wth) * std::pow(100, nbarns[lambda])) * (2.*J + 1.)) * sign;
  nucl.SetMatrixElement(lambda, initialState, finalState, mat);
}

void RelMatWidthElement::Update(const double *res, int &parct) {
  if (GetFixed()) { return; } 
  SetMatRel(res[parct]);
  parct += 1;
}

void RelMatWidthElement::Populate(std::vector<double> &parameters,
                               std::vector<double> &par_LL,
                               std::vector<double> &par_UL) {
    if (!GetFixed()) {
      parameters.push_back(GetMatRel());
      par_LL.push_back(GetMatRelLowerLimit());
      par_UL.push_back(GetMatRelUpperLimit());
    }
}

void RelMatWidthElement::Print(std::ostream &out) const {
  out << type << " " << name << "  :  " << matRel << std::endl;
}

void RelMatWidthElement::Kick(TRandom3 &rand) {
  if (GetFixed()) { return; }
  double newmat;
  if (GetMatRelUpperLimit() * GetMatRelLowerLimit() < 0) {    
    newmat = rand.Gaus(GetMatRel(), std::max(std::abs(GetMatRelLowerLimit()), std::abs(GetMatRelUpperLimit()))/4.);
    double newsign = rand.Uniform() - 0.5;
    if (newsign != 0) {
    newmat = newmat * std::abs(newsign)/newsign;
    }
  }
  else {
    newmat = rand.Gaus(GetMatRel(), (GetMatRelLowerLimit() - GetMatRelUpperLimit())/4.);
  }
  
  if (newmat > GetMatRelUpperLimit()) { newmat = GetMatRelUpperLimit(); }
  if (newmat < GetMatRelLowerLimit()) { newmat = GetMatRelLowerLimit(); }
  
  SetMatRel(newmat);
}
