#ifndef MatrixElement_h
#define MatrixElement_h

#include <iostream>
#include <iomanip>

#include "Nucleus.h"

#include "TRandom3.h"

class FittingElement {
public:
  FittingElement() {;}
  FittingElement(std::string t, std::string n,
                 bool f, int np) {
    type = t;
    name = n;
    fixed = f;
    npars = np;
  }
  virtual ~FittingElement() {};

  virtual void Propagate(Nucleus &nucl,
                         const double *par, int &parct, int opt) = 0; //this implements the parameters from *this or par into the nucleus
  virtual void Populate(std::vector<double> &parameters,
                        std::vector<double> &par_LL,
                        std::vector<double> &par_UL) = 0;  //this puts the parameters from *this into the vectors (to then pass to the TMinuit2)
  virtual void Update(const double *res, int &parct) = 0; //this updates the values of *this from those contained in res
  virtual void Kick(TRandom3 &rand) = 0; //randomizes/perturbs somehow
  virtual void Print(std::ostream &out) const = 0;
  bool GetFixed() const { return fixed; }
  void SetFixed(bool f) { fixed = f; }
  int GetIndex() const { return indx; }
  void SetIndex(int i) { indx = i; }

  std::string GetType() const { return type; }
  std::string GetName() const { return name; }
  int GetNPars() const { return npars; }
protected:
  int indx;
  bool fixed;
  std::string type;
  std::string name; // unique identifier
  int npars;
};

///
///	\class MatrixElement
///
///	\brief Holder class for matrix element information for use in fitting routines
///

class MatrixElement : public FittingElement {

public:
  MatrixElement()	{;}
  MatrixElement (std::string name, int i, int l, int s1, int s2, double me, double mell, double meul, bool fx=false)
    : FittingElement("ME", name, false,1) {
									 	index 			= i; 
										lambda 			= l; 
										initialstate 		= s1; 
										finalstate 		= s2; 
										matrixElement 		= me;	
										matrixElement_ll 	= mell;
										matrixElement_ul 	= meul;
                    fixed=fx;
		}							/*!< Construct matrix element me, of index i, multipolarity l, between states s1 and s2 with lower and upper limits of mell and meul  */
		~MatrixElement() {;}
		MatrixElement& operator = (const MatrixElement& m);	/*!< Assignment operator */
		MatrixElement(const MatrixElement& m);			/*!< Copy constructor */

		void	SetupME(int i, int l, int s1, int s2, double me, double mell, double meul){
									 	index 			= i; 
										lambda 			= l; 
										initialstate 		= s1; 
										finalstate 		= s2; 
										matrixElement 		= me;	
										matrixElement_ll 	= mell;
										matrixElement_ul 	= meul;
		}							/*!< Define matrix element me, of index i, multipolarity l, between states s1 and s2 with lower and upper limits of mell and meul  */
  
		void	SetMatrixElement(double ME)				{ matrixElement = ME;		}	/*!< Define matrix element value */
		void	SetMatrixElementLowerLimit(double ME_LL)		{ matrixElement_ll = ME_LL;	}	/*!< Define matrix element lower limit */
		void	SetMatrixElementUpperLimit(double ME_UL)		{ matrixElement_ul = ME_UL;	}	/*!< Define matrix element upper limit */
    void SetFixed(bool f) { fixed = f; }
		double 	GetMatrixElement() const				{ return matrixElement;		}	/*!< Return matrix element value */
		double	GetMatrixElementLowerLimit() const			{ return matrixElement_ll;	}	/*!< Return matrix element lower limit */
		double 	GetMatrixElementUpperLimit() const			{ return matrixElement_ul;	}	/*!< Return matrix element upper limit */
		int	GetIndex() const					{ return index;			}	/*!< Return matrix element index */
		int 	GetLambda() const 					{ return lambda;		}	/*!< Return matrix element mulitpolarity */
		int 	GetInitialState() const					{ return initialstate;		}	/*!< Return initial state index */
		int	GetFinalState()	const					{ return finalstate;		}	/*!< Return final state index */


  void	Print(std::ostream &out) const;		/*!< Print matrix element information */
  void Kick(TRandom3 &rand);

  void Propagate(Nucleus &nucl,
                 const double *par, int &parct, int opt=0);
  void Update(const double *res, int &parct);

  void Populate(std::vector<double> &parameters,
                        std::vector<double> &par_LL,
                          std::vector<double> &par_UL);

	private:
		int 		index;			/*!< Matrix element index */
		int 		lambda;			/*!< Matrix element multipolarity */
		int 		initialstate;		/*!< Initial state (initial and final are arbitrary) */
		int 		finalstate;		/*!< Final state (initial and final are arbitrary)*/
		double		matrixElement;		/*!< Matrix element value */
		double		matrixElement_ll;	/*!< Matrix element lower limit */
		double		matrixElement_ul;	/*!< Matrix element upper limit */

};

///
///	\class RelativeMatrixElement
///
///	\brief Holder class for matrix element information for use in fitting routines, where this matrix element is fixed relative to some other matrix element
///


class RelativeMatrixElement : public FittingElement {

	public:
		RelativeMatrixElement()	{;}
  RelativeMatrixElement(std::string name, int i, int l, int s1, int s2, int l2, int s12, int s22, double r, double rll, double rul, bool fx=false, bool at=false) :
    FittingElement("RelME",name, fx,1) {
      index 			= i; 
      lambda 			= l; 
      initialstate 		= s1; 
      finalstate 		= s2;
      lambdaRel 			= l2; 
      initialstateRel 		= s12; 
      finalstateRel 		= s22; 
      rel = r;
      relll = rll;
      relul = rul;
      fixed = fx;
      arctan = at;
		}							/*!< Construct matrix element me, of index i, multipolarity l, between states s1 and s2 with lower and upper limits of mell and meul  */
		~RelativeMatrixElement() {;}
		RelativeMatrixElement& operator = (const RelativeMatrixElement& m);	/*!< Assignment operator */
		RelativeMatrixElement(const RelativeMatrixElement& m);			/*!< Copy constructor */

  void	SetupME(int i, int l, int s1, int s2, int l2, int s12, int s22, double r, double rll, double rul, bool fx=false, bool at = false){
      index 			= i; 
      lambda 			= l; 
      initialstate 		= s1; 
      finalstate 		= s2;
      lambdaRel 			= l2; 
      initialstateRel 		= s12; 
      finalstateRel 		= s22; 
      rel = r;
      relll = rll;
      relul = rul;
      fixed = fx;
      arctan = at;
		}							/*!< Define matrix element me, of index i, multipolarity l, between states s1 and s2 with lower and upper limits of mell and meul  */
		void	SetRelativeElement(double r)				{ rel = r;		}	/*!< Define matrix element value */
    void	SetRelativeElementLowerLimit(double ME_LL)		{ relll = ME_LL;	}	/*!< Define matrix element lower limit */
		void	SetRelativeElementUpperLimit(double ME_UL)		{ relul = ME_UL;	}	/*!< Define matrix element upper limit */

		double 	GetRelativeElement() const				{ return rel;		}	/*!< Return matrix element value */
    double 	GetRelativeElementLowerLimit() const				{ return relll;		}	/*!< Return matrix element value */
    double 	GetRelativeElementUpperLimit() const				{ return relul;		}	/*!< Return matrix element value */
		int	    GetIndex() const					{ return index;			}	/*!< Return matrix element index */
		int 	  GetLambda() const 					{ return lambda;		}	/*!< Return matrix element mulitpolarity */
		int 	  GetInitialState() const					{ return initialstate;		}	/*!< Return initial state index */
		int	    GetFinalState()	const					{ return finalstate;		}	/*!< Return final state index */
  	int 	  GetLambdaRel() const 					{ return lambdaRel;		}	/*!< Return matrix element mulitpolarity */
		int 	  GetInitialStateRel() const					{ return initialstateRel;		}	/*!< Return initial state index */
		int	    GetFinalStateRel()	const					{ return finalstateRel;		}	/*!< Return final state index */
    void    SetFixed(bool fx) { fixed = fx; }
    bool    GetFixed() const { return fixed; };
  bool GetArcTan() const { return arctan; };
  void	  Print(std::ostream &out) const;		/*!< Print matrix element information */
  void Kick(TRandom3 &rand);
    void Propagate(Nucleus &nucl,
                   const double *par, int &parct, int opt = 0);
    void Update(const double *res, int &parct);
  void Populate(std::vector<double> &parameters,
                        std::vector<double> &par_LL,
                std::vector<double> &par_UL);

	private:
		int 		  index;			/*!< Matrix element index */
		int 		  lambda;			/*!< Matrix element multipolarity */
		int 		  initialstate;		/*!< Initial state (initial and final are arbitrary) */
		int 		  finalstate;		/*!< Final state (initial and final are arbitrary)*/
    int 		  lambdaRel;			/*!< Matrix element multipolarity */
		int 		  initialstateRel;		/*!< Initial state (initial and final are arbitrary) */
		int 		  finalstateRel;		/*!< Final state (initial and final are arbitrary)*/
		double		rel;		/*!< Matrix element value */
    double		relll;		/*!< Matrix element value */
    double		relul;		/*!< Matrix element value */
  bool arctan;

};

class LifetimeMixingElement : public FittingElement {
public:
  LifetimeMixingElement() {}
  LifetimeMixingElement(std::string name, int init, int fina, int l1, int l2,
                        double wth, double wth_ll, double wth_ul,
                        double mix, double mix_ll, double mix_ul,
                        bool fx = false) :
    FittingElement("LtMixE", name, fx, 2)
  {
    initialState = init;
    finalState = fina;
    lambda1 = l1;
    lambda2 = l2;
    width = wth;
    width_ll = wth_ll;
    width_ul = wth_ul;
    mixingRatio = mix;
    mixingRatio_ll = mix_ll;
    mixingRatio_ul = mix_ul;
    fixed = false;          
  }
  ~LifetimeMixingElement() {}

  
  void	  Print(std::ostream &out) const;		/*!< Print matrix element information */
  void Kick(TRandom3 &rand);
  void Propagate(Nucleus &nucl,
                 const double *par, int &parct, int opt = 0);
  void Update(const double *res, int &parct);
  void Populate(std::vector<double> &parameters,
                std::vector<double> &par_LL,
                std::vector<double> &par_UL);

  double GetWidth() { return width; }
  void SetWidth(double wth) { width = wth; }
  double GetWidthUpperLimit() { return width_ul; }
  double GetWidthLowerLimit() { return width_ll; }
  
  double GetMixingRatio() { return mixingRatio; }
  void SetMixingRatio(double mr) { mixingRatio = mr; }
  double GetMixingRatioUpperLimit() { return mixingRatio_ul; }
  double GetMixingRatioLowerLimit() { return mixingRatio_ll; }

private:
  int initialState;
  int finalState;
  int lambda1;
  int lambda2;
  double width;
  double width_ll;
  double width_ul;
  double mixingRatio;
  double mixingRatio_ll;
  double mixingRatio_ul;

  std::vector<int> nbarns = {1,2,3,4,5,6,0,1};
  std::vector<double> multfactor = {0.62887e-15, 816.24e-12, 1752.1e-6, 5894.45, 2.89e10, 1.95e17, 56.870e-15, 0.738144e-7};
  std::vector<int> power = {3,5,7,9,11,13,3,5};

};

class RelLtMixElement : public FittingElement {
public:
  RelLtMixElement() {}
  RelLtMixElement(std::string name, int init, int fina, int l1, int l2,
                        double rel_wth, double rel_wth_ll, double rel_wth_ul,
                        double mix, double mix_ll, double mix_ul,
                  int ref_init, int ref_fina,                  
                  bool fx = false) :
    FittingElement("RelLtMx", name, fx, 2) {
    initialState = init;
    initialState = init;
    finalState = fina;
    lambda1 = l1;
    lambda2 = l2;
    widthRel = rel_wth;
    widthRel_ll = rel_wth_ll;
    widthRel_ul = rel_wth_ul;
    mixingRatio = mix;
    mixingRatio_ll = mix_ll;
    mixingRatio_ul = mix_ul;
    initRef = ref_init;
    finalRef = ref_fina;
    fixed = fx;
  }

  void	  Print(std::ostream &out) const;		/*!< Print matrix element information */
  void Kick(TRandom3 &rand);
  void Propagate(Nucleus &nucl,
                 const double *par, int &parct, int opt = 0);
  void Update(const double *res, int &parct);
  void Populate(std::vector<double> &parameters,
                std::vector<double> &par_LL,
                std::vector<double> &par_UL);

  double GetWidthRel() { return widthRel; }
  void SetWidthRel(double wth) { widthRel = wth; }
  double GetWidthRelUpperLimit() { return widthRel_ul; }
  double GetWidthRelLowerLimit() { return widthRel_ll; }

  double GetMixingRatio() { return mixingRatio; }
  void SetMixingRatio(double mr) { mixingRatio = mr; }
  double GetMixingRatioUpperLimit() { return mixingRatio_ul; }
  double GetMixingRatioLowerLimit() { return mixingRatio_ll; }

private:
  int initialState;
  int finalState;
  int lambda1;
  int lambda2;
  double widthRel;
  double widthRel_ll;
  double widthRel_ul;
  double mixingRatio;
  double mixingRatio_ll;
  double mixingRatio_ul;
  int initRef;
  int finalRef;

  std::vector<int> nbarns = {1,2,3,4,5,6,0,1};
  std::vector<double> multfactor = {0.62887e-15, 816.24e-12, 1752.1e-6, 5894.45, 2.89e10, 1.95e17, 56.870e-15, 0.738144e-7};
  std::vector<int> power = {3,5,7,9,11,13,3,5};
};

class RelMatWidthElement : public FittingElement {
public:
  RelMatWidthElement() {}
  ~RelMatWidthElement() {}
  RelMatWidthElement(std::string name, int l, int init, int fina,
                     double relmat, double relmat_ll, double relmat_ul,
                     int init_ref, int final_ref,
                     bool fx=false) :
    FittingElement("RelMatW", name, fx, 1) {
    lambda = l;
    initialState = init;
    finalState = fina;
    matRel = relmat;
    matRel_ll = relmat_ll;
    matRel_ul = relmat_ul;
    initRef = init_ref;
    finalRef = final_ref;        
  }

  void	  Print(std::ostream &out) const;		/*!< Print matrix element information */
  void Kick(TRandom3 &rand);
  void Propagate(Nucleus &nucl,
                 const double *par, int &parct, int opt = 0);
  void Update(const double *res, int &parct);
  void Populate(std::vector<double> &parameters,
                std::vector<double> &par_LL,
                std::vector<double> &par_UL);

  double GetMatRel() { return matRel; }
  double GetMatRelUpperLimit() { return matRel_ul; }
  double GetMatRelLowerLimit() { return matRel_ll; }
  void SetMatRel(double mr) { matRel = mr; }

private:
  int initialState;
  int finalState;
  int lambda;
  double matRel;
  double matRel_ll;
  double matRel_ul;
  int initRef;
  int finalRef;

  std::vector<int> nbarns = {1,2,3,4,5,6,0,1};
  std::vector<double> multfactor = {0.62887e-15, 816.24e-12, 1752.1e-6, 5894.45, 2.89e10, 1.95e17, 56.870e-15, 0.738144e-7};
  std::vector<int> power = {3,5,7,9,11,13,3,5};

};

#endif
