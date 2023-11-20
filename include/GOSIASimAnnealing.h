#ifndef GOSIASimAnnealing_h
#define GOSIASimAnnealing_h

#include "Math/IFunctionfwd.h"
#include "Math/IFunction.h"

#include "TRandom.h"
#include "TRandom3.h"

#include <vector>

class GOSIASimAnFunc {
public:

  GOSIASimAnFunc(const ROOT::Math::IMultiGenFunction & func, const double * x);

  GOSIASimAnFunc(const ROOT::Math::IMultiGenFunction & func, const double * x, const double * scale);

protected:

  GOSIASimAnFunc() :
    fFunc(nullptr)
  {}

public:


  /// virtual destructor (no operations)
  virtual ~GOSIASimAnFunc() { } //


  virtual GOSIASimAnFunc & FastCopy(const GOSIASimAnFunc & f);


  virtual GOSIASimAnFunc * Clone() const {
    return new GOSIASimAnFunc(*this);
  }

  virtual double Energy(const double *x) const;

  virtual void Step(TRandom3 & r, double *x, double *widths, double scale);

  virtual double Distance(const GOSIASimAnFunc & func) const;

  virtual void Print();

  void SetX(const double * x) {
    std::copy(x, x+ fX.size(), fX.begin() );
  }

  template <class IT>
  void SetX(IT begin, IT end) {
    std::copy(begin, end, fX.begin() );
  }

  unsigned int NDim() const { return fX.size(); }

  double X(unsigned int i) const { return fX[i]; }

  const std::vector<double> &  X() const { return fX; }
  
  std::vector<double>  GetXCopy() { return fX; }

  double Scale(unsigned int i) const { return fScale[i]; }

  void SetX(unsigned int i, double x) { fX[i] = x; }

  void GetGradient(double *x, double *g) {
    const ROOT::Math::IMultiGradFunction *gradfunc = dynamic_cast<const ROOT::Math::IMultiGradFunction *>(fFunc);
    gradfunc->Gradient(x,g);
  }

  // use compiler generated  copy ctror and assignment operators

private:

  std::vector<double>  fX;
  std::vector<double>  fScale;
  const ROOT::Math::IMultiGenFunction * fFunc;

};

struct GOSIASimAnParams {

  // constructor with some default values
  GOSIASimAnParams() {
    n_tries =    200;
    iters_fixed_T =  10;
    step_size =   10;
    // the following parameters are for the Boltzmann distribution */
    k = 1.0;
    t_initial =  0.002;
    mu_t =  1.005;
    t_min = 2.0E-6;
    max_corr_evals = 200;
    corr_E_drop = 0.5;
    verbosity = 1;
  }


  int n_tries;            // number of points to try for each step
  int iters_fixed_T;      // number of iterations at each temperature
  double step_size;       // max step size used in random walk
  /// parameters for the Boltzman distribution
  double k;
  double t_initial;
  double mu_t;
  double t_min;
  int max_corr_evals;
  double corr_E_drop;
  int verbosity;
};

class GOSIASimAnnealing {
public:
  GOSIASimAnnealing();
  ~GOSIASimAnnealing() {}

private:
  // usually copying is non trivial, so we make this unaccessible

  GOSIASimAnnealing(const GOSIASimAnnealing &) {}

  GOSIASimAnnealing& operator=(const GOSIASimAnnealing & rhs)  {
    if (this == &rhs) return *this;  // time saving self-test
    return *this;
  }

public:

  int Solve(const ROOT::Math::IMultiGenFunction & func, const double * x0, const double * scale, double * xmin, bool debug = false);

  int Solve(GOSIASimAnFunc & func, bool debug = false);


  GOSIASimAnParams & Params() { return fParams; }
  const GOSIASimAnParams & Params() const { return fParams; }
  void SetParams(const GOSIASimAnParams & params) { fParams = params; }
  int GetRecorr() { return reCorr; }
  double Boltzmann(double E, double new_E, double T);

protected:


private:
  int ndim;
  
  GOSIASimAnParams fParams; // parameters for GOSIASimAnnealing

  std::vector<double>x, new_x, best_x, grad, widths;
  
  double E;
  double new_E;
  double best_E;
  double T;
  double T_factor;
  int n_evals = 1;
  int n_iter = 0;
  int n_accepts, n_rejects, n_eless;
  TRandom3 rand;

  int nlast_evals;
  double last_E;

  int reCorr;

};


#endif 
