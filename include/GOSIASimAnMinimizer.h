#ifndef GOSIASimAnMinimizer_h
#define GOSIASimAnMinimizer_h

#include "GOSIASimAnnealing.h"

#include "Math/BasicMinimizer.h"

#include "Math/IFunctionfwd.h"

#include "Math/IParamFunctionfwd.h"

class GOSIASimAnMinimizer : public ROOT::Math::BasicMinimizer {
public:
  GOSIASimAnMinimizer(int t = 0);

  ~GOSIASimAnMinimizer() override;

private:
  // usually copying is non trivial, so we make this unaccessible

  GOSIASimAnMinimizer(const GOSIASimAnMinimizer &) : ROOT::Math::BasicMinimizer() {}

  GOSIASimAnMinimizer &operator=(const GOSIASimAnMinimizer &rhs)
  {
    if (this == &rhs)
      return *this; // time saving self-test
    return *this;
  }

public:
  /// method to perform the minimization
  bool Minimize() override;

  /// number of calls
  unsigned int NCalls() const override;

  /// Get current minimizer option parameteres
  const GOSIASimAnParams &MinimizerParameters() const { return fSolver.Params(); }

  /// set new minimizer option parameters using directly the GOSIASimAnParams structure
  void SetParameters(const GOSIASimAnParams &params)
  {
    fSolver.SetParams(params);
    DoSetMinimOptions(params); // store new parameters also in MinimizerOptions
  }
  int GetRecorr() { return fSolver.GetRecorr(); } //get whether we need to recalculat the correction factors

protected:
  /// set minimizer option parameters from stored ROOT::Math::MinimizerOptions (fOpt)
  void DoSetSimAnParameters(const ROOT::Math::MinimizerOptions &opt);

  /// Set the Minimizer options from the simulated annealing parameters
  void DoSetMinimOptions(const GOSIASimAnParams &params);  

private:
  GOSIASimAnnealing fSolver;


};

#endif 
