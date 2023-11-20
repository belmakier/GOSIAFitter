#include "GOSIASimAnMinimizer.h"
#include "Math/WrappedParamFunction.h"
#include "Math/Error.h"

#include "Math/MinimTransformFunction.h"
#include "Math/MultiNumGradFunction.h"   // needed to use transformation function
#include "Math/FitMethodFunction.h"
#include "Math/GenAlgoOptions.h"

#include <iostream>
#include <cassert>

GOSIASimAnMinimizer::GOSIASimAnMinimizer(int t) :
  ROOT::Math::BasicMinimizer()
{
   // Constructor implementation : create GSLMultiFit wrapper object

   fOptions.SetPrintLevel(0);
   fOptions.SetMinimizerType("GOSIASimAn");
   // set dummy values since those are not used
   fOptions.SetTolerance(-1);
   fOptions.SetMaxIterations(-1);
   fOptions.SetMaxFunctionCalls(0);
   fOptions.SetStrategy(-1);
   fOptions.SetErrorDef(0);
   fOptions.SetPrecision(0);
   fOptions.SetMinimizerAlgorithm("");

   // set extra options
   DoSetMinimOptions(MinimizerParameters());
}

GOSIASimAnMinimizer::~GOSIASimAnMinimizer () {
}


bool GOSIASimAnMinimizer::Minimize() {
   // set initial parameters of the minimizer
   int debugLevel = PrintLevel();

   if (debugLevel >= 1)
      std::cout << "Minimize using GOSIASimAnMinimizer " << std::endl;


   const ROOT::Math::IMultiGenFunction * function = ObjFunction();
   if (function == 0) {
      MATH_ERROR_MSG("GOSIASimAnMinimizer::Minimize","Function has not been set");
      return false;
   }

   // set Sim. An. parameters from existing minimizer options
   DoSetSimAnParameters(fOptions);
   if (debugLevel >= 1) {
      std::cout << "Parameters for simulated annealing: " << std::endl;
      auto simanOpt = fOptions.ExtraOptions();
      if (simanOpt)
         simanOpt->Print();
      else
         std::cout << "no simulated annealing options available" << std::endl;
   }

   // vector of internal values (copied by default)
   unsigned int npar = NPar();
   std::vector<double> xvar;
   std::vector<double> steps(StepSizes(),StepSizes()+npar);

   // needed for the transformation
   std::unique_ptr<ROOT::Math::MultiNumGradFunction> gradFunc(new ROOT::Math::MultiNumGradFunction( *function ));
   //gradFunc->SetOwnership();

   std::unique_ptr<ROOT::Math::MinimTransformFunction>  trFunc(CreateTransformation(xvar, gradFunc.get() ));

   if (trFunc) {
      // transform also the step sizes
      trFunc->InvStepTransformation(X(), StepSizes(), &steps[0]);
      steps.resize( trFunc->NDim() );
   }

   assert (xvar.size() == steps.size() );

   // output vector
   std::vector<double> xmin(xvar.size() );


   int iret = fSolver.Solve( (trFunc) ? *trFunc : *function, xvar.data(), steps.data(), xmin.data(), (debugLevel > 1) );

   SetFinalValues(xmin.data(), trFunc.get());
   SetMinValue( (*ObjFunction())( X() ) );


   if (debugLevel >=1 ) {
      if (iret == 0)
         std::cout << "GOSIASimAnMinimizer: Minimum Found" << std::endl;
      else
         std::cout << "GOSIASimAnMinimizer: Error in solving" << std::endl;

      int pr = std::cout.precision(18);
      std::cout << "FVAL         = " << MinValue() << std::endl;
      std::cout.precision(pr);
      for (unsigned int i = 0; i < NDim(); ++i)
         std::cout << VariableName(i) << "\t  = " << X()[i] << std::endl;
   }


   return ( iret == 0) ? true : false;
}


unsigned int GOSIASimAnMinimizer::NCalls() const {
   // return number of function calls
  auto f = dynamic_cast<const ROOT::Math::FitMethodFunction*>(ObjFunction());
   if (f) return f->NCalls();
   auto gf = dynamic_cast<const ROOT::Math::FitMethodGradFunction*>(ObjFunction());
   if (gf) return gf->NCalls();
   return 0;
}


void GOSIASimAnMinimizer::DoSetMinimOptions(const GOSIASimAnParams &params)
{
   // set the extra minimizer options from GOSIASimAnParams
   ROOT::Math::GenAlgoOptions simanOpt;
   simanOpt.SetValue("n_tries",params.n_tries);
   simanOpt.SetValue("iters_fixed_T",params.iters_fixed_T);
   simanOpt.SetValue("step_size",params.step_size);
   simanOpt.SetValue("k",params.k);
   simanOpt.SetValue("t_initial",params.t_initial);
   simanOpt.SetValue("mu_t",params.mu_t);
   simanOpt.SetValue("t_min",params.t_min);
   simanOpt.SetValue("max_corr_evals",params.max_corr_evals);
   simanOpt.SetValue("corr_E_drop",params.corr_E_drop);
   simanOpt.SetValue("verbosity",params.verbosity);

   fOptions.SetExtraOptions(simanOpt);
}


void GOSIASimAnMinimizer::DoSetSimAnParameters(const ROOT::Math::MinimizerOptions & options)
{
   // get the specific simulated annealing options from MinimizerOptions
   const ROOT::Math::IOptions *simanOpt = options.ExtraOptions();

   if (!simanOpt) {
      return;
   }

   // get the various options. In case they are not defined the default parameters
   // will not be modified and will be used
   GOSIASimAnParams params;
   simanOpt->GetValue("n_tries", params.n_tries);
   simanOpt->GetValue("iters_fixed_T", params.iters_fixed_T);
   simanOpt->GetValue("step_size", params.step_size);
   simanOpt->GetValue("k", params.k);
   simanOpt->GetValue("t_initial", params.t_initial);
   simanOpt->GetValue("mu_t", params.mu_t);
   simanOpt->GetValue("t_min", params.t_min);
   simanOpt->GetValue("max_corr_evals", params.max_corr_evals);
   simanOpt->GetValue("corr_E_drop", params.corr_E_drop);
   simanOpt->GetValue("verbosity", params.verbosity);

   SetParameters(params);
}
