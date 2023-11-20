#include "GOSIASimAnnealing.h"

#include "Math/IFunction.h"
#include "Math/GSLRndmEngines.h"

#include <cassert>
#include <iostream>
#include <cmath>
#include <vector>


GOSIASimAnFunc::GOSIASimAnFunc(const ROOT::Math::IMultiGenFunction & func, const double * x) :
   fX( std::vector<double>(x, x + func.NDim() ) ),
   fScale( std::vector<double>(func.NDim() )),
   fFunc(&func)
{
   // set scale factors to 1
   fScale.assign(fScale.size(), 1.);
}

GOSIASimAnFunc::GOSIASimAnFunc(const ROOT::Math::IMultiGenFunction & func, const double * x, const double * scale) :
   fX( std::vector<double>(x, x + func.NDim() ) ),
   fScale( std::vector<double>(scale, scale + func.NDim() ) ),
   fFunc(&func)
{}


double GOSIASimAnFunc::Energy(const double *x) const {
   // evaluate the energy
   return   (*fFunc)(x);
}

void GOSIASimAnFunc::Step(TRandom3 & random, double *x, double *widths, double scale) {
   // x  -> x + Random[-step,step]     for each coordinate
   unsigned int ndim = NDim();
   for (unsigned int i = 0; i < ndim; ++i) {
     double urndm = random.Uniform();
     double sgnrndm = random.Uniform();     
     double sstep = scale * widths[i] * fScale[i]  ;
     x[i] +=  sstep * std::tan(3.141/2. * (2.0*urndm - 1.0));
     if (sgnrndm < 0.005) {//0.5% chance of sign flip?
       x[i] = -x[i];
     }     
   }
}


double GOSIASimAnFunc::Distance(const GOSIASimAnFunc & f) const {
   // calculate the distance with respect onother configuration
   const std::vector<double> & x = fX;
   const std::vector<double> & y = f.X();
   unsigned int n = x.size();
   assert (n == y.size());
   if (n > 1) {
      double d2 = 0;
      for (unsigned int i = 0; i < n; ++i)
         d2 += ( x[i] - y[i] ) * ( x[i] - y[i] );
      return std::sqrt(d2);
   }
   else
      // avoid doing a sqrt for 1 dim
      return std::abs( x[0] - y[0] );
}

void GOSIASimAnFunc::Print() {
   // print the position  x in standard std::ostream
   // GOSIA prints also niter-  ntrials - temperature and then the energy and energy min value (from 1.10)
   std::cout << "\tx = ( ";
   unsigned n = NDim();
   for (unsigned int i = 0; i < n-1; ++i) {
      std::cout << fX[i] << " , ";
   }
   std::cout << fX.back() << " )\t";
   // energy us printed by GOSIA (and also end-line)
   std::cout << "E  / E_best = ";   // GOSIA print then E and E best
}

GOSIASimAnFunc &  GOSIASimAnFunc::FastCopy(const GOSIASimAnFunc & rhs) {
   // copy only the information which is changed during the process
   // in this case only the x values
   std::copy(rhs.fX.begin(), rhs.fX.end(), fX.begin() );
   return *this;
}



// implementation of GOSIASimAnnealing class


GOSIASimAnnealing::GOSIASimAnnealing()
{
   // Default constructor implementation.
}



// function for solving (from a Genfunction interface)

int GOSIASimAnnealing::Solve(const ROOT::Math::IMultiGenFunction & func, const double * x0, const double * scale, double * xmin, bool debug) {
   // solve the simulated annealing problem given starting point and objective function interface


   // initial conditions
   GOSIASimAnFunc   fx(func, x0, scale);

   int iret =  Solve(fx, debug);

   if (iret == 0) {
      // copy value of the minimum in xmin
      std::copy(fx.X().begin(), fx.X().end(), xmin);
   }
   return iret;

}

double GOSIASimAnnealing::Boltzmann(double E, double new_E, double T) {
  double x = -(new_E - E) / (fParams.k * T);
  return (x < std::log(DBL_MIN)) ? 0.0 : std::exp(x);
}

int GOSIASimAnnealing::Solve(GOSIASimAnFunc & fx, bool debug) {
  int verbosity = fParams.verbosity;
  x = fx.GetXCopy();
  new_x = fx.GetXCopy();
  best_x = fx.GetXCopy();
  grad = fx.GetXCopy();
  widths = fx.GetXCopy();
  
  E = fx.Energy(x.data());
  best_E = E;
  last_E = E;
  nlast_evals = 0;

  if (verbosity > 0) {  
  std::cout << "k = " << fParams.k << " T = " << fParams.t_initial << std::endl;
  std::cout << "max_corr_evals = " << fParams.max_corr_evals << std::endl;
  }
  T = fParams.t_initial;
  T_factor = 1.0/fParams.mu_t;

  double widthscale = 0.5;
  //calculate gradients
  //fx.GetGradient(x.data(),grad.data());
  //make widths array
  double av_width = 0;
  for (int i=0; i<grad.size(); ++i) {
    grad[i] = 1.;
    if (std::abs(grad[i])>1e-1){
      widths[i] = 1.;
      av_width += std::abs(widths[i]);
    }
    else {
      widths[i] = 10;
    }
  }
  av_width = av_width/grad.size();

  int n_cycles = 0;
  while(true) {
    n_accepts = 0;
    n_rejects = 0;
    n_eless = 0;
    ++n_cycles;

    for (int i=0; i<fParams.iters_fixed_T; ++i) {
      std::copy(x.begin(),x.end(), new_x.begin());
      
      //step new_x 
      fx.Step(rand, new_x.data(), widths.data(), widthscale/av_width); //-> pass widths here

      if (verbosity > 0) {
        std::cout << n_cycles << " : " << T << std::endl;
        std::cout << i << "   " << n_evals << "  " << nlast_evals << "   " << last_E << "  " << best_E << std::endl;
      }
      if (n_evals - nlast_evals >= fParams.max_corr_evals || best_E < last_E*fParams.corr_E_drop) {
        reCorr = 1;
        nlast_evals = n_evals;
        //calculate gradients
        //fx.GetGradient(x.data(),grad.data());
        //make widths array
        double av_width = 0;
        for (int i=0; i<grad.size(); ++i) {
          grad[i] = 1.;
          if (std::abs(grad[i])>1e-1){
            widths[i] = 1.;
            av_width += std::abs(widths[i]);
          }
          else {
            widths[i] = 10;
          }
        }
        av_width = av_width/grad.size();

      }
      else {
        reCorr = 0;
      }

      if (reCorr == 1) {
        //need to recalculate correction factors, don't need the step
        std::copy(x.begin(), x.end(), new_x.begin()); //use x for recalculation
        std::copy(x.begin(), x.end(), best_x.begin()); //use x for recalculation
        new_E = fx.Energy(new_x.data()); //recalculation correction factors in here
        last_E = new_E;
        best_E = new_E;
        E = new_E;
        ++n_evals;
      }
      else {
        new_E = fx.Energy(new_x.data());      
        ++n_evals;

        if (new_E < E) {
          if (new_E < best_E) {
            best_E = new_E;
            std::copy(new_x.begin(), new_x.end(), best_x.begin());
          }
        
          E = new_E;
          std::copy(new_x.begin(), new_x.end(), x.begin());
        
          ++n_eless;
        }
        else if (rand.Uniform() < Boltzmann(E, new_E, T)) {
          E = new_E;
          std::copy(new_x.begin(), new_x.end(), x.begin());
        
          ++n_accepts;
        }
        else {
          ++n_rejects;
        }
      }      
    }

    if (verbosity > 0) {
    std::cout << "n_accepts = " << n_accepts << std::endl;
    std::cout << "n_rejects = " << n_rejects << std::endl;
    std::cout << "n_eless = " << n_eless << std::endl;
    }
    double frac = (double)(n_eless + n_accepts)/(double)(n_eless + n_accepts + n_rejects);
    if (verbosity > 0) {
      std::cout << "frac = " << frac << std::endl;
    }
    
    if (frac < 0.4) {
      widthscale *= frac/0.5+0.1;
    }
    else if (frac > 0.6) {
      widthscale *= frac/0.5+0.1;
    }
    widthscale = std::max(widthscale, 0.1);
    if (verbosity > 0) {
    std::cout << "widthscale = " << widthscale << std::endl;
    }
    
    //update scaling factors?

    T *= T_factor;
    ++n_iter;
    if (T < fParams.t_min) {
      break;
    }
  }

  fx.SetX(best_x.data());
  
  return 0;

}

