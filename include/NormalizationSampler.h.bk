#ifndef NormalizationSampler_h
#define NormalizationSampler_h

#include "ScalingParameter.h"
#include "Nucleus.h"
//#include "Minuit2/FCNBase.h"
#include "Literature.h"
#include "ExperimentalInput.h"
#include "MatrixElement.h"
#include "TransitionRates.h"
#include "GOSIAReader.h"
#include "ScalingFitFCN.h"

//#include "Minuit2/Minuit2Minimizer.h"
#include "Math/Functor.h"
#include "Math/Factory.h"
#include "Math/Minimizer.h"
//#include "Minuit2/FunctionMinimum.h"
//#include "Minuit2/MnUserParameterState.h"
//#include "Minuit2/MnMigrad.h" 
//#include "Minuit2/MnMinos.h"
//#include "Minuit2/MnContours.h"
//#include "Minuit2/MnPlot.h"

#include <ctime>
#include <thread>
#include <iomanip>
#include <vector>
#include <chrono>
#include <stdio.h>
#include <stdlib.h>

#include "GOSIASimFitter.h"
#include "Nucleus.h"

class	NormalizationSampler : public GOSIASimFitter {

	public:
		void	SetupCalculation();	/*!< Prepare the calculation */
		void	Run(double*);	

		void	ClearNormalizations()					{ Normalizations.clear();	}
		std::vector<std::vector<double>>	GetNormalizations()	{ return Normalizations;	}

	private:
		std::vector<std::vector<double>>	Normalizations;

		std::vector<int>	exptIndex;


};

#endif
