#ifndef GOSIAReader_h
#define GOSIAReader_h

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include "Gosia.h"
#include "Nucleus.h"
#include "ExperimentalInput.h"


class GOSIAReader {

	public:

		GOSIAReader(Nucleus*,const char*);
  GOSIAReader(Nucleus*,out_yields &gosia_yields);	
		~GOSIAReader()				{ ;}

		void	Clear()				{ gosiaData.clear();		}
		void	ReadGOSIAFile(const char*);			/*!< Read formatted GOSIA file and extract yields */
  void ReadGOSIAFile(out_yields &gosia_yields);

		std::vector<ExperimentData>		GetGOSIAData()		{ return gosiaData;	}	/*!< Return vector of ExperimentData containing GOSIA yields*/

		void	Print()	const;
  std::vector<double> GetAverageAngles() { return angle_av; }
  double GetSolidAngleWidth(int i) { return std::abs((std::cos(angle_high[i]*3.14159265/180.0) - std::cos(angle_low[i]*3.14159265/180.0))); }
  double GetAngleWidth(int i) { return angle_high[i]*3.14159265/180.0 - angle_low[i]*3.14159265/180.0; }
  std::vector<double> GetRutherfords() { return rutherfords; }

	private:

		Nucleus					fNucleus;	/*!< Nucleus object used to define level assignments */
		std::vector<ExperimentData>		gosiaData;	/*!< Vector of ExperimentData objects to hold the data that is read in from the file */
  std::vector<double> rutherfords;
  std::vector<double> energy_low;
  std::vector<double> energy_high;
  std::vector<double> angle_low;
  std::vector<double> angle_high;
  std::vector<double> angle_av;

};

#endif
