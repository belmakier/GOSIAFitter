#ifndef GOSIAFITTER_GOSIA_HH
#define GOSIAFITTER_GOSIA_HH

#include <string>
#include <vector>
#include <iostream>

//Please note that this changes only the C++ side; corresponding definitions in
//GOSIA source file must match this definition
#define INPUTFILE_LINELEN 512

/*
struct out_det {
  int number;
  int nyields;
  int *ni;
  int *nf;
  double *yield;

  out_det() {
    ni = new int[256];
    nf = new int[256];
    yield = new double[256];
  }
  ~out_det() {
    delete ni;
    delete nf;
    delete yield;
  }
};
*/

struct out_exp {
  int number;

  //these are actually what GOSIA calls "clusters"
  int ndets;
  
  //these are total across all detectors
  int nyields;
  int ni[999];
  int nf[999];

  double ruth;
  double en_low;
  double en_high;
  double theta_low;
  double theta_high;

  double yield[999];


  //out_det dets[32];
};

struct out_yields {
  int nexp;
  out_exp experiment[50];
};

struct bricc_lookup {
  int initialized = 0;
  int nenergies;
  double energy[1500];
  double icc[5][1500];
};

struct detector {
  float radius;
  float energy0;
  float C1[8];
  float C2[8];
  float Qk[8];
};

struct gdet {
  int ndets;
  detector det[32];
};

struct input_file {
  int nlines;
  char lines[1000][INPUTFILE_LINELEN];

  input_file() { 
    for (int i=0; i<1000; ++i) {
      for (int j=0; j<INPUTFILE_LINELEN; ++j) {
        lines[i][j] = (char)0;
      }
    }
  }
};

extern "C" {
  void gosia_(input_file *infile, 
              const char dir[], int &dirlen,
              gdet *all_dets,
              double *bst_me,
              out_yields *output,
              int &vflag);
}

void RunGosia(input_file &inputfile,
              std::string dir,
              gdet &all_dets,
              std::vector<double> &matels,
              out_yields &output,
              int verbosity);

#endif
