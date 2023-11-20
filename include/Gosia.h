#ifndef GOSIAFITTER_GOSIA_HH
#define GOSIAFITTER_GOSIA_HH

#include <string>
#include <vector>
#include <iostream>

//Please note that this changes only the C++ side; corresponding definitions in
//GOSIA source file must match this definition
#define INPUTFILE_LINELEN 256

struct out_exp {
  double ruth;
  double en_low;
  double en_high;
  double theta_low;
  double theta_high;
  int number;

  int nyields;
  int ni[999];
  int nf[999];
  double yield[999];
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
