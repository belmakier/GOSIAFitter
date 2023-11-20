#include "Gosia.h"

void RunGosia(input_file &inputfile,
              std::string dir,
              gdet &all_dets,
              std::vector<double> &matels,
              out_yields &output,
              int verbosity) {
  int dirlen =  (int)dir.size();
  //system(("cd "+dir+"; ./gosia < "+filename+" > "+outfile).c_str());

  double bst_me[999];
  for (int i=0; i<matels.size(); ++i) {
    bst_me[i] = matels[i];
  }

  gosia_(&inputfile, 
         dir.c_str(), dirlen,
         &all_dets,
         &bst_me[0],
         &output,
         verbosity);
}
