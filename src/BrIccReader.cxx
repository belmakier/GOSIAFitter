#include <iostream>
#include <cmath>
#include <cstdio>
#include <vector>

#include "BrIccReader.h"
#include "MiscFunctions.h"

BrIccReader::BrIccReader(std::string idxpath, std::string iccpath) {
  idxfile = fopen(idxpath.c_str(), "r");
  iccfile = fopen(iccpath.c_str(), "r");
  std::cout << "Conversion electron coefficients from BrIcc database. \nPlease cite T. Kibedi et al. NIM A589 (2009) 202-229 for the conversion coefficients!" << std::endl;
}

void BrIccReader::Open(std::string idxpath, std::string iccpath) {
  idxfile = fopen(idxpath.c_str(), "r");
  iccfile = fopen(iccpath.c_str(), "r");
  std::cout << "Conversion electron coefficients from BrIcc database. \nPlease cite T. Kibedi et al. NIM A589 (2009) 202-229 for the conversion coefficients!" << std::endl;
}

void BrIccReader::Close() {
  if (idxfile) {
    fclose(idxfile);
  }
  if (iccfile) {
    fclose(iccfile);
  }
  idxfile=NULL;
  iccfile=NULL;
}

double BrIccReader::GetTotalCC(int Z, double egamma, int mult) { //note that this mult is defined differently from GOSIAReader mult
  if (!idxfile || !iccfile) { std::cout << "Open BrIcc data base before getting conversion coefficients!" << std::endl; return -1; }
  fseek(idxfile, 2048*(Z-1), SEEK_SET);
  elem e;
  fread(&e, 12, 1, idxfile);

  if (egamma > 6000) { return 0; }

  std::vector<BrIccReader::eshell> shell;
  for (int i=0; i<37; ++i) {
    eshell rec;
    fread(&rec, 32, 1, idxfile);
    if ((int)rec.exist > 0) { rec.exist = true; }
    else { rec.exist = false; }
    shell.push_back(rec);
  }

  if (egamma < 1100) { 
    shell[36].exist = false; //no IPC below 1100 keV
  }
  if ((mult > 2 && mult < 5) || (mult > 7 && mult < 10)) {
    shell[36].exist = false; //no IPC for E4,E5,M4,M5
  }

  float energy[37][200];
  float cc[37][200][10];
  for (int i = 0; i < 37; ++i) {
    if (egamma < shell[i].be) { shell[i].exist = false; }
    if (!shell[i].exist) { continue; }
    int iEn = 0;
    for (int rec = shell[i].nrec - 1; rec < shell[i].nrec + shell[i].nmesh - 1; ++rec) {
      fseek(iccfile, rec*44, SEEK_SET);
      icc coef;
      fread(&coef, 44, 1, iccfile);
      energy[i][iEn] = coef.energy;
      if (iEn > 200) { std::cout << "SEVERE ERROR! " << iEn << std::endl; }
      for (int imult=0; imult<10; ++imult) { cc[i][iEn][imult] = coef.icc[imult]; }
      ++iEn;
    }        
  }

  double total_cc = 0;
  for (int i=0; i<37; ++i) {
    double this_cc = 0;
    if (!shell[i].exist) { continue; }
    if (egamma < energy[i][0]) {
      this_cc = cc[i][0][mult];
      std::cout << "Warning! Egamma = " << egamma << " is in the regime where solid state effects dominate, conversion coefficients for shell ";
      for (int j=0; j<8; ++j) { std::cout << shell[i].shell[j]; }
      std::cout << std::endl;      
    }
    else if (egamma > energy[i][shell[i].nmesh-1]) {
      this_cc = cc[i][shell[i].nmesh-1][mult];
      std::cout << "Warning! Egamma = " << egamma << " exceeds the range of conversion coefficients table for shell ";
      for (int j=0; j<8; ++j) { std::cout << shell[i].shell[j]; }
      std::cout << std::endl;      
    }
    else {
      int n = shell[i].nmesh;
      double *x = new double[n]();
      double *y = new double[n]();

      for (int iEn = 0; iEn < n; ++iEn) {
        x[iEn] = std::log(energy[i][iEn]);
        y[iEn] = cc[i][iEn][mult];
      }
      MiscFunctions::Spline(x,y,n,std::log(egamma),this_cc);

      delete[] x;
      delete[] y;
    }    
    total_cc += this_cc;
  }

  return total_cc;
 
}


double BrIccReader::GetTotalOmg(int Z, double egamma) { //Gets Omega for E0 transitions
  if (!idxfile || !iccfile) { std::cout << "Open BrIcc data base before getting conversion coefficients!" << std::endl; return -1; }
  fseek(idxfile, 2048*(Z-1), SEEK_SET);
  elem e;
  fread(&e, 12, 1, idxfile);

  if (egamma > 6000) { return 0; }

  std::vector<BrIccReader::eshell> shell;
  for (int i=0; i<41; ++i) {
    eshell rec;
    fread(&rec, 32, 1, idxfile);
    if ((int)rec.exist > 0) { rec.exist = true; }
    else { rec.exist = false; }
    shell.push_back(rec);
  }

  if (egamma < 1100) { 
    shell[40].exist = false; //no IPC below 1100 keV
  }

  float energy[4][200];
  float omg[4][200];
  for (int i = 37; i < 41; ++i) {
    if (egamma < shell[i].be) { shell[i].exist = false; }
    if (!shell[i].exist) { continue; }
    int iEn = 0;
    for (int rec = shell[i].nrec - 1; rec < shell[i].nrec + shell[i].nmesh - 1; ++rec) {
      fseek(iccfile, rec*44, SEEK_SET);
      icc coef;
      fread(&coef, 44, 1, iccfile);
      energy[i-37][iEn] = coef.energy;
      if (iEn > 200) { std::cout << "SEVERE ERROR! " << iEn << std::endl; }
      omg[i-37][iEn] = coef.icc[0];
      ++iEn;
    }        
  }

  double total_omg = 0;
  for (int i=37; i<41; ++i) {
    double this_omg = 0;
    if (!shell[i].exist) { continue; }
    if (egamma < energy[i-37][0]) {
      this_omg = omg[i-37][0];
      std::cout << "Warning! Egamma = " << egamma << " is in the regime where solid state effects dominate, conversion coefficients for shell ";
      for (int j=0; j<8; ++j) { std::cout << shell[i].shell[j]; }
      std::cout << std::endl;      
    }
    else if (egamma > energy[i-37][shell[i].nmesh-1]) {
      this_omg = omg[i-37][shell[i].nmesh-1];
      std::cout << "Warning! Egamma = " << egamma << " exceeds the range of conversion coefficients table for shell ";
      for (int j=0; j<8; ++j) { std::cout << shell[i].shell[j]; }
      std::cout << std::endl;      
    }
    else {
      int n = shell[i].nmesh;
      double *x = new double[n]();
      double *y = new double[n]();

      for (int iEn = 0; iEn < n; ++iEn) {
        x[iEn] = std::log(energy[i-37][iEn]);
        y[iEn] = omg[i-37][iEn];
      }
      MiscFunctions::Spline(x,y,n,std::log(egamma),this_omg);

      delete[] x;
      delete[] y;
    }    
    total_omg += this_omg;
  }

  return total_omg;
 
}
