#include "MiscFunctions.h"

double MiscFunctions::c 	= 299792458; 
double MiscFunctions::hbar	= 1.05457e-34;

double MiscFunctions::GetMaxMatrix(TMatrixD mat){

	double max = 0;
	for(int i=0;i<mat.GetNcols();i++){
		for(int j=0;j<mat.GetNrows();j++){
			if(mat[j][i] > max) max = mat[j][i];
		}
	}
	return max;

}

double MiscFunctions::GetMaxAbsMatrix(TMatrixD mat){

	double max = 0;
	for(int i=0;i<mat.GetNcols();i++){
		for(int j=0;j<mat.GetNrows();j++){
			if(mat[j][i] > TMath::Abs(max)) max = TMath::Abs(mat[j][i]);
		}
	}
	return max;

}

void MiscFunctions::WriteMatrixNucleus(std::ofstream& outfile, TMatrixD tmpMat, Nucleus nucl){

	outfile		<< std::setw(12) << std::left << "Index"
			<< std::setw(12) << " ";
	for(int s = 0; s < nucl.GetNstates(); s++)
		outfile 	<< std::setw(12) << std::left << s;
	outfile		<< std::endl;	

	outfile 	<< std::setw(12) << std::left << " "
			<< std::setw(12) << std::left << "State";
	for(int s = 0; s < nucl.GetNstates(); s++)
		outfile		<< std::setw(12) << std::left << nucl.GetLevelJ().at(s);
	outfile		<< std::endl;

	for(int y = 0; y < tmpMat.GetNrows(); y++){
		outfile 	<< std::setw(12) << std::left << y
				<< std::setw(12) << std::left << nucl.GetLevelJ().at(y);
		for(int x = 0; x < tmpMat.GetNcols(); x++){
			outfile		<< std::setw(12) << std::left << tmpMat[y][x];
		}
		outfile		<< std::endl;
	}

}

void MiscFunctions::PrintMatrixNucleus(TMatrixD tmpMat, Nucleus nucl){

	std::cout	<< std::setw(12) << std::left << "Index"
			<< std::setw(12) << " ";
	for(int s = 0; s < nucl.GetNstates(); s++)
		std::cout 	<< std::setw(12) << std::left << s;
	std::cout	<< std::endl;	

	std::cout 	<< std::setw(12) << std::left << " "
			<< std::setw(12) << std::left << "State";
	for(int s = 0; s < nucl.GetNstates(); s++)
		std::cout	<< std::setw(12) << std::left << nucl.GetLevelJ().at(s);
	std::cout	<< std::endl;

	for(int y = 0; y < tmpMat.GetNrows(); y++){
		std::cout 	<< std::setw(12) << std::left << y
				<< std::setw(12) << std::left << nucl.GetLevelJ().at(y);
		for(int x = 0; x < tmpMat.GetNcols(); x++){
			std::cout	<< std::setw(12) << std::left << tmpMat[y][x];
		}
		std::cout	<< std::endl;
	}

}

void MiscFunctions::PrintVectorNucleus(TVectorD tmpVec, Nucleus nucl, const char *var){

	std::cout	<< std::setw(12) << std::left << "Index" 
			<< std::setw(12) << std::left << "State"
			<< std::setw(20) << std::left << var
			<< std::endl;
	for(int y = 0; y < tmpVec.GetNrows(); y++){
		std::cout 	<< std::setw(12) << std::left << y
				<< std::setw(12) << std::left << nucl.GetLevelJ().at(y)
				<< std::setw(20) << std::left << tmpVec[y]
				<< std::endl;
	}

}

double MiscFunctions::RotationFunction(double beta, int k, int kpp, int kp){ // Beta, K, X, X'

	double cosB 	= TMath::Cos(beta/2.);
	double sinB	= TMath::Sin(beta/2.);
	double ctb	= TMath::Power(cosB/sinB,2);

	int ja	= k + kp;
	int jb	= k - kp;
	int jc 	= k + kpp;
	int jd  = k - kpp;
	double b1	= TMath::Factorial(ja) * TMath::Factorial(jb) * TMath::Factorial(jc) * TMath::Factorial(jd);


	ja	= kp + kpp;
	jb	= 2*k - kp - kpp;
	double f = TMath::Power(-1,k-kp) * TMath::Power(cosB,ja) * TMath::Power(sinB,jb) * TMath::Sqrt(b1);
	
	int mis = 0;
	if(ja < 0)
		mis = -ja;
	int mas = k - kpp;
	if(kpp < kp)
		mas = k - kp;
	ja 	= kp + kpp + mis;
	jb 	= k - kpp - mis;
	jc	= k - kp - mis;
	jd	= mis;
	double b2	= TMath::Factorial(ja) * TMath::Factorial(jb) * TMath::Factorial(jc) * TMath::Factorial(jd);

	double g	= TMath::Power(-ctb,mis)/b2;
	double DJMM	= g;
	ja = mis + 1;
	if(mas >= ja){
		for(int j = ja; j<=mas; j++){
			g *= -ctb*(k-kpp-j+1)*(k-kp-j+1)/((kp+kpp+j)*j);
			DJMM += g;
		}
	}
	DJMM *= f;

	double sum = DJMM;	

	return sum;

}

unsigned int MiscFunctions::doublefactorial(unsigned int n){
	
	if( n ==0 || n == 1)
		return 1;
	return	n * doublefactorial(n-2);

}

//
//	Spherical harmonics tabulated at: https://en.wikipedia.org/wiki/Table_of_spherical_harmonics
//
double MiscFunctions::SphericalHarmonics(double theta, int l, int m, bool sym){

	if(sym){
		if(m != 0)
			return 0;
		switch (l){
			case 0:
				return 	0.50 * TMath::Sqrt(1/TMath::Pi());
			case 2:
				return	0.25 * TMath::Sqrt(5 / TMath::Pi()) * (3 * TMath::Power(TMath::Cos(theta),2) - 1);
			case 4:
				return	(3./16.) * TMath::Sqrt(1 / TMath::Pi()) * (35 * TMath::Power(TMath::Cos(theta),4) - 30 * TMath::Power(TMath::Cos(theta),2) + 3);
			case 6:
				return 	(1. / 32.) * TMath::Sqrt(13. / TMath::Pi()) * (231 * TMath::Power(TMath::Cos(theta),6) - 315 * TMath::Power(TMath::Cos(theta),4) + 105 * TMath::Power(TMath::Cos(theta),2) - 5);
			return 1;
		}
	}
	else{
		switch (l){
			case 0:
				return 	0.50 * TMath::Sqrt(1/TMath::Pi());
			case 2:
				switch (m){
					case -2:
						return	0.25 * TMath::Sqrt(15 / (2 * TMath::Pi())) * (TMath::Power(TMath::Sin(theta),2));
					case -1:
						return	0.50 * TMath::Sqrt(15 / (2 * TMath::Pi())) * (TMath::Sin(theta) * TMath::Cos(theta));
					case 0 :
						return	0.25 * TMath::Sqrt(5 / TMath::Pi()) * (3 * TMath::Power(TMath::Cos(theta),2) - 1);
					case 1 :
						return	-0.50 * TMath::Sqrt(15 / (2 * TMath::Pi())) * (TMath::Sin(theta) * TMath::Cos(theta));
					case 2 :
						return	0.25 * TMath::Sqrt(15 / (2 * TMath::Pi())) * (TMath::Power(TMath::Sin(theta),2));
					return 0;
				}
			case 4:
				switch (m){
					case -4:
						return	(3./16.) * TMath::Sqrt(35. / (2 * TMath::Pi())) * TMath::Power(TMath::Sin(theta),4);
					case -3:
						return	(3./8.) * TMath::Sqrt(35. / (TMath::Pi())) * TMath::Cos(theta) * TMath::Power(TMath::Sin(theta),3);
					case -2:
						return	(3./8.) * TMath::Sqrt(5. / (2 * TMath::Pi())) * TMath::Power(TMath::Sin(theta),2) * (7 * TMath::Power(TMath::Cos(theta),2) - 1);
					case -1:
						return	(3./8.) * TMath::Sqrt(5. / (TMath::Pi())) * TMath::Sin(theta) * (7 * TMath::Power(TMath::Cos(theta),3) - 3 * TMath::Cos(theta));
					case 0 :
						return	(3./16.) * TMath::Sqrt(1 / TMath::Pi()) * (35 * TMath::Power(TMath::Cos(theta),4) - 30 * TMath::Power(TMath::Cos(theta),2) + 3);
					case 1 :
						return	-(3./8.) * TMath::Sqrt(5. / (TMath::Pi())) * TMath::Sin(theta) * (7 * TMath::Power(TMath::Cos(theta),3) - 3 * TMath::Cos(theta));
					case 2 :
						return	(3./8.) * TMath::Sqrt(5. / (2 * TMath::Pi())) * TMath::Power(TMath::Sin(theta),2) * (7 * TMath::Power(TMath::Cos(theta),2) - 1);
					case 3 :
						return	-(3./8.) * TMath::Sqrt(35. / (TMath::Pi())) * TMath::Cos(theta) * TMath::Power(TMath::Sin(theta),3);
					case 4 :
						return	(3./16.) * TMath::Sqrt(35. / (2 * TMath::Pi())) * TMath::Power(TMath::Sin(theta),4);
					return 0;
				}
			case 6:
				switch (m){
					case -6:
						return	(1./64.) * TMath::Sqrt(3003/TMath::Pi()) * TMath::Power(TMath::Sin(theta),6);
					case -5:
						return	(3./32.) * TMath::Sqrt(1001/TMath::Pi()) * TMath::Power(TMath::Sin(theta),5) * TMath::Cos(theta);
					case -4:
						return	(3./32.) * TMath::Sqrt(91/(2*TMath::Pi())) * TMath::Power(TMath::Sin(theta),4) * (11. * TMath::Power(TMath::Cos(theta),2) - 1);
					case -3:
						return	(1./32.) * TMath::Sqrt(1365/TMath::Pi()) * TMath::Power(TMath::Sin(theta),3) * (11 * TMath::Power(TMath::Cos(theta),3) - 3 * TMath::Cos(theta));
					case -2:
						return	(1./64.) * TMath::Sqrt(1365/TMath::Pi()) * TMath::Power(TMath::Sin(theta),2) * (33 * TMath::Power(TMath::Cos(theta),4) - 18 * TMath::Power(TMath::Cos(theta),2) + 1);
					case -1:
						return	(1./16.) * TMath::Sqrt(273/(2*TMath::Pi())) * TMath::Sin(theta) * (33 * TMath::Power(TMath::Cos(theta),5) - 30 * TMath::Power(TMath::Cos(theta),3) + 5 * TMath::Cos(theta));
					case 0 :
						return 	(1. / 32.) * TMath::Sqrt(13. / TMath::Pi()) * (231 * TMath::Power(TMath::Cos(theta),6) - 315 * TMath::Power(TMath::Cos(theta),4) + 105 * TMath::Power(TMath::Cos(theta),2) - 5);
					case 1 :
						return	-(1./16.) * TMath::Sqrt(273/(2*TMath::Pi())) * TMath::Sin(theta) * (33 * TMath::Power(TMath::Cos(theta),5) - 30 * TMath::Power(TMath::Cos(theta),3) + 5 * TMath::Cos(theta));
					case 2 :
						return	(1./64.) * TMath::Sqrt(1365/TMath::Pi()) * TMath::Power(TMath::Sin(theta),2) * (33 * TMath::Power(TMath::Cos(theta),4) - 18 * TMath::Power(TMath::Cos(theta),2) + 1);
					case 3 :
						return	-(1./32.) * TMath::Sqrt(1365/TMath::Pi()) * TMath::Power(TMath::Sin(theta),3) * (11 * TMath::Power(TMath::Cos(theta),3) - 3 * TMath::Cos(theta));
					case 4 :
						return	(3./32.) * TMath::Sqrt(91/(2*TMath::Pi())) * TMath::Power(TMath::Sin(theta),4) * (11. * TMath::Power(TMath::Cos(theta),2) - 1);
					case 5 :
						return	-(3./32.) * TMath::Sqrt(1001/TMath::Pi()) * TMath::Power(TMath::Sin(theta),5) * TMath::Cos(theta);
					case 6 :
						return	(1./64.) * TMath::Sqrt(3003/TMath::Pi()) * TMath::Power(TMath::Sin(theta),6);
					return 0;
				}
			return 0;
		}

	}

	std::cout	<< l << "\t"
			<< m << std::endl;

	return 0;

}

void MiscFunctions::SplineFit(double *x, double *y, int n, double dy1, double dyn, double* ddy) {
  double y2[n];
  double u[1500];
  double qn;
  double un;

  if (dy1 > 1e30) {
    y2[0] = 0;
    u[0] = 1500;
  }
  else {
    y2[0] = -0.5;
    u[0] = (3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-dy1);
  }
  for (int i=1; i<n-1; ++i) {
    double sig = (x[i] - x[i-1])/(x[i+1]-x[i-1]);
    double p = sig * ddy[i-1] + 2;
    y2[i] = (sig-1)/p;
    u[i] = (6 * ((y[i+1]-y[i])/(x[i+1] - x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]))/(x[i+1]-x[i-1]) - sig*u[i-1])/p;    
  }
  if (dyn > 1e30) {
    qn = 0;
    un = 0;
  }
  else {
    qn = 0.5;
    un = (3/(x[n-1]-x[n-2])) * (dyn - (y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
  }
  ddy[n-1] = (un-qn*u[n-2])/(qn*ddy[n-2] + 1.);
  for (int k=n-2; k>=0; --k) {
    ddy[k] = ddy[k]*ddy[k+1] + u[k];
  }
}

void MiscFunctions::SplineEval(double *x, double *y, double *ddy, int n, double &xval, double &yval) {
  int klo = 0;
  int khi = n-1;
  int k;
  while (khi-klo > 1) {
    k = (khi + klo)/2;
    if (x[k] > xval) { khi = k; }
    else { klo = k; }
  }
  double h = x[khi] - x[klo];

  double a;
  double b;
  if (std::abs(h) < 1e-9) {
    if ( xval < x[0] ) { h = 1; }
    if ( xval > x[n-1] ) { k = n-2; }
    a = (y[k] - y[k+1])/(x[k] - x[k+1]);
    b = (y[k] + y[k+1] - a*(x[k] + x[k+1]))*0.5;
    yval = a*xval + b;
    std::cout << "SplineEval " << xval << " " << yval << " extrapolation" << std::endl;
    return;
  }
  a = (x[khi] - xval)/h;
  b = (xval - x[klo])/h;
  yval = a*y[klo] + b*y[khi] + ((std::pow(a,3) - a)*ddy[klo] + (std::pow(b,3) - b)*ddy[khi])*(h*h/6.);
                                                     
  return;
}

void MiscFunctions::Spline(double *x, double *y, int n, double xval, double &yval) {
  double *y_loc = new double[n]();
  for (int i=0; i<n; ++i) {
    y_loc[i] = std::log(y[i]);
  }
  double dy1 = (y_loc[1]-y_loc[0])/(x[1]-x[0]);
  double dyn = (y_loc[n-1]-y_loc[n-2])/(x[n-1]-x[n-2]);

  double *ddy = new double[n]();
  SplineFit(x, y_loc, n, dy1, dyn, ddy);

  SplineEval(x, y_loc, ddy, n, xval, yval);

  yval = std::exp(yval);

  delete[] y_loc;
  delete[] ddy;
  return;
}

