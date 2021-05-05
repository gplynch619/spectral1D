#include "sc.h"
#include "quadrature.h"
//#include "quadrature.h"

#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>

#include <cassert>

#define DIMENSION 3

inline int MOD(int x, int y) { return (x-y*(x/y)); }

void output(std::string fname, double* col1, std::complex<double> * col2, int N){
	std::ofstream OutFile;
	int i;
	OutFile.open(fname, std::ios::out | std::ios::app);
	for(i=0; i<N; ++i){
		OutFile << col1[i] <<" "<<pow(std::abs(col2[i]),2.0)
			<<" "<<std::arg(col2[i])<<std::endl;
	}
}

void output(std::string fname, double* col1, double* col2, int N){
	std::ofstream OutFile;
	int i;
	OutFile.open(fname, std::ios::out | std::ios::app);
	for(i=0; i<N; ++i){
		OutFile << col1[i] <<" "<<col2[i] <<std::endl;
	}
}

double vecnorm(std::vector<std::complex<double>>& a, 
		std::vector<std::complex<double>>& b, double lower, double upper, 
		std::vector<double>& grid){

	assert(a.size()==b.size());

	int N = a.size();

	std::vector<double> prod;
	prod.resize(N);

	double tmp = 0.0;
	for(int i=0; i<N; ++i){
		tmp = std::real(std::conj(a.at(i))*b.at(i));
		prod.at(i) = tmp;
	}

	InterpFunc integrand = InterpFunc(grid.data(), prod.data(), N);

	return qtrap(integrand, lower, upper);
}

//watch out this clears whatever data was in wf
void gaussian_ic(double width, double center, double phase,
		std::vector<std::complex<double>>& wf, std::vector<double>& grid){

	double tpi = 2.0*4.0*atan(1.0);
	double coeff = pow(tpi*width*width, -0.25);
	int N = grid.size();
	
	wf.clear();
	wf.resize(N);
	for(int i=0; i<N; ++i){
		double magnitude, arg, exponent;
		exponent = -1.0*pow((grid.at(i) - center)/(2.0*width), 2.0);
		magnitude = coeff*exp(exponent);
		arg = phase*grid.at(i);
		wf.at(i) = std::polar(magnitude, arg);
	}
	
}
