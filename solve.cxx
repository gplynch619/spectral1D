#include "spectral.h"
#include "initread.h"

#include <string>
#include <fstream>
#include <sstream>
#include <vector>

#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>

#include <math.h>

//if(k_i>nq){k_i = -MOD(ng-k_i, ng);}

void output(std::string fname, double* col1, std::complex<double> * col2, int N){
	std::ofstream OutFile;
	int i;
	OutFile.open(fname, std::ios::out | std::ios::app);
	for(i=0; i<N; ++i){
		OutFile << col1[i] <<" "<<col2[i].real()<<" "<<col2[i].imag()<<std::endl;
	}
}

void assign_delta(std::complex<double>* arr, int N){
	for(int i=0; i<N; ++i){
		if(i==0){ 
			arr[i].real(1.0); 
			arr[i].imag(0.0); 
		}
		else {
			arr[i].real(0.0);
			arr[i].imag(0.0);
		}
	}
}

int main(int argc, char *argv[]){

//parameter read in
//
	int argvi = optind;
	std::string paramfile;

	paramfile = argv[argvi];

	Parameters params(paramfile);

	std::string paramout = params.getParamsString();

	std::cout << paramout << std::endl;

//allocation

	int Ng = params.ng();
	double L = params.rL();
	std::vector<std::complex<double>> rspace;
	std::vector<std::complex<double>> kspace;

	rspace.resize(Ng);
	kspace.resize(Ng);

	SpectralFreeParticle* solver = new SpectralFreeParticle(Ng, L, 
			params.hbar(), rspace, kspace);

	std::cout<<"Free particle solver successfully initialized"<<std::endl;

	assign_delta(rspace.data(), Ng);

	std::ostringstream outname;
	outname<<"initial.txt";

	output(outname.str(), solver->real_grid(), rspace.data(), Ng);

	solver->fft_forward(rspace.data(), kspace.data());

	outname.str("");
	outname<<"forward.txt";

	output(outname.str(), solver->k_grid(), kspace.data(), Ng);

	solver->fft_backward(kspace.data(), rspace.data());

	outname.str("");
	outname<<"backward.txt";

	output(outname.str(), solver->real_grid(), rspace.data(), Ng);
	
	return 0;
}

