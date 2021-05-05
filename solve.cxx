#include "spectral.h"
#include "initread.h"
#include "TimeStepper.h"
#include "sctypes.h"
#include "sc.h"

#include <string>
#include <fstream>
#include <sstream>
#include <vector>

#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>

#include <math.h>

#define EPS 1e-8

//if(k_i>nq){k_i = -MOD(ng-k_i, ng);}

void print_grid(double* grid, int N){
	std::ostringstream ss;

	for(int i=0; i<N; ++i){
		ss<<grid[i]<<" ";
	}

	std::cout<<ss.str()<<std::endl;

}

void assign_gaussian(std::vector<std::complex<double>>& arr, 
		SpectralSolverBase* solver){
	
	int ng = solver->Ng();
	const double pi = 4.0*atan(1.0);
	const double tpi = 2.0*pi;

	double physical_center = solver->dx()*ng/2.;

	double width=1.0;
	double coeff = pow(pi, -0.25)*pow(1.0, -0.5);

	double* xgrid = solver->real_grid();

	for(int i=0; i<ng; ++i){
		double physx = xgrid[i];
		double shift = ((physx<=physical_center) ? physx : solver->rL()-physx);	
	
		shift=shift-physical_center;
		double mag = coeff*exp(-(shift*shift)/2.);

		arr.at(i) = std::polar(mag, 0.0);
	
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

	print_grid(solver->k_grid(), Ng);
	
	std::cout<<"Free particle solver successfully initialized"<<std::endl;

	TimeStepper ts(params.zin(), 
			params.zfin(), 
			params.nsteps(), 
			params.zpr(), 
			params.nzpr());

	std::cout<<"Time stepper successfully initialized"<<std::endl;

	//Initialize
	//assign_gaussian(rspace, solver);

	gaussian_ic(1.0, L/2.0, 0.0, rspace, solver->real_grid_vec());
	
	REAL previous_step=ts.tau();
	bool just_printed = false;


	//Actual stepping routine
	for(int step=0; step<params.nsteps(); step++){
	
		std::cout<<"TIMESTEP: "<<step<<std::endl;
		std::cout<<"z: "<<ts.zz()<<std::endl;

		previous_step = ts.tau();
		ts.set_step(); //sets tau to the default (afin-ain)/nstep

		//note, some of this code is convoluting, but those portions are mostly
		//just to get the printing right

		if(ts.m_zpr.front() - ts.zz() < ts.tau()){
			ts.set_step(ts.m_zpr.front() - ts.zz());
			std::cout<<"Adjusting step"<<std::endl;
		}

		if(just_printed){
			ts.set_step();
			ts.set_step(2.0*ts.tau() - previous_step);
			just_printed=false;
		}

		std::cout<<"tau: "<<ts.tau()<<std::endl;

		//////////////////////////////
		//
		// STREAM STEP
		//
		//////////////////////////////

		solver->fft_forward(rspace.data(), kspace.data());

		solver->map1(kspace, ts.tau2());

		solver->fft_backward(kspace.data(), rspace.data());

		ts.advanceHalfStep();
		//////////////////////////////
		//
		// KICK STEP GOES HERE
		//
		//////////////////////////////
		
		//////////////////////////////
		//
		// STREAM STEP
		//
		//////////////////////////////
	
		solver->fft_forward(rspace.data(), kspace.data());

		solver->map1(kspace, ts.tau2());

		solver->fft_backward(kspace.data(), rspace.data());

		ts.advanceHalfStep();
	
		if((ts.m_zpr.front()-ts.zz())<EPS || step+1==params.nsteps()){
			std::cout<<"Printing on step "<<step<<std::endl;
			std::ostringstream outname;
			outname<<params.outBase()<<".z"<<ts.zz();
			output(outname.str(), solver->real_grid(), rspace.data(), Ng);	
		
			ts.m_zpr.erase(ts.m_zpr.begin()); //updates zpr list
			just_printed=true;
		}

		std::cout<<"Step "<<step<<" done!"<<std::endl;
	
	}

	return 0;
}
