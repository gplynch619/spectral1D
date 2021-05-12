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
	std::vector<std::complex<double>> initial_copy;
	std::vector<double> corrfunc;
	std::vector<double> tgrid;

	rspace.resize(Ng);
	kspace.resize(Ng);
	initial_copy.resize(Ng);

	double omega=1.0;

	//SpectralFreeParticle* solver = new SpectralFreeParticle(Ng, L, 
	//		params.hbar(), rspace, kspace);
	
	SpectralQHO* solver = new SpectralQHO(Ng, L, 
			params.hbar(), rspace, kspace, omega);

	
	std::ostringstream potname;
	potname<<params.outBase()<<".potential";
	output(potname.str(), solver->real_grid(), solver->get_potential(), Ng);	

	std::cout<<"Free particle solver successfully initialized"<<std::endl;

	TimeStepper ts(params.zin(), 
			params.zfin(), 
			params.nsteps(), 
			params.zpr(), 
			params.nzpr());

	corrfunc.resize(params.nsteps());
	tgrid.resize(params.nsteps());

	std::cout<<"Time stepper successfully initialized"<<std::endl;

	//Initialize
	//assign_gaussian(rspace, solver);

	gaussian_ic(1.0, L/2.0, 0.0, rspace, solver->real_grid_vec());
	
	initial_copy.assign(rspace.begin(), rspace.end());

	REAL previous_step=ts.tau();
	bool just_printed = false;

	//check norm
	double norm = vecnorm(rspace, rspace, 0, L, solver->real_grid_vec());

	std::cout<<"Norm of wf: "<<norm<<std::endl;
	
	std::ostringstream ost;
	ost<<params.outBase()<<".ini";
	output(ost.str(), solver->real_grid(), rspace.data(), Ng);	
	
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
		
		solver->map2(rspace, ts.tau());		

		ts.advanceHalfStep();
		//////////////////////////////
		//
		// STREAM STEP
		//
		//////////////////////////////
	
		solver->fft_forward(rspace.data(), kspace.data());

		solver->map1(kspace, ts.tau2());

		solver->fft_backward(kspace.data(), rspace.data());

		if((ts.m_zpr.front()-ts.zz())<EPS || step+1==params.nsteps()){
			std::cout<<"Printing on step "<<step<<std::endl;
			std::ostringstream outname;
			outname<<params.outBase()<<".z"<<ts.zz();
			output(outname.str(), solver->real_grid(), rspace.data(), Ng);	
		
			ts.m_zpr.erase(ts.m_zpr.begin()); //updates zpr list
			just_printed=true;
		}

		double corr; 
		
		corr = vecnorm(initial_copy, rspace, 0, L, solver->real_grid_vec());

		corrfunc.at(step)=corr;
		tgrid.at(step) = ts.currentHalfStep()*ts.tau2();	

		std::cout<<"Step "<<step<<" done!"<<std::endl;
	
	}

	std::cout<<"test"<<std::endl;
	write_energy_spectrum(corrfunc, tgrid, solver, params.nsteps(), params.zfin(), params.outBase());

	return 0;
}
