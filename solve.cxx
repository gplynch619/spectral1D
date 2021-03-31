#include <fftw3.h>
#include <initread.h>

#include <string>
#include <fstream>
#include <sstream>

#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>

#include <math.h>

inline int MOD(int x, int y) { return (x - y*(x/y)); }
//if(k_i>nq){k_i = -MOD(ng-k_i, ng);}

void output(std::string fname, double* col1, fftw_complex * col2, int N){
	std::ofstream OutFile;
	int i;
	OutFile.open(fname, std::ios::out | std::ios::app);
	for(i=0; i<N; ++i){
		OutFile << col1[i] <<" "<<col2[i][0]<<" "<<col2[i][1]<<std::endl;
	}
}

void assign_delta(fftw_complex* arr, int N){
	fftw_complex zero,one;
	
	for(int i=0; i<N; ++i){
		if(i==0){ 
			arr[i][0] = 1.0; 
			arr[i][1] = 0.0; 
		}
		else {
			arr[i][0] = 0.0;
			arr[i][1] = 0.0;
		}
	}
}

int main(int argc, char *argv[]){

//parameter read in

//allocation 
	int N = 64;
	int nq=N/2;
	const double pi = 4.0*atan(1.0);
	const double tpi = 2.0*pi;
	
	double t[N]; //time grid
	fftw_complex signal[N]; //signal, sampeled at discrete t
	double freq[N];	
	fftw_complex out[N]; //array to hold output
	fftw_plan planf;
	fftw_plan planb;

	planf = fftw_plan_dft_1d(N, signal, out,  FFTW_FORWARD, FFTW_ESTIMATE); 
	planb = fftw_plan_dft_1d(N, out, signal, FFTW_BACKWARD, FFTW_ESTIMATE); 

	double T = tpi;
	double DeltaT = T/N;
	double omega = tpi*1.0;


	assign_delta(signal, N);

	for(int i=0; i<N; ++i){
		t[i] = i*DeltaT; //t holds evenly spaced values from 0 to T;
		freq[i]=i/T;
		if(i>=nq){ freq[i] = -MOD(N-i, N)/T; }
		//puts frequency array in correct format. The i=0 term is the 0 frequency
		//then 1<=i<N/2 - 1 are positive frequencies, and N/2+1<=i<N-1 are neg
		//with i=N/2 being f_c
	}

	std::ostringstream outname;
	outname<<"before.txt";

	output(outname.str(), t, signal, N);

	check_space(signal, N, 0);
	
	fftw_execute(planf);

	outname.str("");
	outname<<"forward.txt";

	output(outname.str(), freq, out, N);

	check_space(out, N, 1);

	fftw_execute(planb);

	outname.str("");
	outname<<"backward.txt";

	output(outname.str(), t, signal, N);

	check_space(signal, N, 0);
	
	fftw_destroy_plan(planf);
	fftw_destroy_plan(planb);

	return 0;
}

