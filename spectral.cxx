#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>

#include <math.h>

inline int MOD(int x, int y) { return (x - y*(x/y)); }

void initialize_rspace_grid(double* rgrid, int N, double L){


}	

void output(std::string fname, double* col1, fftw_complex * col2, int N){
	std::ofstream OutFile;
	int i;
	OutFile.open(fname, std::ios::out | std::ios::app);
	for(i=0; i<N; ++i){
		OutFile << col1[i] <<" "<<col2[i][0]<<" "<<col2[i][1]<<std::endl;
	}
}
