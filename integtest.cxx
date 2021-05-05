#include "quadrature.h"

#include <string>
#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>

#include <math.h>

int main(int argc, char *argv[]){
	int N = 1024;

	double xg[N];
	double yg[N];

	double L = 5.0;
	double dx = L/N;

	for(int i=0; i<N; ++i){
		xg[i] = i*dx;
		yg[i] = i*dx*i*dx;
	}

	InterpFunc* myFunc = new InterpFunc(xg, yg, N);

	double result;

	result = qtrap( *myFunc, 0.0, L);

	std::cout<<"Result: "<<result<<std::endl;
	

	return 0;
}

