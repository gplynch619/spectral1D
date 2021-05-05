#ifndef SC_H 
#define SC_H 

#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <complex>

#include "sctypes.h"

#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>

void output(std::string fname, double* col1, std::complex<double> * col2, int N);
void output(std::string fname, double* col1, double* col2, int N);

double vecnorm(std::vector<std::complex<double>>& a, std::vector<std::complex<double>>& b, double lower, double upper, std::vector<double>& grid);

void gaussian_ic(double width, double center, double phase,
		std::vector<std::complex<double>>& wf, std::vector<double>& grid);
#endif
