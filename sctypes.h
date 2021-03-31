#ifndef SC_TYPES_H
#define SC_TYPES_H

#include <stdint.h>

typedef double REAL; 

#ifdef GRID_64
typedef double GRID_T;
#else
typedef long GRID_T;
#endif

typedef struct{
	REAL mag;
	REAL phase;
} complex_p;

#endif

