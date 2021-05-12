#ifndef QUADRATURE_H
#define QUADRATURE_H
#include <math.h>
#include <stdio.h>

struct Quadrature{

	int n;
	virtual double next() = 0;

};

template<class T>
struct Trapzd : Quadrature {
	double a,b,s;

	T &func;
	Trapzd() {};
	Trapzd(T &funcc, const double aa, const double bb) :
		func(funcc), a(aa), b(bb) {n=0;}

	double next() {
		double x, tnm, sum, del;
		int it,j;
		n++;
		if(n==1){
			return(s=0.5*(b-a)*(func(a) + func(b)));
		} else {
			for(it=1, j=1; j<n-1;j++) it <<= 1; //bitwise shift
			tnm=it;
			del=(b-a)/tnm;
			x=a+0.5*del;
			for(sum=0.0,j=0;j<it;j++,x+=del) sum += func(x);
			s=0.5*(s+(b-a)*sum/tnm);
			return s;
		}
	}
};

template<class T>
double qtrap(T& func, const double a, const double b, const double eps=1e-10){

	const int JMAX=20;
	double s,olds=0.0;
	Trapzd<T> t(func, a, b);
	for(int j=0; j<JMAX;j++){
		s=t.next();
		if(j>5)
			if(abs(s-olds) < eps*abs(olds) ||
				(s == 0.0 && olds == 0.0)) return s;
		olds=s;
	}
	throw("Too many steps in routine qtrap");
}

class InterpFunc {

	public:

		InterpFunc (double* xx, double* yy, unsigned long n) : 
			xarr(xx),
			yarr(yy),
			_n ( n ) {}
		double interpolate(double xx[], double yy[], unsigned long n, double x);
		void locate(double xx[], unsigned long n, double x, unsigned long *j);	
		
		double operator() (double x) { return interpolate(xarr, yarr, _n, x); }

	private:

		double* xarr;
		double* yarr;
		unsigned long _n;

};

double InterpFunc::interpolate(double xx[], double yy[], unsigned long n, double x){
  double y, dx, dy;
  unsigned long jlo;

  //locate(xx, n, x, &jlo);
  // Linear interpolation:
  
  /*
  dx = xx[jlo] - xx[jlo+1];
  dy = log10(yy[jlo]) - log10(yy[jlo+1]);
  y = dy/dx*(x-xx[jlo]) + log10(yy[jlo]);
  y=pow(10.0,y);
  */
  locate(xx,n,x,&jlo);

  dx = xx[jlo] - xx[jlo+1];
  dy = yy[jlo]-yy[jlo+1];
  y = dy/dx*(x-xx[jlo]) + yy[jlo];
  
  return(y);
}

void InterpFunc::locate(double xx[], unsigned long n, double x, unsigned long *j){
  unsigned long ju,jm,jl;
  int ascnd;

  jl=0;
  ju=n-1;
  ascnd=(xx[n-1] >= xx[0]);
  while (ju-jl > 1) {
    jm=(ju+jl)/2;
    if ((x >= xx[jm]) == ascnd)
      jl=jm;
    else
      ju=jm;
  }
  if (x == xx[0]) *j=0;
  else if(x == xx[n-1]) *j=n-2;
  else *j=jl;
  return;
}
#endif
