#ifndef TIMESTEPPER_H
#define TIMESTEPPER_H

#include <cmath>
#include <vector>
#include <iostream>

typedef double REAL;

class TimeStepper {
	public:

		TimeStepper(REAL zin_,
				REAL zfin_,
				int nsteps_,
				REAL* zpr_,
				size_t nzpr_);
		~TimeStepper();

		void advanceHalfStep();
		void advanceFullStep();

		REAL zz() { return m_zz; }
		REAL tau() { return m_tau; }
		REAL tau2() { return m_tau2; }
		REAL zin() { return m_zin; }
		REAL zfin() {return m_zfin; }
		int nsteps() { return m_nsteps; }
		int currentHalfStep() {return m_currentHalfStep; }
		
		std::vector<REAL> m_zpr;
		std::vector<REAL> m_apr;

		void set_step(REAL da=0);

	private:

		TimeStepper();
		TimeStepper( const TimeStepper& );
		TimeStepper& operator = (const TimeStepper& );

		REAL m_zz;
		REAL m_tau;
		REAL m_tau2;
		REAL m_zin;
		REAL m_zfin;
		int m_nsteps;
		int m_currentHalfStep;
};

#endif
