#include "TimeStepper.h"

TimeStepper::TimeStepper(REAL zin_,
		REAL zfin_,
		int nsteps_,
		REAL* zpr_,	
		size_t nzpr_) :
	m_zz(-1.0),
	m_tau(-1.0),
	m_tau2(-1.0),
	m_zin(-1.0),
	m_zfin(-1.0),
	m_nsteps(-1),
	m_currentHalfStep(0)
{
	m_nsteps = nsteps_;
	m_zin = zin_; 
	m_zfin = zfin_; 
	m_tau = (m_zfin - m_zin)/(1.0*m_nsteps);
	m_tau2 = 0.5*m_tau;
	m_zz = m_zin;

	m_zpr.assign(zpr_, zpr_ + nzpr_);
}

TimeStepper::~TimeStepper() {}

void TimeStepper::set_step(REAL da){
	if(da==0){
		m_tau = (m_zfin - m_zin)/(1.0*m_nsteps);
		m_tau2 = 0.5*m_tau;
	} else {
		m_tau = da;
		m_tau2 = 0.5*m_tau;
	}
}

void TimeStepper::advanceHalfStep() {
	m_zz += m_tau2; 
	m_currentHalfStep++;
	return;
}

/*
void TimeStepper::reverseHalfStep() {
	m_pp -= m_tau2;
	m_aa = pow(m_pp, t.0/m_alpha);
	m_zz = 1.0/m_aa - 1.0;
	set_adot();
	set_scal();
	m_currentHalfStep--;
	return;
}
*/
void TimeStepper::advanceFullStep() {
	advanceHalfStep();
	advanceHalfStep();
	return;
}
/*
void TimeStepper::reverseFullStep() {
	reverseHalfStep();
	reverseHalfStep();
	return;
}*/
