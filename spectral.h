#ifndef SPECTRAL_H
#define SPECTRAL_H

#include <algorithm>
#include <vector>
#include <cassert>
#include <complex>
#include <fftw3.h>


//////////////////////////////////////////////////
//
//Base class for spectral solver. Specific implementations
//(i.e. specific must be defined as derived classes.
//
//
//////////////////////////////////////////////////

inline int MOD(int x, int y) { return (x-y*(x/y)); }

class SpectralSolverBase {

	public:

		//rsarray = array that will hold real space data
		//ksarray = array that will hold momentum space data
		//will probably pass a pointer to timestepper here in order to implement
		//time dependent problems

		SpectralSolverBase(int N, double rL, double hbar,
				std::vector<std::complex<double>>& rsarray, 
				std::vector<std::complex<double>>& ksarray)
			: m_Ng(N),
			m_rL(rL),
			m_planck(hbar),
			m_rspointer(rsarray.data()),
			m_kspointer(ksarray.data())

		{
			initialize();
		}

		virtual ~SpectralSolverBase() {
			//free allocations
			fftw_destroy_plan(m_fft_forward);
			fftw_destroy_plan(m_fft_backward);
		}

		//Forward
		
		void fft_forward(std::complex<double>* in, std::complex<double>* out){
			fftw_execute_dft(m_fft_forward, reinterpret_cast<fftw_complex*>(in),
					reinterpret_cast<fftw_complex*>(out));
		}
		
		void fft_backward(std::complex<double>* in, std::complex<double>* out){
			fftw_execute_dft(m_fft_backward, reinterpret_cast<fftw_complex*>(in),
					reinterpret_cast<fftw_complex*>(out));
			double scale = 1.0/m_Ng;
			for(int i=0; i<m_Ng; ++i){
				std::complex<double> norm(scale, 0.0);
				out[i] *= norm;
			}//accounts for fftw normalizaiton conventions
		}

		////////////////////////////////////////////////
		//              Momentum kick map
		//	This function compute the momentum kick in a 
		//	split operator method, i.e. the term containing
		//	exp(-i*hbar/2 * k^2 * deltaT/2). 
		//
		//	inputs:
		//		in - vector (pass by ref) containing wave function
		//			in k-space values, psi(k)
		//		dt - the duration of the kick	
		//		out - vector (pass by ref) that will contain the
		//			psi(k) after the kick. In the future this might
		//			be done in-place but for now we are doing it
		//			in a memory inefficient way.
		//
		//input: kspace array of values
		//output: rspace array after applied kick and ifft
		void map1(std::vector<std::complex<double>>& in, double dt){
			
			for(int i=0; i<m_Ng; ++i){
				double kk = m_kgrid.at(i);
				double phase = -0.5*m_planck*kk*kk*dt;
				std::complex<double> stream = std::polar(1.0, phase); 
			
				in.at(i)*=stream;
			}
		}
		
		//Initialization
		
		void initialize(){

			m_rgrid.resize(m_Ng);
			m_kgrid.resize(m_Ng);

			m_dx = m_rL/m_Ng;

			//
			m_fft_forward = fftw_plan_dft_1d(m_Ng, 
					reinterpret_cast<fftw_complex*>(m_rspointer), 
					reinterpret_cast<fftw_complex*>(m_kspointer), 
					FFTW_FORWARD, FFTW_ESTIMATE); 
	
			m_fft_backward = fftw_plan_dft_1d(m_Ng, 
					reinterpret_cast<fftw_complex*>(m_kspointer), 
					reinterpret_cast<fftw_complex*>(m_rspointer), 
					FFTW_BACKWARD, FFTW_ESTIMATE); 
		
			for(int i=0; i<m_Ng; ++i){
				m_rgrid.at(i) = i*m_dx;
			}

			int nq = m_Ng/2;

			for(int i=0; i<m_Ng; ++i){
				double kk = i;
				if(i>=nq){ kk = -MOD(m_Ng-i, m_Ng); }
				m_kgrid.at(i) = (tpi/m_rL)*kk;
				//m_kgrid.at(i) = (tpi/m_Ng)*kk; 
				//sets kspace grid in same format as fftw
			}			

		}


		double* real_grid(){ return m_rgrid.data(); }
		std::vector<double>& real_grid_vec() { return m_rgrid; }
		double* k_grid(){ return m_kgrid.data(); }

		double dx(){ return m_dx; }
		double rL(){ return m_rL; }
		int Ng(){ return m_Ng; }

		virtual void map2(std::vector<std::complex<double>>& in, double dt){};

	protected:			

		std::vector<double> m_rgrid;
		std::vector<double> m_kgrid;

		//std::vector<std::complex<double>> m_buf1;
		//std::vector<std::complex<double>> m_buf2;

		std::complex<double>* m_rspointer; 
		//internal class point which should point to the array holding rs data
		std::complex<double>* m_kspointer;
		//internal class point which should point to the array holding ks data

		fftw_plan m_fft_forward;
		fftw_plan m_fft_backward;

		const double pi = 4.0*atan(1.0);
		const double tpi = 2.0*pi;
		
		int m_Ng;
		double m_planck;
		double m_rL;
		double m_dx;
};

class SpectralFreeParticle : public SpectralSolverBase {

	public:
		SpectralFreeParticle(int N, double rL, double hbar,
				std::vector<std::complex<double>>& rsarray,
				std::vector<std::complex<double>>& ksarray) 
			: SpectralSolverBase(N,rL,hbar,rsarray,ksarray)
		{
		}

		void map2(std::vector<std::complex<double>>& in, 
				std::vector<std::complex<double>>& out, double dt){
			int i=0;
			//normally this would have the stream map but a free particle
			//has no potential so this wouldn't do anything.
			//instead of actually looping over and multiplying by
			//one I have a function that does nothin
		}

};

class SpectralQHO : public SpectralSolverBase {

	public:
		SpectralQHO(int N, double rL, double hbar,
				std::vector<std::complex<double>>& rsarray,
				std::vector<std::complex<double>>& ksarray,
				double omega) 
			: SpectralSolverBase(N,rL,hbar,rsarray,ksarray),
			m_omega(omega)
		{
			m_center = rL/2.0;
		}


		//we will start with a HO w/ time independent freq
		//and non-forced
		void map2(std::vector<std::complex<double>>& in, double dt){
			for(int i=0; i<m_Ng; ++i){
				double xx = m_rgrid.at(i);
				double phase = -0.5*m_omega*(xx - m_center)*(xx - m_center)*dt/m_planck;
				std::complex<double> kick = std::polar(1.0, phase); 
			
				in.at(i)*=kick;
			}	
		}

		double* get_potential(){
			std::vector<double> v;
			v.resize(m_Ng);
			for(int i=0; i<m_Ng; ++i){
				double xx = m_rgrid.at(i);
				v.at(i) = -0.5*m_omega*(xx-m_center)*(xx - m_center)/m_planck;
			}
			return v.data();
		}

	private:

		double m_omega;
		double m_center;
};
#endif
