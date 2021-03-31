////////////////////////////////////
//
// initread.h
// This is the basic class that reads in the raw parameter file. Based essentially on Basedata.h 
//
/////////////////////////////////////

#ifndef INITREAD_H
#define INITREAD_H

#include <unistd.h>
#include <string>
#include <stdlib.h>
#include <map>
#include <math.h> //needed for pow

#include "sctypes.h"

#include <iostream>
#include <fstream>
#include <sstream>

class Parameters {
	
	public:
		
		Parameters(std::string pName);
		~Parameters() {}

		
		void readNewParams(std::string &fn);
		std::string getParamsString();
		
		REAL hbar() const{ return m_hbar; }
		REAL zin() const { return m_zin; }
		REAL zfin() const { return m_zfin; }
		int ng() const { return m_ng; }
		REAL rL() const { return m_rL; }
		std::string outBase() const { return m_outBase; }
		
		size_t nzpr() { return zpr_size; }
		
		REAL* zpr() { return m_zpr; }
		REAL zpr(int i) const { return m_zpr[i]; }

		int nsteps() const { return m_nsteps; }

	private:

		std::map<std::string, std::string> m_params; //map that holds parameter:value pairs from read in

		//Initializer set up
		REAL m_hbar;
		REAL m_zin;
		REAL m_zfin;
		REAL* m_zpr;
		size_t zpr_size;

		//Code Parameters
		int m_ng;
		REAL m_rL;

		int m_nsteps;
		std::string m_outBase;
};

inline bool getRank0Stream(const char *filename, std::stringstream &ss){

	ss.clear();
	std::ifstream fs(filename);
	if(!fs.is_open()){
			fprintf(stderr, "ERROR: failed to open file \"%s\"\n", filename);
			return false;
	}

	ss.str(std::string(""));
	ss << fs.rdbuf();
	fs.close();

	std::string s = ss.str();
	int sz = (int) s.size()+1;

	return true;
}
#endif
