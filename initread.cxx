#include "initread.h"


#include <cstdio>
#include <cstring>
#include <stdexcept>
#include <iomanip>
#include <limits>

inline std::string& ltrim(std::string& s, const char* t = " \t\n\r\f\v")
{
    s.erase(0, s.find_first_not_of(t));
    return s;
}

// trim from right

inline std::string& rtrim(std::string& s, const char* t = " \t\n\r\f\v")
{
    s.erase(s.find_last_not_of(t) + 1);
    return s;
}

// trim from left & right

inline std::string trim(std::string s, const char* t = " \t\n\r\f\v")
{
    return ltrim(rtrim(s, t), t);
}

using namespace std;

Parameters::Parameters(string pName) :
	zpr_size(0),
	m_zin(-1.0),
	m_zfin(-1.0),
	m_ng(-1),
	m_rL(-1.0),
	m_nsteps(-1)
	{
		readNewParams(pName);
	}

void Parameters::readNewParams(string &fn) {
	m_params.clear();

	std::stringstream myfile;
	if (!getRank0Stream(fn.c_str(), myfile)) {
		std::ostringstream ost;
		ost << "initread cannot open '"
			<< fn
			<< "'";
		throw std::runtime_error( ost.str() );
	}

	while (myfile.good()){
		string line;
		getline(myfile, line);
		size_t pos = line.find_first_not_of(" \t");
		if(pos < line.size()){
			line = line.substr(pos);
		}
		
		if(line.empty()) {
			continue;
		} else if (line[0]== '#'){
			continue;
		}

		pos = line.find_first_of(" \t");
		string key = line.substr(0, pos);

		if(pos < line.size()){
			line = line.substr(pos);
		} else {
			line = "";
		}

		pos = line.find_first_not_of(" \t");
		if(pos < line.size()){
			line = line.substr(pos);
		} else {
			line = "";
		}

		string value = line;

		m_params.insert(std::make_pair(key, value));
	}
	
	m_hbar = atof(m_params["HBAR"].c_str());
	m_rL = atof(m_params["RL"].c_str());
	m_ng = atof(m_params["NG"].c_str());
	m_zin = atof(m_params["Z_IN"].c_str());
	m_zfin = atof(m_params["Z_FIN"].c_str());
	m_nsteps = atof(m_params["NSTEPS"].c_str());

	/*
	if(m_params["TRANS"] == "CMB"){
		m_trans = 0;
	} else if (m_params["TRANS"] == "KH"){
		m_trans = 1;
	} else if (m_params["TRANS"] == "AXCAMB"){
		m_trans = 2;
	}*/

	std::istringstream instream1(m_params["Z_PR"].c_str());
	while(!instream1.eof()){
		REAL z;
		if(instream1 >> z){
			zpr_size++;
		}
	}
	m_zpr = (REAL *)malloc(sizeof(REAL)*zpr_size);
	std::istringstream instream(m_params["Z_PR"].c_str());
	for(int i=0; i<zpr_size; i++){
		REAL z;
		instream >> z;
		m_zpr[i]=z;
	}
/*
	for(int i=0; i<zpr_size; i++){
		m_apr[i] = 1.0/(1.0 + m_zpr[i]);
	}
*/
	std::stringstream ss;
	ss.str(m_params["OUTBASE"].c_str());
	m_outBase  = trim(ss.str()); 
	ss.clear();
}

string Parameters::getParamsString(){
	std::stringstream ss;
    ss << "OUTBASE " << m_outBase << "\n";
	ss << "\n";

    ss << "NG " << m_ng << "\n";
    ss << "Z_IN " << m_zin << "\n";
    ss << "Z_FIN " << m_zfin << "\n";
    ss << "N_STEPS " << m_nsteps << "\n";
    ss << "HBAR " << m_hbar << "\n";
    ss << "RL " << m_rL << "\n";
	ss << "Z_PR ";
	for(int i=0; i<zpr_size; i++){
		ss<<m_zpr[i];
		ss<<" ";
	}
	ss<<"\n";

    return ss.str();
}
