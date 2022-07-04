#ifndef MFILEHANDLINGVARS_H_
#define MFILEHANDLINGVARS_H_

#include "global.h"
#include <vector>
#include <string>

class MFileHandlingVars
{
public:
	MFileHandlingVars();
	virtual ~MFileHandlingVars();
	
	bool FileExists(const std::string &filename);
	int OpenInputFile();
	int OpenLogFile();
	int OpenCaseFiles(int status);
	int OpenOutputFile(const std::string &ext, FILE **file, int status);
	int CloseLastOpenedFile();
	int CloseAllButInputAndLog();
	
public:
	std::string sph_filein, sph_fileout;
	std::string sph_filerestart;
	std::string sph_title;
	
	int sph_filelen[3];
	// The following vector contains pointers to all open files
	std::vector<FILE **> sph_openfiles;
	
	// shortcuts to file variables
	FILE *f_inputfile;
	FILE *f_logfile;
	FILE *f_case;
	FILE *f_geo;
	FILE *f_sig;
	FILE *f_dis;
	FILE *f_vel;
	FILE *f_acc;
	FILE *f_eps;
	FILE *f_rho;
	FILE *f_egy;
	FILE *f_snd;
	FILE *f_mas;
	FILE *f_bnd;
	FILE *f_nbr;
	FILE *f_tmp;
	FILE *f_time;
};

#endif /*MFILEHANDLINGVARS_H_*/
