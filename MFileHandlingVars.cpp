#include "MFileHandlingVars.h"
#include <stdio.h>

MFileHandlingVars::MFileHandlingVars()
{
	f_inputfile = NULL;
	f_logfile = NULL;
	f_case = NULL;
	f_geo = NULL;
	f_sig = NULL;
	f_dis = NULL;
	f_vel = NULL;
	f_acc = NULL;
	f_eps = NULL;
	f_rho = NULL;
	f_egy = NULL;
	f_snd = NULL;
	f_mas = NULL;
	f_bnd = NULL;
	f_nbr = NULL;
	f_tmp = NULL;
	f_time = NULL;
}

MFileHandlingVars::~MFileHandlingVars()
{
	// Close all open files
	for (int i=sph_openfiles.size()-1; i>0; i--)
	{
		if (*sph_openfiles[i]!=NULL) fclose(*sph_openfiles[i]);
		*sph_openfiles[i] = NULL;
		sph_openfiles.pop_back();
	}
}

bool MFileHandlingVars::FileExists(const std::string &filename)
{
    FILE *file;

    if ( (file = fopen(filename.c_str(), "r")) )
    {
        fclose(file);
        return true;
    }
    return false;
}

int MFileHandlingVars::OpenInputFile()
{
	// sph_openfiles[0] == input file
	std::string filename = sph_filein + ".mcm";

 	// check that input file exists, if it does not then error termination
 	if ( !FileExists(filename) )
 	{
		printf("\nError - input file %s does not exist.\n", filename.c_str());
  		return ERROR;
 	}
 	// Open input file
	if ( (f_inputfile=fopen(filename.c_str(), "rt"))==NULL )
	{
		printf("\tError: Couldn't open input file.\n");
		return ERROR;
	}
	sph_openfiles.push_back(&f_inputfile);

 	return OK;
}

int MFileHandlingVars::OpenLogFile()
{
	// sph_openfiles[1] == log file

	std::string filename = sph_filein + ".log";

	if ( (f_logfile=fopen(filename.c_str(), "wt"))==NULL )
	{
		printf("\tError: Couldn't open log file.\n");
		return ERROR;
	}
	sph_openfiles.push_back(&f_logfile);

 	return OK;
}

// Open all Ensight CASE files
int MFileHandlingVars::OpenCaseFiles(int status)
{
	if (OpenOutputFile("_sig.ens",&f_sig, status)!=OK) return ERROR;
	if (OpenOutputFile("_dis.ens",&f_dis, status)!=OK) return ERROR;
	if (OpenOutputFile("_vel.ens",&f_vel, status)!=OK) return ERROR;
	if (OpenOutputFile("_acc.ens",&f_acc, status)!=OK) return ERROR;
	if (OpenOutputFile("_eps.ens",&f_eps, status)!=OK) return ERROR;
	if (OpenOutputFile("_rho.ens",&f_rho, status)!=OK) return ERROR;
	if (OpenOutputFile("_egy.ens",&f_egy, status)!=OK) return ERROR;
	if (OpenOutputFile("_snd.ens",&f_snd, status)!=OK) return ERROR;
	if (OpenOutputFile("_mas.ens",&f_mas, status)!=OK) return ERROR;
	if (OpenOutputFile("_bnd.ens",&f_bnd, status)!=OK) return ERROR;
	if (OpenOutputFile("_nbr.ens",&f_nbr, status)!=OK) return ERROR;
	if (OpenOutputFile("_tmp.ens",&f_tmp, status)!=OK) return ERROR;

	return OK;
}

// Status should be true for first time step, false otherwise
int MFileHandlingVars::OpenOutputFile(const std::string &ext, FILE **file, int status)
{
	std::string filename;
	std::string symbol;

	filename = sph_filein + ext;

	// To open new or to append to a file, new file at zero time step
	if (status==0)
		symbol="wt";
	else
		symbol="at";

	if ( (*file=fopen(filename.c_str(), symbol.c_str()))==NULL )
	{
		printf("\tERROR: Couldn't open Ensight CASE file for writting.\n");
		return ERROR;
	}

	sph_openfiles.push_back(file);

	return OK;
}

int MFileHandlingVars::CloseLastOpenedFile()
{
	int ilast=sph_openfiles.size()-1;

	if (*sph_openfiles[ilast]!=NULL) fclose(*sph_openfiles[ilast]);
	*sph_openfiles[ilast] = NULL;
	sph_openfiles.pop_back();

	return OK;
}

// Routine all files except input and log files 0 and 1
int MFileHandlingVars::CloseAllButInputAndLog()
{
	for (int i=sph_openfiles.size()-1; i>1; i--)
	{
		if (*sph_openfiles[i]!=NULL) fclose(*sph_openfiles[i]);
		*sph_openfiles[i] = NULL;
		sph_openfiles.pop_back();
	}

	return OK;
}
