#ifndef MOUTPUTVARS_H_
#define MOUTPUTVARS_H_

#include "global.h"
#include <string>
#include <vector>

class MOutputVars
{
public:
	MOutputVars();
	virtual ~MOutputVars();
	
public:
	int sph_state_opt;
	int sph_thnode, sph_num_transducer;
	real sph_stpltime, sph_thpltime;
	
	int sph_istate;
	real sph_nextsttime, sph_nextthtime;
	
	int *sph_thnodeid;
	
	real sph_plottime;
	RealVector sph_plot_start_x, sph_plot_end_x;
	
	// Ensight case file variables
	std::string sph_filename_sig, sph_filename_rho, sph_filename_snd, sph_filename_egy;
	std::string sph_filename_mas, sph_filename_bnd, sph_filename_nbr, sph_filename_tmp;
	std::vector<real> sph_timenum;
};

#endif /*MOUTPUTVARS_H_*/
