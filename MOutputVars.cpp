#include "MOutputVars.h"

MOutputVars::MOutputVars()
{
	sph_state_opt = 2;            // Ensight case is the default option
	sph_thnode = 0;               // default no time history particles
	sph_num_transducer = 0;       // default no pressure transducers
	sph_istate = 0;			  	  // to write new file at the beginning

	sph_plottime = 0.0;
}

MOutputVars::~MOutputVars()
{
}
