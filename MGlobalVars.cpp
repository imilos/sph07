#include "MGlobalVars.h"

MGlobalVars::MGlobalVars()
{
	sph_starttime = time((time_t *)NULL);

	sph_disctype = 0;			   // default discretization type ==0: basic SPH
	sph_max_np = 0;                // initialise max number of particles
	sph_itss = 0.0;                // default initial time step size
	sph_status_interval = 100;     // default 100 steps between status reports to screen
	sph_restart_interval = 0;      // default no restart fie written
	sph_run_restart = 0;           // default no running restart file
	sph_maxnbr = 40;               // default maximum number of neighbours
	sph_ptime = 0.0;			   // start time of the simulation
	sph_timestep = 0;

	//sph_coord_minmax.resize(2,3);

	for (int i=0; i<3; i++)
	{
		sph_coord_minmax(0,i) =  1.0e20;
		sph_coord_minmax(1,i) = -1.0e20;
	}
}

MGlobalVars::~MGlobalVars()
{
}
