#ifndef MGLOBALVARS_H_
#define MGLOBALVARS_H_

#include "global.h"
#include <time.h>


//
//---------------------------------------------------------------------------
// Global simulation data
//
class MGlobalVars
{
public:
	MGlobalVars();
	virtual ~MGlobalVars();

public:
	time_t sph_starttime;
	int sph_np, sph_nummat, sph_ndim, sph_nstressp, sph_nvelocp, sph_max_np;
	int sph_ssp, sph_esp, sph_svp, sph_evp;
	int sph_timestep, sph_status_interval;
	int sph_restart_interval, sph_run_restart;
	int sph_axopt, sph_disctype;
	int sph_maxnbr;
	int sph_next_restart, sph_next_run_restart;

	real sph_endtime, sph_tssf, sph_itss;
	real sph_ptime, sph_dt, sph_dtold, sph_critts, sph_init_ts;
	real lennardJonesD, lennardJonesR0;
	RealMatrix sph_coord_minmax;
};

#endif /*MGLOBALVARS_H_*/
