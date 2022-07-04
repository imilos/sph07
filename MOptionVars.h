#ifndef MOPTIONVARS_H_
#define MOPTIONVARS_H_

#include "global.h"

class MOptionVars
{
public:
	MOptionVars();
	virtual ~MOptionVars();

public:
	bool sph_repulsive_force;
	bool sph_boundary;
	bool sph_drelax;
	bool sph_lennardjones;
	bool sph_surface_tension;
	bool sph_rigid_body;
	int sph_contacttype, sph_krtype, sph_init_h_opt, sph_init_v_opt, sph_init_rhoe;
	int sph_massopt, sph_tcrit_opt, sph_veloc_opt, sph_h_opt;
	int sph_nlcur, sph_nthpx, sph_nthpy, sph_nthpz;
	int sph_ncontmats, sph_cont_opt;
	real sph_drelax_scale;

	// For cylinder example, angle velocities
	bool bConstrainVelocityCylinder;
	double InnerW, OuterW;
};

#endif /*MOPTIONVARS_H_*/
