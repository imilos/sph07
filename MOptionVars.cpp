#include "MOptionVars.h"

MOptionVars::MOptionVars()
{
	sph_init_v_opt=0;
	sph_massopt = 0;                // default mass given for material
	sph_init_h_opt = 0;             // default smoothing length given for material
	sph_init_rhoe = 0;              // default initial density and energy from material cards
	sph_contacttype = 0;            // default contact type = 0 kernel contact
	sph_tcrit_opt = 0;              // default crtitical timestep calculation
	sph_veloc_opt = 0;              // default particles move with own velocity
	sph_boundary = false;           // no boundary planes
	sph_nlcur = 0;                  // number of  load curves
	sph_h_opt = 0;                  // default constant smoothing length
	sph_krtype = 1;                 // default kernel type =1 b-spline
	sph_nthpx = 0;                  // no base acceleration in X-direction
	sph_nthpy = 0;                  // no base acceleration in Y-direction
	sph_nthpz = 0;                  // no base acceleration in Z-direction
	sph_drelax = false;             // Default standard analysis
	sph_lennardjones = false;		// No Lennard-Jones particles
	sph_surface_tension = false;	// Surface tension disabled by default
	sph_rigid_body = false;			// Rigid body disabled by default

	bConstrainVelocityCylinder = false; // temporary variable
}

MOptionVars::~MOptionVars()
{
}
