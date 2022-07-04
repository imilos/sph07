#ifndef MSYMPERVARS_H_
#define MSYMPERVARS_H_

#include "global.h"

// Symmetry + periodic
class MSymPerVars
{
public:
	MSymPerVars();
	virtual ~MSymPerVars();

public:
	// New master arrays
	//
	IntVector sph_boundary_type;
	IntMatrix sph_boundary_code;
	RealMatrix sph_boundary_x;
	//
	int sph_g_maxnbr;
	IntMatrixL sph_g_nbrlist;
	//
	// New multiply and add arrays
	//
	RealMatrix sph_mincoord_mult, sph_maxcoord_mult;
	RealMatrix sph_mincoord_add, sph_maxcoord_add;
};

#endif /*MSYMPERVARS_H_*/
