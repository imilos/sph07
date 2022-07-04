#include "MSymPerVars.h"

MSymPerVars::MSymPerVars()
{
	sph_boundary_type.resize(3);
	sph_boundary_type(0)=sph_boundary_type(1)=sph_boundary_type(2)=1;
	sph_boundary_code.resize(3,3); sph_boundary_code.clear();
	sph_boundary_x.resize(3,3); sph_boundary_x.clear();

	sph_mincoord_mult.resize(3,3); sph_mincoord_mult.clear();
	sph_maxcoord_mult.resize(3,3); sph_maxcoord_mult.clear();
	sph_mincoord_add.resize(3,3); sph_mincoord_add.clear();
	sph_maxcoord_add.resize(3,3); sph_maxcoord_add.clear();
}

MSymPerVars::~MSymPerVars()
{
}
