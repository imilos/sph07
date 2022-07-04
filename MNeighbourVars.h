#ifndef MNEIGHBOURVARS_H_
#define MNEIGHBOURVARS_H_

#include "global.h"
#include <vector>

class MNeighbourVars
{
public:
	MNeighbourVars();
	virtual ~MNeighbourVars();

	void InitGrid(int ndim, IntVector &boundary_type, RealMatrix &coord_maxmin, RealMatrix &boundary_x);
	int& GetGridValue(int i, int j, int k);
	int SetGridValue(int i, int j, int k, int value);
	int GetBoxCoords(RealVector &x, IntVector &boxcrd);
	int GetLoopLimits(int ndim, IntVector &boxcrd, IntMatrix &looplim);

public:
	real sph_hmax;
	real m_factor;

	RealVector sph_gridsize;
	IntVector sph_gridlim;
	RealVector sph_gridmin;

	//IntCubeL sph_llgrid;
	int *sph_llgrid;
};

#endif /*MNEIGHBOURVARS_H_*/
