#include "MNeighbourVars.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

MNeighbourVars::MNeighbourVars() {
	sph_hmax = 0.0;
	m_factor = 2.0;
	sph_gridmin.resize(3);
	sph_gridlim.resize(3);
	sph_gridsize.resize(3);
	sph_llgrid = NULL;
}

MNeighbourVars::~MNeighbourVars() {
}

void MNeighbourVars::InitGrid(int ndim, IntVector &boundary_type, RealMatrix &coord_minmax, RealMatrix &boundary_x)
{
	// setup underlying grid
	// cell size is "m_factor" times the maximum smoothing length
	// This is done at the beginning of each time step

	sph_gridmin.clear();
	sph_gridlim(0) = sph_gridlim(1) = sph_gridlim(2) = 1;
	sph_gridsize(0) = sph_gridsize(1) = sph_gridsize(2) = 1.0
			/ (m_factor * sph_hmax);

	for (int i = 0; i < ndim; i++) {
		switch (boundary_type(i)) {
		// Both free
		case 1:
			sph_gridmin(i) = coord_minmax(0, i) - 0.01 * sph_hmax;
			sph_gridlim(i) = (int) (sph_gridsize(i)
					* (coord_minmax(1, i) - sph_gridmin(i))) + 1;
			break;
			// Min fixed, max free, to be implemented
		case 2:
			break;
			// Min free, max fixed, to be implemented
		case 3:
			break;
			// Both fixed
		case 4:
			// Adjust gridsize to ensure an integer number of cells
			sph_gridlim(i) = (int) (sph_gridsize(i)
					* (boundary_x(1, i) - boundary_x(0, i)));
			sph_gridlim(i) = MAX(sph_gridlim(i), 0);
			sph_gridsize(i) = sph_gridlim(i)
					/ (boundary_x(1, i) - boundary_x(0, i));
			sph_gridmin(i) = boundary_x(0, i) - 1.0 / sph_gridsize(i);
			sph_gridlim(i) += 2; // two cell rows for periodic boundary
			break;
		}
	}

	// Initialise linked list, sph_llgrid is an array of ints
	// Zero means empty cell
	if (sph_llgrid!=NULL) free(sph_llgrid);

	// Memory allocation for sph_llgrid
	sph_llgrid = (int *) malloc(sph_gridlim(0)*sph_gridlim(1)*sph_gridlim(2)*sizeof(int));

	for (int i=0; i<sph_gridlim(0); i++) for (int j=0; j<sph_gridlim(1); j++) for (int k=0; k<sph_gridlim(2); k++)
		SetGridValue(i,j,k,LINKEDLISTEND);
}

// Get value (first particle id) in grid on grid coordinates i,j,k - zero based

int& MNeighbourVars::GetGridValue(int i, int j, int k)
{
	return sph_llgrid[sph_gridlim(2) * (sph_gridlim(1)*i+j) + k];
}

// Set value (first particle) in grid on grid coordinates i,j,k - zero based
int MNeighbourVars::SetGridValue(int i, int j, int k, int value)
{
	sph_llgrid[sph_gridlim(2) * (sph_gridlim(1)*i+j) + k] = value;
	return OK;
}


int MNeighbourVars::GetBoxCoords(RealVector &x, IntVector &boxcrd)
{
	RealVector temp = x - sph_gridmin;
	for (int i = 0; i < 3; i++)
		boxcrd(i) = sph_gridsize(i) * temp(i);
	return OK;
}

int MNeighbourVars::GetLoopLimits(int ndim, IntVector &boxcrd, IntMatrix &looplim)
{
	for (int m = 0; m < 3; m++)
	{
		looplim(0, m) = MAX(boxcrd(m)-1, 0);
		looplim(1, m) = MIN(boxcrd(m)+1, sph_gridlim(m)-1);
	}

	for (int m = ndim; m < 3; m++)
		looplim(0, m) = looplim(1, m) = 0;

	return OK;
}

