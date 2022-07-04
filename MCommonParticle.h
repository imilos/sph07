#ifndef MCOMMONPARTICLE_H_
#define MCOMMONPARTICLE_H_

#include "global.h"
#include <vector>

class MCommonParticle
{
public:
	MCommonParticle();
	virtual ~MCommonParticle();

	void CalculateVabs(int ndim);
	real Distance(const MCommonParticle &other, int ndim);
	real DistanceSquared(const MCommonParticle &other, int ndim);
	int AddNeighbour(int id);
	int RemoveNeighbours();
	int ReserveMemForNbrs(int size);
	inline bool ToCompute()
		{ return active && !fromOtherProc; }

public:
	bool active;                     // Particle deletion flag
	bool fromOtherProc;				 // If particle is taken from another MPI process
	bool lennardJones;				 // If particle is of Lennard-Jones type (Monaghan)
	int mat;                         // Material model id
  	int nnbr;						 // Number of neighbours
  	int llpointer;                   // Linked list pointer
  	int color;						 // Color denotes material phase used for edge detection

  	real N;							 // For normals correction
  	real mass;                       // Mass
  	real h;                          // Smoothing length
  	real hold;                       // Smoothing length for previous timestep
  	real rho;                        // Density
  	real vabs;						 // Absolute velocity
  	real p;							 // Pressure
  	real c;                          // speed of sound
  	real normalizedNormalDiv;		 // Normalized surface normal divergence

  	// vectors
  	RealVector x;            		 // Current particle coordinates
  	RealVector xzero;        		 // 'Initial' particle coordinates
  	RealVector v;            		 // Current particle velocity
  	RealVector normal;				 // Surface normal, computed according to Morris 2000.
  	RealVector normalizedNormal;	 // Normalized surface normal

  	// tensors
  	RealMatrix sigma;      		     // stress tensor
  	RealMatrix q;          		     // artificial viscosity stress tensor

  	// Neighbourhood list
  	IntVectorL m_nbrlist;

};

#endif /*MCOMMONPARTICLE_H_*/
