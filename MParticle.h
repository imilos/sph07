#ifndef MPARTICLE_H_
#define MPARTICLE_H_

#include "global.h"
#include "MCommonParticle.h"
#include <vector>

class MParticle : public MCommonParticle
{
public:
	MParticle();
	virtual ~MParticle();
	void CalculateAabs(int ndim);
	void CalcTraceROD(int ndim);

public:
	int n;
  	int dispbc;                     // Displacement boundary condition
  	int delpointer;                 // inactive particle linked list pointer

  	int g_nnbr;                     // number of ghost neighbours
  	int boundary;                   // flag for boundary particles
  	int nsym;                       // flag for symmetry planes

  	real h0;                         // Initial smoothing length (fixed during initialisation)
  	real rho0;                       // Initial density (fixed during initialisation)
  	real rhoold;
  	real e;                          // TOTAL particle internal energy
  	real etry;                       // trial particle internal energy
  	real einc;                       // particle internal energy increment
  	real pcut;                       // Cutoff pressure
  	real p_cut;                      // Cutoff pressure at different time
  	real efps;                       // Effective plastic strain
  	real fail;                       // Failure, varies between undamaged (1.0) and failed (0.0)
  	real temper;                     // Temperature, unused at this time
  	real mindist;                    // Distance to nearest neighbour
  	real critts;                     // Particle critical timestep
  	real aabs;						 // Absolute value of acceleration
  	real tracerod;					 // Trace of ROD
  	RealVector epx;         		 // material history variables
  	RealMatrix alfa;				 // back stress (kinematic hardening)
  	// vectors
  	RealVector smooth_v;     		 // Interpolated particle velocity, currently used only for XSPH
  	RealVector a;            		 // Particle acceleration
  	RealVector bndnorm;      		 // Unit normal vector to surface for boundary particles
  	RealVector repulsion;			 // contact acceleration
  	// tensors
  	RealMatrix rod;        		     // Rate-of-deformation tensor
  	RealMatrix spin;       		     // Spin tensor
  	RealMatrix s;	        		 // deviatoric stress tensor
  	RealMatrix qold;       		     // artificial viscosity stress tensor from previous timestep
  	RealMatrix qq;         		     // global to current configuration rotation matrix
  	RealMatrix qr;         		     // material axes rotation matrix for orthotropic materials
};

#endif /*MPARTICLE_H_*/
