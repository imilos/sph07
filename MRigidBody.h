#ifndef MRIGIDBODY_H_
#define MRIGIDBODY_H_

#include "global.h"
#include "MParticle.h"
#include "MParticleData.h"

class MRigidBody
{
public:
	MRigidBody(int ID, real totalMass, real inertiaMoment, real particleMass, real deltap);
	virtual ~MRigidBody();

	real B(real x, real y, real c, real h);
	int CalculateCenterOfMass(MParticleData& par);

// Attributes
public:
	// Rigid body ID is a negative material number in input file starting from -2
	// since -1 is reserved for Lennard-Jones particles
	int ID;
	real totalMass;
	real inertiaMoment;
	real particleMass;
	real deltap; 		// boundary particle spacing
	RealVector V;		// rigid body velocity
	RealVector OMEGA;	// rigid body angle velocity
	RealVector R;		// rigid body center of mass
	RealVector A;		// rigid body linear acceleration
	RealVector ALPHA;	// rigid body angle acceleration
};

#endif /* MRIGIDBODY_H_ */
