#include "MRigidBody.h"

MRigidBody::MRigidBody(int ID, real totalMass, real inertiaMoment, real particleMass, real deltap)
{
	this->ID = ID;
	this->totalMass = totalMass;
	this->inertiaMoment = inertiaMoment;
	this->particleMass = particleMass;
	this->deltap = deltap;

	V.resize(3); V.clear();
	R.resize(3); R.clear();
	OMEGA.resize(3); OMEGA.clear();
	A.resize(3); A.clear();
	ALPHA.resize(3); ALPHA.clear();
}

MRigidBody::~MRigidBody()
{
}

// According to Monaghan 2003, sound speed taken into accounting of factor "beta"
real MRigidBody::B(real x, real y, real c, real h)
{
	real hi = x<deltap ? 1-x/deltap : 0;
	real q = y/h;
	real beta = 0.02*SQR(c)/(y+0.01*h);

	if (q<2./3) return hi * 2/3*beta;
	else if (q<1) return hi * beta*(2*q-1.5*SQR(q));
	else if (q<2) return hi * 0.5*beta*SQR(2-q);
	else return 0;
}

// Method calculates rigid body mass center
// does not work with periodic boundary conditions
int MRigidBody::CalculateCenterOfMass(MParticleData& par)
{
	int i, count;

	R.clear();

	for (i=1,count=0; i<=par.GetSize(); ++i)
		if (par[i].mat==ID) { R += par[i].x; count++; }

	R *= 1./count;

	return OK;
}
