#ifndef MRIGIDBODYDATA_H_
#define MRIGIDBODYDATA_H_

#include "MRigidBody.h"

/**
 * Class MRigidBodyData is intended to contain all data about rigid bodies
 * So far, only one rigid body can be handled, planned to be multiple.
 * This structure is 1 based, same as particle and material.
 * Milos Ivanovic, September 2008.
 */

class MRigidBodyData
{
// Attributes
public:
	std::vector<MRigidBody> rigidBodyList;
	int ID;
	real totalMass;
	real inertiaMoment;
	real particleMass;
	real deltap;

// Methods
public:
	MRigidBodyData();
	virtual ~MRigidBodyData();

	int AddRigidBody(const MRigidBody &rb);
	MRigidBody& GetRigidBody(int ID);
	MRigidBody& operator[](int ID);
	int DeleteContents();
	int GetSize();
};

#endif /* MRIGIDBODYDATA_H_ */
