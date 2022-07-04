#include "MRigidBodyData.h"

/**
 * Class MRigidBodyData is intended to contain all data about rigid bodies
 * So far, only one rigid body can be handled, planned to be multiple.
 * This structure is 1 based, same as particle and material.
 * Milos Ivanovic, September 2008.
 */

MRigidBodyData::MRigidBodyData() {
	// TODO Auto-generated constructor stub

}

MRigidBodyData::~MRigidBodyData() {
	// TODO Auto-generated destructor stub
}

int MRigidBodyData::AddRigidBody(const MRigidBody &rb)
{
	rigidBodyList.push_back(rb);
	return rigidBodyList.size()-1;
}

MRigidBody& MRigidBodyData::GetRigidBody(int ID)
{
	return rigidBodyList[ID-1];
}

MRigidBody& MRigidBodyData::operator[](int ID)
{
	return rigidBodyList[ID-1];
}

int MRigidBodyData::GetSize()
{
	return rigidBodyList.size();
}

int MRigidBodyData::DeleteContents()
{
	rigidBodyList.clear();
	return OK;
}
