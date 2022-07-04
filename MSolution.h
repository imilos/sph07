#ifndef MSOLUTION_H_
#define MSOLUTION_H_

#include "MSimulationData.h"
#include "MOutput.h"
#include "MKernelBSpline.h"
#include "MConstitutive.h"
#include <vector>
#include <string>
#include "MRigidBody.h"

class MSolution
{
protected:
	MSolution(MSimulationData *simdata, MOutput *output);
	virtual ~MSolution();

public:
	static MSolution* Instance(MSimulationData *simdata, MOutput *output);
	void DestroyInstance();
private:
	static MSolution *_instance;

public:
	int Solution();
	void AViscosity(MParticle &part);
	int CalcCriticalTimestep();
	void CalcTotalEnergy();
	int Neighbours();
	int SetupLinkedList();
	int LinkedListNeighbours();
	int FindNewMaxNeighbour();
	int UpdateH();
	void PrintScreenLog(const std::string &msg);
	int Strain();
	int ClassicStrain();
	MKernel *GenerateKernel(real havg);
	int Constitutive();
	int RhoUpdate();
	int Momentum();
	int ClassicAcceleration();
	int MorrisAcceleration();
	int LennardJonesAcceleration();
	int UpdateVelocity();
	int BoundCond();
	int ConstrainQuantity(RealVector &quantity, int code);
	int MoveParticles();

	int GhostSetup();
	int GenerateGhostParticles(IntVector &down_limit, IntVector &up_limit, RealVector &addVec);

	// Surface tension handling
	int CalculateSurfaceNormals();

	// Rigid body methods, names of the methods start with 'RigidBody'
	int RigidBodyForce();
	int RigidBodyAcceleration();
	int RigidBodyVelocity();
	inline bool InRigidBody(const MCommonParticle &p)
		{ return m_simdata->m_rigidbodydata[1].ID==p.mat; }

	// Temporary functions, hardwired stuff
	int PlotOutput(const std::string &filename);
	int ConstrainVelocityCylinder(MParticle &p);
	int PlotDroppingCylinder(const std::string &filename);

	// MPI variables and methods
	int InitializeMpiProcess();
	int HandleMpiMemory(bool reserving);
	int TransferParticlesWithPeriodicProcess();
	int ExchangeData(int phase);
	int TransferParticles(int phase);
	int TransferParticlesToRoot(bool rigidBodyParticles);
	int PackMpiQuantities(bool packing, int &count, MParticle &part);
	int TransferMpiBuffer(bool sending, int &count, int proc, int tag);
	real *mpiExchangeTable;
	IntVectorL box_procmin;
	IntVectorL box_procmax;

// Attributes
protected:
	MSimulationData *m_simdata;
	MOutput *m_output;
};

#endif /*MSOLUTION_H_*/
