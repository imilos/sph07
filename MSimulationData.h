#ifndef MSIMULATIONDATA_H_
#define MSIMULATIONDATA_H_

#include "global.h"
#include "MParticleData.h"
#include "MGhostParticleData.h"
#include "MMaterialData.h"
#include "MGlobalVars.h"
#include "MHistoryVars.h"
#include "MOptionVars.h"
#include "MOutputVars.h"
#include "MNeighbourVars.h"
#include "MSymPerVars.h"
#include "MFileHandlingVars.h"
#include "MLoadCurveVars.h"
#include "MContactVars.h"
#include "MBaseAccelVars.h"
#include "MRigidBodyData.h"

class MSimulationData
{
protected:
	MSimulationData();
	virtual ~MSimulationData();

public:
	static MSimulationData* Instance();
	void DestroyInstance();
private:
	static MSimulationData *_instance;

public:
	// Constants
	real pi;
	// Global variables
	MGlobalVars m_globvars;
	// History variables
	MHistoryVars m_histvars;
	// Options
	MOptionVars m_optvars;
	// Output
	MOutputVars m_outvars;
	// Neighbours and Neighbour searching
	MNeighbourVars m_neighbourvars;
	// Symmetry + periodic planes
	MSymPerVars m_sympervars;
	// File handling
	MFileHandlingVars m_filevars;
	// Load Curves
	MLoadCurveVars m_loadcurvevars;
	// Contact Data
	MContactVars m_contactvars;
	// Base accelerations
	MBaseAccelVars m_baseaccelvars;
	// Rigid body data
	MRigidBodyData m_rigidbodydata;
	//
	//----------------------------------------------------------------------------
	// Dynamic data
	//
	MParticleData par;
	MGhostParticleData gpar;
	MMaterialData mat;
	//---------------------------------------------------------------------------
	// MPI variables
	int my_id, numprocs;

// Member Functions
public:
	int Startup(int argc, char **argv, bool &newproblem);
	void Shutdown(int exitcode);
	void Fatal(const std::string &message);
	int GetInput();
	int PassInputComments(std::string &ret);

protected:
	int ParseCommandLine(int argc, char **argv);
	void PrintVersion();
	int ControlInput(int version_number);
	int MaterialInput();
	int PrintMaterialProperties(MMaterial &tmp);
	int NodeInput();
	int VelocityInput();
	int BoundaryInput();
	int BaseAccelerationInput();
	int PlotInput();
	int LennardJonesInput();
	int RigidBodyInput();
};

#endif /*MSIMULATIONDATA_H_*/
