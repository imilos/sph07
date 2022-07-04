#ifndef MSIMULATIONINIT_H_
#define MSIMULATIONINIT_H_

#include "MSimulationData.h"
#include "MSolution.h"
#include "MOutput.h"

class MSimulationInit
{
// Member functions
protected:
	MSimulationInit(MSimulationData *simdata, MSolution *solution, MOutput *output);
	virtual ~MSimulationInit();

public:
	static MSimulationInit* Instance(MSimulationData *simdata, MSolution *solution, MOutput *output);
	void DestroyInstance();
private:
	static MSimulationInit *_instance;

public:
	int Initialize();
	int InitMass();
	int InitVariables();
	int InitG();
	int InitFirstTimeStep();
	int InitNeighbours();
	int InitBoundaryMultipliersAndAddtions();
	int InitRigidBodyData();

private:
	void PrintScreenLog(const std::string &msg);
	int PrintMass();
	int CalcMassFromTotal();
	int InitRho_P_C();
	int InitSigma();
	int InitH();
	int InitOld();
	int EOScalc(MParticle &part);

// Attributes
public:
	MSimulationData *m_simdata;
	MSolution *m_solution;
	MOutput *m_output;
};

#endif /*MSIMULATIONINIT_H_*/
