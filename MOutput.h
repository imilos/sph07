#ifndef MOUTPUT_H_
#define MOUTPUT_H_

#include "MSimulationData.h"

class MOutput
{
protected:
	MOutput(MSimulationData *simdata);
	virtual ~MOutput();

public:
	static MOutput* Instance(MSimulationData *simdata);
	void DestroyInstance();
private:
	static MOutput *_instance;

public:
	void PrintScreenLog(const std::string &msg);
	int InitStateOutput();
	int StateOutput();
	int WriteCase();
	int WriteCase1D();
	int WriteCase2D();
	int WriteCase3D();
	int WriteCaseGeo();
	void IncreaseStateCounter();

private:
	MSimulationData *m_simdata;
};

#endif /*MOUTPUT_H_*/
