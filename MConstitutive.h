#ifndef MCONSTITUTIVE_H_
#define MCONSTITUTIVE_H_

#include "global.h"
#include "MParticle.h"
#include "MMaterial.h"
#include <string>
#include "MSimulationData.h"
#include <vector>



class MConstitutive
{
public:
	MConstitutive() {};
	virtual ~MConstitutive() {};
	static void PrintScreenLog(const std::string &msg, FILE ** logfile);
	static int Constitutive(MParticle &p, MMaterial &mater, FILE **logfile, real gvdt);
	static int F3DM9(MParticle &p, MMaterial &mater);
	static int F3DM1(MParticle &p, MMaterial &mater, real);
	static void AViscosity(MParticle &p, MMaterial &mater);
	static int EquationOfState(MParticle &p, MMaterial &mater, FILE **logfile);
	static int HieUpdate(MParticle &p, real gvdt);
	static int EOS28(MParticle &p, MMaterial &mater);
	static int EOS29(MParticle &p, MMaterial &mater);
	static int EOS13(MParticle &p, MMaterial &mater);
	static int StressUp(MParticle &p);
	
	// Attributes
    protected:
	MSimulationData *m_simdata;
	
};


#endif /*MCONSTITUTIVE_H_*/
