#ifndef MLOADCURVEVARS_H_
#define MLOADCURVEVARS_H_

#include "global.h"
#include <vector>

class MLoadCurveVars
{
public:
	MLoadCurveVars();
	virtual ~MLoadCurveVars();
	
public:
	RealVector sph_npc, sph_plc;
};

#endif /*MLOADCURVEVARS_H_*/
