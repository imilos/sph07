#ifndef MBASEACCELVARS_H_
#define MBASEACCELVARS_H_

#include "global.h"

class MBaseAccelVars
{
public:
	MBaseAccelVars();
	virtual ~MBaseAccelVars();

public:
	bool sph_baseaccel;
	RealVector sph_base_a;
};

#endif /*MBASEACCELVARS_H_*/
