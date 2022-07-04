#ifndef MHISTORYVARS_H_
#define MHISTORYVARS_H_

#include "global.h"

class MHistoryVars
{
public:
	MHistoryVars();
	virtual ~MHistoryVars();
public:
	real sph_thermale, sph_kinetice, sph_internale, sph_totale;	
};

#endif /*MHISTORYVARS_H_*/
