#ifndef MCONTACTVARS_H_
#define MCONTACTVARS_H_

#include "global.h"
#include <vector>

class MContactVars
{
public:
	MContactVars();
	virtual ~MContactVars();

public:
	std::vector<real> *sph_k_cont, *sph_n_cont;
};

#endif /*MCONTACTVARS_H_*/
