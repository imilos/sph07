#ifndef MGHOSTPARTICLE_H_
#define MGHOSTPARTICLE_H_

#include "global.h"
#include "MCommonParticle.h"

class MGhostParticle : public MCommonParticle
{
public:
	MGhostParticle();
	virtual ~MGhostParticle();
	
public:
	int par;                         // Real particle that ghost is created from
};

#endif /*MGHOSTPARTICLE_H_*/
