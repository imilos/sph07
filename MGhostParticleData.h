#ifndef MGHOSTPARTICLEDATA_H_
#define MGHOSTPARTICLEDATA_H_

#include "global.h"
#include "MGhostParticle.h"
#include "MParticle.h"
#include <vector>

class MGhostParticleData
{
public:
	MGhostParticleData();
	virtual ~MGhostParticleData();

public:
	int AddParticle(const MGhostParticle &gpar);
	MGhostParticle& GetParticle(int ID);
	int CreateGhostFromReal(MParticle &realpar, int id, RealVector &addVec, int &llgrid_input);
	int RemoveParticle(int ID);
	MGhostParticle& operator[](int ID);
	int DeleteContents();
	int GetSize();

public:	
	int sph_ngp;
	
private:
	std::vector<MGhostParticle> m_GhostParList;
};

#endif /*MGHOSTPARTICLEDATA_H_*/
