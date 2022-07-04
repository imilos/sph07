#ifndef MPARTICLEDATA_H_
#define MPARTICLEDATA_H_

#include "global.h"
#include "MParticle.h"
#include <vector>

class MParticleData
{
public:
	MParticleData();
	virtual ~MParticleData();

public:
	int AddParticle(const MParticle &par);
	int RemoveParticle(int ID);
	MParticle& GetParticle(int ID);
	MParticle& operator[](int ID);
	int DeleteContents();
	int GetSize();
	void CalculateVabs(int ndim);
	void InitializeUnused(int np);
	void InitXzero();
	int InitNeighbours(int maxnbr);

private:
	std::vector<MParticle> m_ParList;
};

#endif /*MPARTICLEDATA_H_*/
