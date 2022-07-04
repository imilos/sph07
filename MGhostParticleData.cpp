#include "MGhostParticleData.h"

MGhostParticleData::MGhostParticleData()
{
	sph_ngp = 0;
}

MGhostParticleData::~MGhostParticleData()
{
}

int MGhostParticleData::CreateGhostFromReal(MParticle &realpar, int id, RealVector &addVec, int &llgrid_cell)
{
	MGhostParticle ghostpar;

	ghostpar.par = id;
	ghostpar.mat = realpar.mat;
	ghostpar.mass = realpar.mass;
	ghostpar.h = realpar.h;
	ghostpar.hold = realpar.hold;
	ghostpar.rho = realpar.rho;
	ghostpar.color = realpar.color;

	ghostpar.x = realpar.x + addVec;
	ghostpar.xzero = realpar.xzero + addVec;
	ghostpar.v = realpar.v;
	//ghostpar.p = realpar.p;
	ghostpar.vabs = realpar.vabs;

	ghostpar.llpointer = llgrid_cell;
	llgrid_cell = - (++sph_ngp);

	m_GhostParList.push_back(ghostpar);

	return OK;
}

int MGhostParticleData::AddParticle(const MGhostParticle &gpar)
{
	m_GhostParList.push_back(gpar);
	return m_GhostParList.size()-1;
}

int MGhostParticleData::DeleteContents()
{
	sph_ngp = 0;
	m_GhostParList.clear();

	return OK;
}

int MGhostParticleData::GetSize()
{
	return m_GhostParList.size();
}

MGhostParticle& MGhostParticleData::GetParticle(int ID)
{
	return m_GhostParList[ID-1];
}

// IMPORTANT: The operator [] uses one based notation, not zero based
MGhostParticle& MGhostParticleData::operator[](int ID)
{
	return m_GhostParList[ID-1];
}
