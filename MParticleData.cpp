#include "MParticleData.h"

MParticleData::MParticleData()
{
}

MParticleData::~MParticleData()
{
}

int MParticleData::AddParticle(const MParticle &par)
{
	m_ParList.push_back(par);
	return m_ParList.size()-1;
}

int MParticleData::RemoveParticle(int ID)
{
	m_ParList[ID-1].active = false;
	return OK;
}

int MParticleData::GetSize()
{
	return m_ParList.size();
}

MParticle& MParticleData::GetParticle(int ID)
{
	return m_ParList[ID-1];
}

// IMPORTANT: The operator [] uses one based notation, not zero based

MParticle& MParticleData::operator[](int ID)
{
	return m_ParList[ID-1];
}

int MParticleData::DeleteContents()
{
	m_ParList.clear();
	return OK;
}

void MParticleData::CalculateVabs(int ndim)
{
	for (int i=0; i<(int)m_ParList.size(); i++)
		m_ParList[i].CalculateVabs(ndim);
}

void MParticleData::InitializeUnused(int np)
{
	int i;
	for (i=np+1; i<=GetSize(); i++)
	{
		(*this)[i].mat = 0;
		(*this)[i].active = false;
		(*this)[i].delpointer = i+1;
		(*this)[i].h = (*this)[1].h;
		(*this)[i].x = (*this)[1].x;
	}
	(*this)[ GetSize() ].delpointer = 0;
}

void MParticleData::InitXzero()
{
	for (int i=0; i<(int)m_ParList.size(); i++)
		m_ParList[i].xzero = m_ParList[i].x;
}

int MParticleData::InitNeighbours(int maxnbr)
{
	// Vector of vectors, one dimension is number of particles
	for (int i=0; i<(int)m_ParList.size(); i++)
	{
		m_ParList[i].nnbr = 0;
		m_ParList[i].g_nnbr = 0;

		m_ParList[i].RemoveNeighbours();
		m_ParList[i].ReserveMemForNbrs(maxnbr);
	}

	return OK;
}
