#include "MCommonParticle.h"

MCommonParticle::MCommonParticle()
{
	llpointer = LINKEDLISTEND;

	active = true;
	fromOtherProc = false;
	lennardJones = false;

	color = 1;
	N = 0;
	nnbr = 0;
	x.resize(3); x.clear();
	xzero.resize(3); xzero.clear();
	v.resize(3); v.clear();
	p = 0.0;
	sigma.resize(3,3); sigma.clear();
	q.resize(3,3); q.clear();
}

MCommonParticle::~MCommonParticle()
{
}

void MCommonParticle::CalculateVabs(int ndim)
{
	vabs = 0.0;
	for (int j=0; j<ndim; j++) vabs+= SQR( v(j) );
	vabs = sqrt(vabs);
}

real MCommonParticle::Distance(const MCommonParticle &other, int ndim)
{
	real dist_sq = 0.0;
	int j;

	for (j=0; j<ndim; j++)
		dist_sq+= SQR(x(j)-other.x(j));

	return sqrt(dist_sq);
}

real MCommonParticle::DistanceSquared(const MCommonParticle &other, int ndim)
{
	real dist_sq = 0.0;
	int j;

	for (j=0; j<ndim; j++)
		dist_sq+= SQR(x(j)-other.x(j));

	return dist_sq;
}


int MCommonParticle::AddNeighbour(int id)
{
	m_nbrlist(nnbr++) = id;
	return OK;
}

int MCommonParticle::RemoveNeighbours()
{
	nnbr = 0;
	m_nbrlist.clear();
	return OK;
}

int MCommonParticle::ReserveMemForNbrs(int size)
{
	m_nbrlist.resize(size);
	return OK;
}
