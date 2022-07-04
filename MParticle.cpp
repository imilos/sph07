#include "MParticle.h"

MParticle::MParticle() : MCommonParticle()
{
	// delete particles linked list
	delpointer = 0;
	efps = 0.0;
	fail = 1.0;
	tracerod = 0.0;

  	epx.resize(3); epx.clear();
  	alfa.resize(3,3); alfa.clear();

  	// vectors
  	smooth_v.resize(3); smooth_v.clear();
	x.resize(3); x.clear();
	v.resize(3); v.clear();
	a.resize(3); a.clear();
  	bndnorm.resize(3); bndnorm.clear();
  	repulsion.resize(3); repulsion.clear();

  	// tensors
  	rod.resize(3,3); rod.clear();
  	spin.resize(3,3); spin.clear();
  	s.resize(3,3); s.clear();
  	qold.resize(3,3); qold.clear();
  	qq.resize(3,3); qq.clear();
  	qr.resize(3,3); qr.clear();
}


MParticle::~MParticle()
{
}

void MParticle::CalculateAabs(int ndim)
{
	aabs = 0.0;
	for (int j=0; j<ndim; j++) aabs+= SQR( a(j) );
	aabs = sqrt(aabs);
}

void MParticle::CalcTraceROD(int ndim)
{
	tracerod = 0.0;
	for (int j=0; j<ndim; j++) tracerod+= rod(j,j);
}

