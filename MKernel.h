#ifndef MKERNEL_H_
#define MKERNEL_H_

#include "global.h"

class MKernel
{
public:
	MKernel(int ndim, real havg);
	virtual ~MKernel();
	
	virtual real Kernel(RealVector &xi, RealVector &xj) = 0;
	virtual void GradW(RealVector &dWdx, real &dWdr, RealVector xi, RealVector xj) = 0;
	
public:
	int m_ndim;
	real m_havg;
};

#endif /*MKERNEL_H_*/
