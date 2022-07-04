#ifndef MKERNELBSPLINE_H_
#define MKERNELBSPLINE_H_

#include "global.h"
#include "MKernel.h"

class MKernelBSpline : public MKernel
{
public:
	MKernelBSpline(int ndim, real havg);
	virtual ~MKernelBSpline();

	real Kernel(RealVector &xi, RealVector &xj);
	void GradW(RealVector &dWdx, real &dWdr, RealVector xi, RealVector xj);
	// static functions
	real static Kernel(RealVector &xi, RealVector &xj, int ndim, real havg);
	void static GradW(RealVector &dWdx, real &dWdr, RealVector &xi, RealVector &xj, int ndim, real havg);
};

#endif /*MKERNELBSPLINE_H_*/
