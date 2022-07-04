#include "MKernelBSpline.h"
#include <math.h>

MKernelBSpline::MKernelBSpline(int ndim, real havg) : MKernel(ndim, havg)
{
}

MKernelBSpline::~MKernelBSpline()
{
}

real MKernelBSpline::Kernel(RealVector &xi, RealVector &xj)
{
	real const1, const2, xij[3], z, rad, W;
	int n;

	switch (m_ndim)
	{
		case 3:
			const1 = 1.0  / ( PI*CUBE(m_havg) );
			break;
		case 2:
			const1 = 10.0 / ( 7.0*PI*SQR(m_havg) );
			break;
		case 1:
			const1 = 2.0  / ( 3.0*m_havg );
			break;
	}

	const2 = 0.25 * const1;
	for (n=0, rad=0.0; n<m_ndim; n++)
	{
		xij[n] = xi(n) - xj(n);
		rad += SQR( xij[n] );
	}
	rad = sqrt(rad);
	z = rad/m_havg;

	if (z<1.0)
		W = const1 * ( ( (0.75 * z - 1.5) * z) * z + 1.0);
	else if (z<2.0)
		W = const2 * CUBE(2.0 - z);
	else
		W = 0.0;

	return W;
}

void MKernelBSpline::GradW(RealVector &dWdx, real &dWdr, RealVector xi, RealVector xj)
{
	real const1, const2, xij[3], z, rad, dW;
	int n;

	switch (m_ndim)
	{
		case 3:
			const1 = 1.0  / ( PI*CUBE(m_havg) );
			break;
		case 2:
			const1 = 10.0 / ( 7.0*PI*SQR(m_havg) );
			break;
		case 1:
			const1 = 2.0  / ( 3.0*m_havg );
			break;
	}

	const2 = 0.25 * const1;
	for (n=0, rad=0.0; n<m_ndim; n++)
	{
		xij[n] = xi(n) - xj(n);
		rad+= SQR( xij[n] );
	}
	rad = sqrt(rad);
	z = rad/m_havg;

	if (z<1.0)
		dW = const1 * ( (9.0/4.0 * z - 3.0) * z );
	else if (z<2.0)
		dW = -3.0 * const2 * SQR(2.0 - z);
	else
		dW = 0.0;

	for (n=0; n<m_ndim; n++)
	{
		if (rad==0.0)
		{
			dWdx(n) = 0.0;
			dWdr = 0.0;
		}
		else
		{
			dWdx(n) = dW*xij[n] / (rad*m_havg);
			dWdr = dW / m_havg;
		}
	}
}

real MKernelBSpline::Kernel(RealVector &xi, RealVector &xj, int ndim, real havg)
{
	real const1, const2, xij[3], z, rad, W;
	int n;
	int m_ndim = ndim;
	real m_havg = havg;

	switch (m_ndim)
	{
		case 3:
			const1 = 1.0  / ( PI*CUBE(m_havg) );
			break;
		case 2:
			const1 = 10.0 / ( 7.0*PI*SQR(m_havg) );
			break;
		case 1:
			const1 = 2.0  / ( 3.0*m_havg );
			break;
	}

	const2 = 0.25 * const1;
	for (n=0, rad=0.0; n<m_ndim; n++)
	{
		xij[n] = xi(n) - xj(n);
		rad += SQR( xij[n] );
	}
	rad = sqrt(rad);
	z = rad/m_havg;

	if (z<1.0)
		W = const1 * ( ( (0.75 * z - 1.5) * z) * z + 1.0);
	else if (z<2.0)
		W = const2 * CUBE(2.0 - z);
	else
		W = 0.0;

	return W;
}

void MKernelBSpline::GradW(RealVector &dWdx, real &dWdr, RealVector& xi, RealVector& xj, int ndim, real havg)
{
	real const1, const2, xij[3], z, rad, dW;
	int n;
	int m_ndim = ndim;
	real m_havg = havg;
	dWdx.clear();

	switch (m_ndim)
	{
		case 3:
			const1 = 1.0  / ( PI*CUBE(m_havg) );
			break;
		case 2:
			const1 = 10.0 / ( 7.0*PI*SQR(m_havg) );
			break;
		case 1:
			const1 = 2.0  / ( 3.0*m_havg );
			break;
	}

	const2 = 0.25 * const1;
	for (n=0, rad=0.0; n<m_ndim; n++)
	{
		xij[n] = xi(n) - xj(n);
		rad+= SQR( xij[n] );
	}
	rad = sqrt(rad);
	z = rad/m_havg;

	if (z<1.0)
		dW = const1 * ( (9.0/4.0 * z - 3.0) * z );
	else if (z<2.0)
		dW = -3.0 * const2 * SQR(2.0 - z);
	else
		dW = 0.0;

	for (n=0; n<m_ndim; n++)
	{
		if (rad==0.0)
		{
			dWdx(n) = 0.0;
			dWdr = 0.0;
		}
		else
		{
			dWdx(n) = dW*xij[n] / (rad*m_havg);
			dWdr = dW / m_havg;
		}
	}
}
