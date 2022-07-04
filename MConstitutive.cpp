#include "MConstitutive.h"
#include <stdio.h>

void MConstitutive::PrintScreenLog(const std::string &msg, FILE** logfile)
{
	printf("%s\n", msg.c_str());
	fprintf(*logfile, "%s\n", msg.c_str());
}

int MConstitutive::Constitutive(MParticle &p, MMaterial &mater, FILE** logfile, real gvdt)
{

char msg[300];

	switch (mater.model)
	{
		// Elastic
		case 1: break;
			F3DM1(p, mater, gvdt);
			AViscosity(p, mater);
			break;
		// Elastic-plastic
		case 3: break;
		// Fluid Hydrodynamic
		case 9:
			F3DM9(p, mater);
			AViscosity(p, mater);
			HieUpdate(p, gvdt);
			if (EquationOfState(p, mater, logfile) != OK) return ERROR;
			StressUp(p);
			break;
		// Elastic-plastic hydrodynamic
		case 10: break;

		default:
			sprintf(msg, "Material type %5d is not supported.", mater.model);
			PrintScreenLog(msg, logfile);
			return ERROR;
		}

	return OK;
}

int MConstitutive::F3DM9(MParticle &p, MMaterial &mater)
{
	real davg, xmu;
	RealMatrix davgMatrix;
	// Dynamic viscosity hardcoded, to be changed
	if (mater.eos==13)
		xmu = mater.strinput(1);
	else if (mater.eos==28 || mater.eos==29)
		xmu = mater.av_l*mater.rho;

	// Deviatoric part of stress tensor
	if ( fabs(xmu)<1.e-10 )
  		p.s.clear();
	else
	{
		davg = 1.0/3.0 * p.tracerod;
		davgMatrix.clear();
		for (int i=0; i<3; i++) davgMatrix(i,i)=davg;

		// Deviatoric part of stress tensor
		p.s = 2 * xmu * (p.rod - davgMatrix);
	}

	return OK;
}
/*
int MConstitutive::LIEUPD(MParticle &p, MMaterial &mater, real gvdt)
{

        real qpm = 0.5;
//??????????????????????????
        p.e = p.e +(p.einc - gvdt*qpm)*p.mass/p.rho 
	return OK;


integer:: nn,nn1,davg
real(kind=real_acc) :: qp,qpm
real(kind=real_acc) :: sss(3,3),ss(9),ppp(3,3),pp(9) 
!
sss = par(nn)%q(1:3,1:3)
ppp = par(nn)%rod(1:3,1:3)
pp  = pack(ppp,.true.)
ss  = reshape(source=sss,shape=shape(ss))
qpm = dot_product(ss,pp)
!	qp=sum(qpm)
!
par(nn)%e = par(nn)%e + (par(nn)%einc - mcm_dt*qpm) * par(nn)%mass/par(nn)%rho
!
end
*/
//}
int MConstitutive::F3DM1(MParticle &p, MMaterial &mater, real gvdt)
{
        real ym = 2.1e10; //mater.strinput(42);   //Young's modulus
       	real pr = 0.3; //mater.strinput(42);   //Poison's ratio
        real gg = 0.0; //mater.strinput(42);   // shear modulus
        gg = ym/(2*(1+pr));
        gvdt = 0.5;

	RealMatrix davgMatrix, p_incrMatrix;

        real gdt, gd2, blk, davg, p_incr;


		gdt = gvdt*gg;
		gd2 = 2.0*gdt;
		blk = gvdt*ym/(1.0-2.0*(pr));
                davg = 1.0/3.0 * p.tracerod;
		p_incr = blk*davg;
		for (int i=0; i<3; i++) 
		{
			davgMatrix(i,i)=davg;
			p_incrMatrix(i,i) = p_incr;
		}

	p.sigma = p.sigma + p_incrMatrix + gd2*(p.rod - davgMatrix); // + gv.sph_dt(p.sigma*p.spin + p.spin*p.sigma)
        return OK;
}

void MConstitutive::AViscosity(MParticle &p, MMaterial &mater)
{
	if (p.tracerod<0.0)
	{
		real q = mater.av_l * p.rho * p.c * p.h * fabs(p.tracerod) + // linear
			mater.av_q * p.rho * SQR(p.h*p.tracerod);   		     // quadratic

		p.q.clear();
		for (int i=0; i<3; i++) p.q(i,i)=q;
	}
}

int MConstitutive::EquationOfState(MParticle &p, MMaterial &mater, FILE **logfile)
{
	char msg[300];

	switch (mater.eos)
	{
		// EOSes to be implemented
		case 1: case 4: case 41:
			sprintf(msg, "EOS %5d not supported yet.", mater.eos);
			PrintScreenLog(msg, logfile);
			return ERROR;
			break;
		// Monaghan quasi-incompressible fluid
		case 28:
			EOS28(p, mater);
			break;
		// Morris quasi-incompressible fluid
		case 29:
			EOS29(p, mater);
			break;
		// Perfect gas
		case 13:
			EOS13(p, mater);
			break;
	}

	return OK;
}

int MConstitutive::EOS28(MParticle &p, MMaterial &mater)
{
	// ASSIGN EOS PARAMETERS PER MATERIAL
	real b      = mater.eosinput(0);
	real gamma  = mater.eosinput(1);

	// CALCULATE pressure and sound speed, Monaghan 94.
	p.p = b * ( pow(p.rho/p.rho0,gamma) - 1.0);
	p.c   = sqrt(b * gamma/p.rho0);

	return OK;
}

// EOS 29 is the equation introduced by Morris 97.
int MConstitutive::EOS29(MParticle &p, MMaterial &mater)
{
	// CALCULATE pressure and sound speed, Morris 97.
	p.c = mater.eosinput(0);
	p.p = SQR(p.c) * (p.rho - p.rho0);
	//p.p = SQR(p.c) * p.rho;

	return OK;
}

int MConstitutive::EOS13(MParticle &p, MMaterial &mater)
{
	real gamma_m1  = mater.eosinput(0)-1;

	// specific 'trial' value for internal energy
	real e1try = p.etry/p.mass;
	real dvol = p.mass/p.rho - p.mass/p.rhoold;
	real vol0 = p.mass/p.rho0;

	// CALCULATE p (according to J ANDERSON: Modern Compressible Flow
	//           James' notes)

	p.p = (p.rho*gamma_m1*e1try) / (1.0+0.5*p.rho*gamma_m1*dvol/p.mass);

	// CALCULATE e(i)
	p.e = p.etry - 0.5*dvol*p.p;

	/***************************************************************/

	return OK;
}

int MConstitutive::StressUp(MParticle &p)
{
	RealMatrix pressure;

	pressure.clear();

	for (int i=0; i<3; i++) pressure(i,i)=p.p;

	p.sigma = -pressure + p.s;


	return OK;
}

//    Purpose: Calculate trial value of e
int MConstitutive::HieUpdate(MParticle &p, real gvdt)
{
	real pold, othird=1./3, trace_q, volold, volnew, vavg, eincr, dvol;
	real mean_incr;
	RealMatrix de;
	int i, j, k, l;

	pold = -othird * (p.sigma(0,0) + p.sigma(1,1) + p.sigma(2,2));
	trace_q = othird*(p.q(0,0) + p.q(1,1) + p.q(2,2));
	volold = p.mass/p.rhoold;
	volnew = p.mass/p.rho;

	for(l=0; l<3; l++)
	for (k=0; k<3; k++)
		if (k==l)
			de(k,l) = p.rod(k,l)*(0.5*(p.sigma(k,l)+pold+p.s(k,l)) - p.q(k,l)+trace_q);
		else
			de(k,l) = p.rod(k,l)*(0.5*(p.sigma(k,l)+p.s(k,l)) - p.q(k,l));

	vavg = 0.5*(volnew + volold);

	for(l=0, eincr=0.0; l<3; l++)
	for (k=0; k<3; k++)
		eincr += de(k,l);

	eincr *= vavg;

	// calculate energy increment due to hydrostatic (mean) terms
	dvol = volnew - volold;
	mean_incr = dvol * (0.5*pold + trace_q);

	// calculate trial value of internal energy
	p.etry = p.e + gvdt*eincr - mean_incr;

	p.qold = p.q;

	return OK;
}

