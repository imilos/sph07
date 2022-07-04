#ifndef MMATERIAL_H_
#define MMATERIAL_H_

#include "global.h"
#include <string>

class MMaterial
{
public:
	MMaterial();
	MMaterial(const MMaterial &mat);
	virtual ~MMaterial();
	MMaterial& operator=(const MMaterial &mat);
	void CheckDefaults();
	int CheckMaterialModel(std::string &msg);

public:
	int n;
	char head[2][12][6];             // Material title
  	int model;                       // Material model number
  	int eos;                         // Material equation of state model
  	int visc_type;                   // Material artificial viscosity type
  	real rho;                        // density
  	real g;                          // shear modulus
  	real av_l;                       // linear artificial viscosity coefficient
  	real av_q;                       // quadratic artificial viscosity coefficient
  	real surfaceTensionCoeff;		 // surface tension coefficient
  	real mass;                       // total material mass, used if sph_massopt = 0
  	real h;                          // smoothin length, used if sph_init_h_opt = 0
  	real rho_min;                    // minimum density limit factor
  	real rho_max;                    // maximum density limit factor
  	real ke;                         // Total kinetic energy of material
  	real ie;                         // Total internal energy of material
  	real xcf;						 // Contact force in X-direction
  	real ycf;						 // Contact force in Y-direction
  	real zcf;						 // Contact force in Z-direction
  	RealVectorL mom;     		     // Momentum vector
  	RealVectorL strinput;   	     // Constitutive coefficients
  	RealVectorL strinput2;  	     // Additional constitutive coefficients
  	RealVectorL eosinput;   	     // EoS coefficients
};

#endif /*MMATERIAL_H_*/
