#include "MMaterial.h"
#include <string.h>
#include <math.h>
#include <stdio.h>

MMaterial::MMaterial()
{
  	mom.resize(3);
  	strinput.resize(48);
  	strinput2.resize(48);
  	eosinput.resize(24);

	model = 0;
	eos = 0;

	for (int i=0; i<2; i++) for (int j=0; j<12; j++)
		strcpy(head[i][j],"");

	eosinput.clear();
	strinput.clear();
	strinput2.clear();
}

MMaterial::~MMaterial()
{
}

MMaterial&  MMaterial::operator=(const MMaterial &mat)
{
  	n = mat.n;
  	model = mat.model;
  	eos = mat.eos;
  	visc_type = mat.visc_type;
  	rho = mat.rho;
  	g = mat.g;
  	av_l = mat.av_l;
  	av_q = mat.av_q;
  	surfaceTensionCoeff = mat.surfaceTensionCoeff;
  	mass = mat.mass;
  	h = mat.h;
  	rho_min = mat.rho_min;
  	rho_max = mat.rho_max;
  	ke = mat.ke;
  	ie = mat.ie;
  	xcf = mat.xcf;
  	ycf = mat.ycf;
  	zcf = mat.zcf;

  	// Vectors
  	mom.resize(3);
  	strinput.resize(48);
  	strinput2.resize(48);
  	eosinput.resize(24);
  	mom = mat.mom;
	eosinput = mat.eosinput;
	strinput = mat.strinput;
	strinput2 = mat.strinput2;

	for (int i=0; i<2; i++) for (int j=0; j<12; j++)
		strcpy(head[i][j],mat.head[i][j]);

	return *this;
}

MMaterial::MMaterial(const MMaterial &mat)
{
	this->operator=(mat);
}


void MMaterial::CheckDefaults()
{
	// set default artificial viscosity if it is zero in the input file
	if (visc_type == 0) visc_type = 1;
	if (av_l == 0.0) av_l = 0.06;
	if (av_q == 0.0) av_q = 1.5;

	// Minimum and maximum densities
	if (rho_min == 0.0) rho_min=1.0e-10;
	if (rho_max == 0.0) rho_max=1.0e+10;
}

// msg is error message
int MMaterial::CheckMaterialModel(std::string &msg)
{
	char tmp[300];

	switch (model)
	{
		// These models do not require EOS, error message if supplied
		case 1: case 3:
			sprintf(tmp, "\nEOS not permitted for material number %5d!\n", n);
			msg = tmp;
			return ERROR;
			break;
		//
		// 9 and 10 materials require EOS
		case 9: case 10:
			switch (eos)
			{
				// Supported EOS
				case 1: case 4: case 13: case 28: case 29: case 41:
				break;
				// Wrong EOS type
				default:
					sprintf(tmp, "Wrong EOS type for material %d !", n);
					msg = tmp;
					return ERROR;
			};break;
		// Wrong material number
		default:
			sprintf(tmp, "Wrong material model for material %d!", n);
			msg = tmp;
			return ERROR;
	} // material model switch

	return OK;
}
