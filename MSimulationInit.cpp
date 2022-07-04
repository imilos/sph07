#include "MSimulationInit.h"
#include <math.h>
#include <stdio.h>

MSimulationInit::MSimulationInit(MSimulationData *simdata, MSolution *solution, MOutput *output)
{
	m_simdata = simdata;
	m_solution = solution;
	m_output = output;
}

MSimulationInit::~MSimulationInit()
{
}

MSimulationInit* MSimulationInit::_instance=0;
MSimulationInit *MSimulationInit::Instance(MSimulationData *simdata, MSolution *solution, MOutput *output)
{
	if (_instance == 0)
		_instance = new MSimulationInit(simdata, solution, output);
	return _instance;
}

void MSimulationInit::DestroyInstance()
{
	delete _instance;
	_instance=0;
}

// Initialises stress, calculates initial timestep, particle masses and everything needed
int MSimulationInit::Initialize()
{
	MSimulationData &sd = *m_simdata;
	MGlobalVars &gv = sd.m_globvars;
	MParticleData &par = sd.par;
	int i;

	PrintScreenLog("\nProblem Initialisation Begun");

	// Apply displacement boundary conditions to initial velocity
	// to be implemented for periodic boundary conditions
	//InitBoundConditions();
	// Include boundary multipliers
	if (InitBoundaryMultipliersAndAddtions()!=OK) return ERROR;
	// Calculate particle mass
	if (InitMass()!=OK) return ERROR;
	// Initialize history and variables
	if (InitVariables()!=OK) return ERROR;
	// Calculate shear modulus
	if (InitG()!=OK) return ERROR;
	// Init density, pressure and speed of sound
	if (InitRho_P_C()!=OK) return ERROR;
	// Initialize stress tensor and rate of deformation tensor
	if (InitSigma()!=OK) return ERROR;
	// Initial smoothing lengths
	if (InitH()!=OK) return ERROR;
	// Initial smoothing lengths
	if (InitRigidBodyData()!=OK) return ERROR;
	// Initial timestep and artificial viscosity calculation
	gv.sph_timestep = 1;

	for (i=gv.sph_ssp; i<=gv.sph_esp; i++)
	{
		// Trace of the rate of deformation
		par[i].CalcTraceROD(gv.sph_ndim);
		// Artificial viscosity calculation
		m_solution->AViscosity(par[i]);
	}

	// Get neighbours and boundary particles
	InitNeighbours();
	m_solution->Neighbours();
	m_solution->LinkedListNeighbours();
	// Calculate critical timestep
	if (m_solution->CalcCriticalTimestep() != OK) return ERROR;
	// Set first time step
	if (InitFirstTimeStep()!=OK) return ERROR;

	// Calculate total energy, to be implemented
	//m_solution->CalcTotalEnergy();
	// Initialize rhoold and qold[3][3]
	if (InitOld()!=OK) return ERROR;

	gv.sph_next_restart = gv.sph_restart_interval;
	gv.sph_next_run_restart = gv.sph_run_restart;

	// Initialise state output variables and plot first state output file
	if (m_output->InitStateOutput()!=OK) return ERROR;
	m_output->StateOutput();
	// Initialise time history output files and write first data line
	// to be implemented in the future
	//m_output->InitTimeHistory();
	//m_output->TimeHistory();

	// Write status report to screen and log file
	PrintScreenLog("\nProblem Initialisation Complete");

	return OK;
}

void MSimulationInit::PrintScreenLog(const std::string &msg)
{
	printf("%s\n", msg.c_str());
	fprintf(m_simdata->m_filevars.f_logfile, "%s\n", msg.c_str());
}

int MSimulationInit::InitMass()
{
	switch (m_simdata->m_optvars.sph_massopt)
	{
		// Mass defined as total material mass
		case 0:
			 CalcMassFromTotal();
			 PrintMass();
			 break;
		// Mass associated to each particle in the input file
		case 1:
			PrintMass();
			break;
		default:
			PrintScreenLog("Error: Mass option not recognized.");
			return ERROR;
	}

	return OK;
}

int MSimulationInit::PrintMass()
{
	int i, j, count;
	real mass;
	MSimulationData &sd = *m_simdata;
	MGlobalVars &gv = sd.m_globvars;
	MFileHandlingVars &fv = sd.m_filevars;

	fprintf(fv.f_logfile, "\tPARTICLE MASS CALCULATION\n");
	fprintf(fv.f_logfile, "\tMaterial     Total Mass  No. particles\n");

	for ( i=1; i<=gv.sph_nummat; i++)
	{
		for (j=gv.sph_svp,count=0,mass=0.; j<=gv.sph_evp; j++)
			if ( sd.par[j].mat==i )
			{
				count++;
				mass += sd.par[j].mass;
			}
		fprintf(fv.f_logfile, "\t%8d   %12.5e       %8d\n", i, mass, count);
	}

	return OK;
}

int MSimulationInit::CalcMassFromTotal()
{
	int i;
	MSimulationData &sd = *m_simdata;
	MGlobalVars &gv = sd.m_globvars;
	MParticleData &par = sd.par;
	std::vector<int> par_per_mat;

	// IMPORTANT: Ignore index zero - one based !!!
	par_per_mat.push_back(0);

	// Initialize par_per_mat
	for (i=1; i<=gv.sph_nummat; i++) par_per_mat.push_back(0);

	// Calculate number of particles per material
	for (i=gv.sph_svp; i<=gv.sph_evp; i++)
		if (!par[i].lennardJones) par_per_mat[par[i].mat]++;

	// Calculate mass of each particle
	for (i=gv.sph_svp; i<=gv.sph_evp; i++)
		par[i].mass = sd.mat[par[i].mat].mass / par_per_mat[par[i].mat];
	// To be implemented: axis symmetry

	return OK;
}

// Initialize variable data
int MSimulationInit::InitVariables()
{
	int i;
	char msg[300];
	MSimulationData &sd = *m_simdata;
	MGlobalVars &gv = sd.m_globvars;
	MParticleData &par = sd.par;
	MMaterialData &mat = sd.mat;

	// Initial particle coordinates
	par.InitXzero();

	// Initialise stress point variables, i.e. cut-off pressure
	// ssp = svp, esp = evp in case of collocated particles
	for (i=gv.sph_ssp; i<=gv.sph_esp; i++)
	{
		MMaterial &mater = mat[par[i].mat];

		switch ( mater.model )
		{
			case 1: case 3:
				par[i].pcut = -9.9e20;
				break;
			case 9:
				par[i].pcut = mater.strinput(0);
				break;
			case 10:
				par[i].pcut = mater.strinput(3);
				break;
			default:
				sprintf(msg,"\nMaterial number %5d does not exist.\n", mater.model);
				PrintScreenLog(msg);
				return ERROR;
		}
	}

	// Initialize velocity point variables
	par.CalculateVabs(gv.sph_ndim);

	// Find maximum and minimum values of the coordinates for linked list underlying grid
	for (i=1; i<=par.GetSize(); i++)
		for (int j=0; j<gv.sph_ndim; j++)
		{
			gv.sph_coord_minmax(0,j) = MIN(gv.sph_coord_minmax(0,j), par[i].x(j));
			gv.sph_coord_minmax(1,j) = MAX(gv.sph_coord_minmax(1,j), par[i].x(j));
		}

	// Initialize unused particles, mat==0 => not drawn
	if (gv.sph_max_np>gv.sph_np) par.InitializeUnused(gv.sph_np);

	return OK;
}

// Calculate or assign shear modulus for each material
int MSimulationInit::InitG()
{
	MSimulationData &sd = *m_simdata;
	MParticleData &par = sd.par;
	MMaterialData &mat = sd.mat;
	int i;
	char msg[300];

	for (i=1; i<=mat.GetSize(); i++)
		switch (mat[i].model)
		{
			// Elastic, elasto-plastic
			case 1: case 3:
				mat[i].g = mat[i].strinput(0) / ( 2.0*(1.0+mat[i].strinput(5)) );
				break;
			// Fluid
			case 9:
				mat[i].g = 0.0;
				break;
			// Hydro-dynamic
			case 10:
				mat[i].g = mat[i].strinput(0);
				break;
			default:
				PrintScreenLog("Error in initialization of shear modulus:");
				sprintf(msg,"\nMaterial number %5d does not exist.\n", mat[par[i].mat].model);
				PrintScreenLog(msg);
				return ERROR;
		}

	return OK;
}

int MSimulationInit::InitRho_P_C()
{
	MSimulationData &sd = *m_simdata;
	MParticleData &par = sd.par;
	MMaterialData &mat = sd.mat;
	MGlobalVars &gv = sd.m_globvars;
	int i;
	real pmax;
	char msg[300];

	for (i=gv.sph_ssp; i<=gv.sph_esp; i++)
	{
		par[i].rho0 = mat[par[i].mat].rho;

		switch (mat[par[i].mat].model)
		{
			// Elastic and elastic-plastic
			case 1: case 3:
				if (sd.m_optvars.sph_init_rhoe!=1)
					par[i].rho = par[i].rho0;
				par[i].p = 0.0;
				par[i].c = sqrt(mat[par[i].mat].strinput(1)/par[i].rho);
				break;
			// Fluid material
			case 9:
				if (EOScalc(par[i])!=OK) return ERROR;
				// Pressure cut-off
				pmax = mat[par[i].mat].strinput(0);
				if (par[i].p > pmax) par[i].p = pmax;
				break;
			// Hydro-dynamic
			case 10:
				if (EOScalc(par[i])!=OK) return ERROR;
				break;
			// Error for wrong material
			default:
				PrintScreenLog("Error in initialization of rho, p & c:");
				sprintf(msg,"\nMaterial number %5d does not exist.\n", mat[par[i].mat].model);
				PrintScreenLog(msg);
				return ERROR;
		}
	}

	return OK;
}

int MSimulationInit::InitSigma()
{
	MSimulationData &sd = *m_simdata;
	MParticleData &par = sd.par;
	MGlobalVars &gv = sd.m_globvars;

	for (int i=gv.sph_ssp; i<gv.sph_esp; i++)
	{
		par[i].sigma(0,0) = -par[i].p;
		par[i].sigma(1,1) = -par[i].p;
		par[i].sigma(2,2) = -par[i].p;
	}

	return OK;
}

int MSimulationInit::InitH()
{
	MSimulationData &sd = *m_simdata;
	MParticleData &par = sd.par;
	MMaterialData &mat = sd.mat;
	MNeighbourVars &nv = sd.m_neighbourvars;
	MOptionVars &ov = sd.m_optvars;
	int i;

	switch ( ov.sph_init_h_opt )
	{
		// h defined with the material
		case 0:
			for (i=1; i<=par.GetSize(); i++) par[i].h = mat[par[i].mat].h;
			// Initialize hmax to maximum
			for (i=1; i<=mat.GetSize(); i++)
				if (nv.sph_hmax<mat[i].h) nv.sph_hmax = mat[i].h;
			break;
		// h defined for each particle separately
		case 1:
			for (i=1; i<=par.GetSize(); i++) par[i].h = mat[par[i].mat].h;
				if (nv.sph_hmax<par[i].h) nv.sph_hmax = par[i].h;
			break;
		default:
			PrintScreenLog("Error in initialization of h, option not recognized.");
			return ERROR;
	}

	// Set history variables
	for (i=1; i<=par.GetSize(); i++)
	{
		par[i].h0 = par[i].h;
		par[i].hold = par[i].h;
	}

	return OK;
}

int MSimulationInit::InitOld()
{
	MSimulationData &sd = *m_simdata;
	MParticleData &par = sd.par;
	MGlobalVars &gv = sd.m_globvars;
	int i,j,k;

	for (i=gv.sph_ssp; i<=gv.sph_esp; i++)
	{
 		par[i].rhoold = par[i].rho;
 		for (j=0; j<3; j++)
 			for (k=0; k<3; k++) par[i].qold(j,k) = par[i].q(j,k);
	}

	return OK;
}

int MSimulationInit::EOScalc(MParticle &part)
{
	MSimulationData &sd = *m_simdata;
	MMaterial &mater = sd.mat[part.mat];
	MOptionVars &ov = sd.m_optvars;
	real B, gamma, e1;
	char msg[300];

	// Only models 28, 29  and 13 implemented, others to be implemented
	switch (mater.eos)
	{
 		// Monaghan incompressible fluid
 		case 28:
 			if (mater.eosinput(3)==0.0) mater.eosinput(3)=1.0;
 			if (ov.sph_init_rhoe!=1)
 				part.rho = part.rho0 / mater.eosinput(3);

 			B = mater.eosinput(0);
 			gamma = mater.eosinput(1);
 			// Calculate p and c
 			part.p = B * ( pow(part.rho/part.rho0, gamma) - 1.0);
 			part.c = sqrt(B*gamma/part.rho0);
 			break;
 	 	// Morris incompressible fluid
 	 	case 29:
 	 		if (mater.eosinput(3)==0.0) mater.eosinput(3)=1.0;
 	 		if (ov.sph_init_rhoe!=1)
 	 			part.rho = part.rho0 / mater.eosinput(3);
 	 		part.c = mater.eosinput(0);
 	 		// Calculate p
 	 		part.p = part.c * (part.rho-part.rho0);
 	 		break;
 	 	// Perfect gas
 	 	case 13:
 	 		if (ov.sph_init_rhoe!=1)
 	 		{
 	 			part.rho = part.rho0 / mater.eosinput(3);
 	 			part.e   = mater.eosinput(2) * part.mass / part.rho;
 	 		}
 	 		// Calculate p and c
 	 		gamma = mater.eosinput(0);
 	 		// EoS requires internal energy per unit mass
 	 		e1 = part.e/part.mass;
 	 		part.p = part.rho * (gamma-1) * e1;
 	 		part.c = sqrt( (gamma-1)*gamma*e1 );
 	 		break;

 		// Error if EOS type not recognizable
 		default:
			PrintScreenLog("Error in initialization of EOS:");
			sprintf(msg,"\nEOS type %5d does not exist.\n", mater.eos);
			PrintScreenLog(msg);
			return ERROR;
	}

	return OK;
}

int MSimulationInit::InitFirstTimeStep()
{
	MSimulationData &sd = *m_simdata;
	MGlobalVars &gv = sd.m_globvars;
	char msg[300];

	if ( fabs(gv.sph_itss) > 1.e-10 )
	{
		// Select minimal
		gv.sph_dt = MIN(gv.sph_itss, gv.sph_critts*gv.sph_tssf);

		if ( gv.sph_itss > gv.sph_critts*gv.sph_tssf )
		{
			PrintScreenLog("WARNING: Initial time step size, greater than crit. size*safety_factor!");
			sprintf(msg,"Initial %11.3e, Critical %11.3e, Safety factor %11.3e.", gv.sph_itss, gv.sph_critts, gv.sph_tssf);
			PrintScreenLog(msg);
			//return ERROR;
		}
	}
	else
		gv.sph_dt = gv.sph_critts*gv.sph_tssf;

	return OK;
}

int MSimulationInit::InitNeighbours()
{
	MSimulationData &sd = *m_simdata;
	MGlobalVars &gv = sd.m_globvars;
	MParticleData &par = sd.par;
	//MSymPerVars &spv = sd.m_sympervars;

	// Symmetry and periodic planes, to be implemented

	par.InitNeighbours(gv.sph_maxnbr);
	//spv.InitGhostNbrList(gv.sph_ngp, spv.sph_g_maxnbr);

	return OK;
}

int MSimulationInit::InitBoundaryMultipliersAndAddtions()
{
	MSimulationData &sd = *m_simdata;
	MSymPerVars &spv = sd.m_sympervars;
	IntVector count_vec(3), cell(3);

	// Boundary additions and multipliers setup
	for (int i=0; i<3; i++)
		for (int j=0; j<3; j++)
		{
			spv.sph_mincoord_mult(i,j) = 1.0; spv.sph_mincoord_add(i,j) = 0.0;
			spv.sph_maxcoord_mult(i,j) = 1.0; spv.sph_maxcoord_add(i,j) = 0.0;
		}

	for (int m=0; m<3; m++)
	{
		switch ( spv.sph_boundary_code(0,m) )
		{
			// No boundary constraints
			case 0: break;
			// Symmetry, to be implemented
			case 1:
				PrintScreenLog("ERROR: Symmetry boundary conditions not implemented yet.");
				return ERROR;
			// Periodic boundary conditions, Xmax also present (required)
			case 2:
				spv.sph_mincoord_add(m,m) =  spv.sph_boundary_x(1,m)-spv.sph_boundary_x(0,m);
				spv.sph_maxcoord_add(m,m) = -spv.sph_boundary_x(1,m)+spv.sph_boundary_x(0,m);
				break;
			default:
				PrintScreenLog("ERROR: Boundary conditions code wrong.");
				return ERROR;
		}
	}

	return OK;
}

int MSimulationInit::InitRigidBodyData()
{
	MRigidBodyData &rbd = m_simdata->m_rigidbodydata;
	rbd.AddRigidBody( MRigidBody(rbd.ID, rbd.totalMass, rbd.inertiaMoment, rbd.particleMass, rbd.deltap) );
	return OK;
}
