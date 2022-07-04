#include "MSolution.h"
#include <math.h>
#include <mpi.h>

MSolution::MSolution(MSimulationData *simdata, MOutput *output)
{
	m_simdata = simdata;
	m_output = output;
}

MSolution::~MSolution()
{
	HandleMpiMemory(false);
}

MSolution* MSolution::_instance=0;
MSolution *MSolution::Instance(MSimulationData *simdata, MOutput *output)
{
	if (_instance == 0)
		_instance = new MSolution(simdata, output);
	return _instance;
}

void MSolution::DestroyInstance()
{
	delete _instance;
	_instance=0;
}
void MSolution::AViscosity(MParticle &part)
{
	MConstitutive::AViscosity(part, m_simdata->mat[ part.mat ] );
}

int MSolution::CalcCriticalTimestep()
{
	MSimulationData &sd = *m_simdata;
	MParticleData &par = sd.par;
	MGlobalVars &gv = sd.m_globvars;
	MOptionVars &ov = sd.m_optvars;
	int i;
	real viscosity, surfaceTensionCoeff;

	for (i=gv.sph_ssp, gv.sph_critts = 1.e20; i<=gv.sph_esp; i++)
	{
		MParticle &part = par[i];

		switch (ov.sph_tcrit_opt)
		{
			// Option 0: DYNA 3D formula
			case 0:
				PrintScreenLog("ERROR: Critical time step option not implemented yet.\n");
				return ERROR;
				break;
			// Option 1: Use h as length
			case 1:
				part.critts = part.h / (part.c + part.vabs);
				if (gv.sph_timestep==1) part.critts = 1.e-09;
				break;
			// Option 2: Use minimum interparticle distance as length
			case 2:
				part.critts = part.mindist / (part.c + part.vabs);
				if (gv.sph_timestep==1) part.critts = 1.e-09;
				break;
			// Option 3: Morris et al 1997. & 2000. (surface tension) time step calculation
			case 3:
				viscosity = sd.mat[part.mat].av_l;
				surfaceTensionCoeff = sd.mat[part.mat].surfaceTensionCoeff;
				part.critts = MIN( 0.25*part.h/part.c, 0.125*SQR(part.h)/viscosity );
				part.CalculateAabs(gv.sph_ndim);
				part.critts = MIN( part.critts, 0.25*sqrt(part.h/part.aabs) );
				// Surface tension CFL condition, according to Morris 2000.
				if (ov.sph_surface_tension)
					part.critts = MIN( part.critts, 0.25*sqrt(part.rho*CUBE(part.h)/2/PI/surfaceTensionCoeff) );
				break;
			// Error message
			default:
				PrintScreenLog("ERROR: Wrong critical time step option.\n");
				break;
		}

		// Check if particle timestep is less than current
		if (part.active && !part.lennardJones && !InRigidBody(part) && part.critts<gv.sph_critts)
			gv.sph_critts = part.critts;
	}

    // Gather minimal time step results for all processes
	double minimum_timestep;
	MPI_Allreduce(&gv.sph_critts, &minimum_timestep, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
	//MPI_Reduce(&gv.sph_critts, &minimum_timestep, 1, MPI_DOUBLE, MPI_MIN, ROOT_PROCESS, MPI_COMM_WORLD);
	//MPI_Bcast(&minimum_timestep, 1, MPI_DOUBLE, ROOT_PROCESS, MPI_COMM_WORLD);
	gv.sph_critts = minimum_timestep;

	return OK;
}

int MSolution::Neighbours()
{
	MSimulationData &sd = *m_simdata;
	MGlobalVars &gv = sd.m_globvars;
	MOptionVars &ov = sd.m_optvars;

	// 0, 1, 2 = Conventional SPH
	// 3    = Total Lagrangian
	switch (gv.sph_disctype)
	{
		case 0: case 1: case 2: case 3:
			if (SetupLinkedList()!=OK) return ERROR;
			// calculate coordinates and values of ghost particles
			//if (LinkedListNeighbours()!=OK) return ERROR;
			break;
	}

	return OK;
}

int MSolution::SetupLinkedList()
{
	MSimulationData &sd = *m_simdata;
	MNeighbourVars &nv = sd.m_neighbourvars;
	MGlobalVars &gv = sd.m_globvars;
	MSymPerVars &spv = sd.m_sympervars;
	MParticleData &par = sd.par;
	IntVector boxcrd(3);
	int i;

	boxcrd.clear();
	// Init underlying grid of size factor*h
	nv.InitGrid(gv.sph_ndim, spv.sph_boundary_type, gv.sph_coord_minmax, spv.sph_boundary_x);
	// Initialize MPI process grid values
	if (InitializeMpiProcess()!=OK) return ERROR;
	// Now set up list, particle loop puts particle to cells
	for (i=1; i<=gv.sph_np; i++)
	{
		// Only consider active particles
		if ( !par[i].ToCompute() ) { par[i].llpointer=LINKEDLISTEND; continue; }
		nv.GetBoxCoords(par[i].x, boxcrd);

		// Set inactive particles which are not interesting for this process
		if ( gv.sph_timestep==1 && (boxcrd(0)<box_procmin(sd.my_id) || boxcrd(0)>box_procmax(sd.my_id)) )
		{	par[i].active = false;	continue;	}

		// Each particle llpointer references to previous particle in that box
		par[i].llpointer = nv.GetGridValue( boxcrd(0), boxcrd(1), boxcrd(2) );
		nv.SetGridValue( boxcrd(0), boxcrd(1), boxcrd(2), i);
	}

	return OK;
}

int MSolution::LinkedListNeighbours()
{
	MSimulationData &sd = *m_simdata;
	MNeighbourVars &nv = sd.m_neighbourvars;
	MGlobalVars &gv = sd.m_globvars;
	MParticleData &par = sd.par;
	MGhostParticleData &gpar = sd.gpar;
	IntVector boxcrd;
	IntMatrix looplim;
	int id_j, r, s, t;
	real h_avg, dist_squared;

	par.InitNeighbours(gv.sph_maxnbr);

	// Loop over all particles associated to corresponding process
	for (int i=1; i<=gv.sph_np; i++)
	{
		if (!par[i].ToCompute()) continue;
		MParticle &par_i = par[i];
		par_i.mindist = 1.0e20;

		// Only do search is particle is inside the sort domain, otherwise leave nnbr = 0
		// continue to next particle
		// Identify which cell the particle lies in, and prevent loop over grid
		// cells which do not exist
		nv.GetBoxCoords(par_i.x, boxcrd);
		nv.GetLoopLimits(gv.sph_ndim, boxcrd, looplim);

		for (r=looplim(0,0); r<=looplim(1,0); r++)
		for (s=looplim(0,1); s<=looplim(1,1); s++)
		for (t=looplim(0,2); t<=looplim(1,2); t++)
		{
		// first particle in grid cell
		id_j = nv.GetGridValue(r, s, t);

		while (id_j != LINKEDLISTEND)
		{
			MCommonParticle &par_j = (id_j>-1) ? *(MCommonParticle *) &par[id_j] :
												 *(MCommonParticle *) &gpar[-id_j];
			
            // Do not count the particle itself
			if (id_j==i) { id_j = par_j.llpointer; continue; }

			dist_squared = par_i.DistanceSquared(par_j, gv.sph_ndim);
			h_avg = nv.m_factor * AVG(par_i.h, par_j.h);

			// Check if particle is within 2h domain around i and same material
			// material comparison commented out for surface tension, July 2008.
			if (dist_squared<SQR(h_avg) /*&& par_j.mat==par_i.mat*/)
			{
				par_i.mindist = MIN(par_i.mindist, sqrt(dist_squared));
				if (par_i.nnbr < gv.sph_maxnbr)
					par_i.AddNeighbour(id_j);
				else {
					PrintScreenLog("ERROR: Neighbour list exceeded maximum.");
					return ERROR;
				}
			}
			// Go further following llpointer
			id_j = par_j.llpointer;
		} // while loop end
		} // 3 for loops end
	} // upper for loop (i)

	return OK;
}

void MSolution::PrintScreenLog(const std::string &msg)
{
	printf("%s\n", msg.c_str());
	fprintf(m_simdata->m_filevars.f_logfile, "%s\n", msg.c_str());
}

int MSolution::FindNewMaxNeighbour()
{
	MSimulationData &sd = *m_simdata;
	MGlobalVars &gv = sd.m_globvars;
	MParticleData &par = sd.par;

	// Find new maxnbr
	gv.sph_maxnbr = 0;

 	for (int i=1; i<=gv.sph_np; i++)
 		gv.sph_maxnbr = MAX(gv.sph_maxnbr, par[i].nnbr);

	return OK;
}

int MSolution::Solution()
{
	MSimulationData &sd = *m_simdata;
	MGlobalVars &gv = sd.m_globvars;
	MOptionVars &ov = sd.m_optvars;
	MOutputVars &outv = sd.m_outvars;
	bool iterate = true, bPlotWritten=false;
	char msg[300];

	// Reserve memory for MPI transfer table
	HandleMpiMemory(true);

	while (iterate)
	{
		// Calculate new neighbour set
		if (Neighbours()!=OK) return ERROR;

		TransferParticles(0);
		TransferParticles(1);
		ExchangeData(0);
		ExchangeData(1);

		// Setup ghost particles
 		if (sd.m_optvars.sph_boundary) GhostSetup();
 		// Linked list neighbours after exchanging data among processes
 		if (LinkedListNeighbours()!=OK) return ERROR;
 		// Rate of deformation
		if (Strain()!=OK) return ERROR;
		// Update smoothing length if necessary, to be implemented
		if (ov.sph_h_opt>0) UpdateH();
		// Calculate normals needed for surface tension effects
		if (sd.m_optvars.sph_surface_tension /*||sd.m_optvars.sph_rigid_body*/) CalculateSurfaceNormals();
		// Material modeling
		if (Constitutive()!=OK) return ERROR;
		// Update time and set new time step
		if (CalcCriticalTimestep()!=OK) return ERROR;
		// Time advance
		gv.sph_dtold = gv.sph_dt;
		gv.sph_dt = MIN(gv.sph_tssf * gv.sph_critts, 1.1*gv.sph_dtold);
		gv.sph_ptime += gv.sph_dt;
		// Solve momentum equation and calculate new acceleration
		if (Momentum()!=OK) return ERROR;
		// Impose accel. boundary conditions
		if (BoundCond()!=OK) return ERROR;
		// Calculate total energy, to be implemented
		//CalcTotalE();
		// Write to state plot files if necessary
		if ( gv.sph_ptime>outv.sph_nextsttime || gv.sph_ptime>gv.sph_endtime )
		{
			TransferParticlesToRoot(false);
			if (sd.my_id==ROOT_PROCESS)
			{
				m_output->StateOutput();
				sprintf(msg, "\tState plot written at time %10.4e", gv.sph_ptime);
				PrintScreenLog(msg);
				PlotDroppingCylinder(sd.m_filevars.sph_filein+".drp");
			}
			outv.sph_nextsttime += outv.sph_stpltime;
		}
		// Time history state output, to be implemented
		// Write problem status to log file
		if (sd.my_id==ROOT_PROCESS && gv.sph_timestep%gv.sph_status_interval==0 )
		{
			sprintf(msg, "\tProblem status for time-step %6d\t Time: %10.4e\t dt: %10.4e",
					gv.sph_timestep, gv.sph_ptime, gv.sph_dt);
			PrintScreenLog(msg);
		}
		// Update particle velocities
		if (UpdateVelocity()!=OK) return ERROR;
		// Update particle positions
		if (MoveParticles()!=OK) return ERROR;
		// Write plot output in right time
		if (gv.sph_ptime>=outv.sph_plottime && !bPlotWritten) {
			TransferParticlesToRoot(false);
			if ( sd.my_id==ROOT_PROCESS) PlotOutput("plot_out.dat");
			bPlotWritten=true;
		}
		// Write CPU time and memory data, to be implemented
		// Increment time step counter
		gv.sph_timestep++;
		// Check for run termination conditions
 		if (gv.sph_ptime-gv.sph_dt >= gv.sph_endtime) iterate=false;
		// Exchange data at the end of the time step
	}

	return OK;
}

// Control function calling the right one
int MSolution::Strain()
{
	MSimulationData &sd = *m_simdata;
	MGlobalVars &gv = sd.m_globvars;
	MParticleData &par = sd.par;

	// Half step positions back in time to that positions and velocities are held at the same time
	for (int i=1; i<=gv.sph_np; i++)
	{
		if (!par[i].ToCompute() || par[i].lennardJones) continue;
		par[i].x -= 0.5 * par[i].v * gv.sph_dt;
	}

	switch (gv.sph_disctype)
	{
		// Basic SPH, Morris 1997. equal strain calculation
		case 0: case 2:
			ClassicStrain();
			break;
		// Mixed correction, to be implemented
		case 1:
			//CorrectedStrain();
			break;
	}

	// Half step return
	for (int i=1; i<=gv.sph_np; i++)
	{
		if (!par[i].ToCompute() || par[i].lennardJones) continue;
		par[i].x += 0.5 * par[i].v * gv.sph_dt;
	}

	return OK;
}

int MSolution::ClassicStrain()
{
	MSimulationData &sd = *m_simdata;
	MGlobalVars &gv = sd.m_globvars;
	MParticleData &par = sd.par;
	MGhostParticleData &gpar = sd.gpar;
	int j, id_j;
	real Volj, havg, dWdr;
	MKernel *pkernel;
	RealVector dWdx;
	RealMatrix grad_v, trans_grad_v;

	//	Uses the 'standard' sph approach
	for (int i=1; i<=gv.sph_np; i++)
	{
		if (!par[i].ToCompute() || par[i].lennardJones /*|| InRigidBody(par[i])*/) continue;
		MParticle &par_i = par[i];
		par_i.rod.clear();
		par_i.spin.clear();
		grad_v.clear();

		for (j=0; j<par_i.nnbr; j++)
		{
			id_j = par_i.m_nbrlist(j);
			MCommonParticle &par_j = (id_j>-1) ? *(MCommonParticle *)&par[id_j] :
												 *(MCommonParticle *)&gpar[-id_j];
			// do not consider Lennard-Jones neighbours when computing strain
			if (par_j.lennardJones) continue;
			Volj = par_j.mass / par_j.rho;
			havg = AVG(par_i.h, par_j.h);
			//pkernel = GenerateKernel(havg);
			//pkernel->GradW(dWdx, dWdr, par_i.x, par_j.x);
			//delete pkernel;
			MKernelBSpline::GradW(dWdx, dWdr, par_i.x, par_j.x, gv.sph_ndim, havg);
			// Calculate velocity gradient
			for (int r=0; r<gv.sph_ndim; r++) for (int s=0; s<gv.sph_ndim; s++)
				grad_v(r,s) += Volj * (par_j.v(r)-par_i.v(r))*dWdx(s);
		}

		for (int i=0; i<3; i++)
			for (int j=0; j<3; j++)
				trans_grad_v(i,j) = grad_v(j,i);

		// Rate of deformation and spin from
		par_i.rod  = 0.5 * (grad_v + trans_grad_v);
		par_i.spin = 0.5 * (grad_v - trans_grad_v);
		// ROD trace
		par_i.CalcTraceROD(gv.sph_ndim);
	} // for loop end

	return OK;
}

// Function used to generate appropriate kernel object according to input
// The kernel object must be deleted afterwards
MKernel* MSolution::GenerateKernel(real havg)
{
	MSimulationData &sd = *m_simdata;
	MGlobalVars &gv = sd.m_globvars;
	MKernel *pkernel=NULL;

	switch (m_simdata->m_optvars.sph_krtype)
	{
		// Classic B-Spline kernel
		case 1:
			pkernel = new MKernelBSpline(gv.sph_ndim, havg);
			break;
		default:
			PrintScreenLog("ERROR: Kernel option not supported.\n");
	}

	return pkernel;
}

int MSolution::Constitutive()
{
	MSimulationData &sd = *m_simdata;
	MGlobalVars &gv = sd.m_globvars;
	MParticleData &par = sd.par;
	MMaterialData &mat = sd.mat;

	// Update rho using Trace(ROD)
	RhoUpdate();

//	for (int i=1; i<=gv.sph_np; i++)
//	{
//		if (!par[i].ToCompute() || par[i].lennardJones /*|| InRigidBody(par[i])*/) continue;
//		MMaterial &mater = mat[ par[i].mat ];
//		MConstitutive::Constitutive(par[i], mater, &sd.m_filevars.f_logfile);
//	}

    for (int i=1; i<=gv.sph_np; i++)
	{
    	if (!par[i].ToCompute() || par[i].lennardJones /*|| InRigidBody(par[i])*/) continue;
		MMaterial &mater = mat[ par[i].mat ];
		MConstitutive::Constitutive(par[i], mater, &sd.m_filevars.f_logfile, gv.sph_dt);
	}

	return OK;
}

int MSolution::RhoUpdate()
{
	MSimulationData &sd = *m_simdata;
	MGlobalVars &gv = sd.m_globvars;
	MParticleData &par = sd.par;

	for (int i=1; i<=gv.sph_np; i++)
	{
		if (!par[i].ToCompute() || par[i].lennardJones /*|| InRigidBody(par[i])*/ ) continue;
		par[i].rhoold = par[i].rho;
		par[i].rho *= ( 1.0 - par[i].tracerod * gv.sph_dt );
	}

	return OK;
}

int MSolution::Momentum()
{
	MSimulationData &sd = *m_simdata;
	//Update Ghost particle stress, to be implemented

	switch (sd.m_globvars.sph_disctype)
	{
		// Basic SPH
		case 0:
			ClassicAcceleration();
			if (sd.m_optvars.sph_lennardjones) LennardJonesAcceleration();
			break;
		// Mixed correction
		case 1:
			//CorrectedAcceleration();
			break;
		// According to Morris et al. 1997.
		case 2:
			MorrisAcceleration();
			if (sd.m_optvars.sph_lennardjones) LennardJonesAcceleration();
			if (sd.m_optvars.sph_rigid_body) { RigidBodyAcceleration(); /*RigidBodyForce();*/ }
			break;
	}

	return OK;
}

int MSolution::ClassicAcceleration()
{
	MSimulationData &sd = *m_simdata;
	MGlobalVars &gv = sd.m_globvars;
	MParticleData &par = sd.par;
	MGhostParticleData &gpar = sd.gpar;
	int i, j, id, m, n;
	real Volj, havg, dWdr;
	RealVector dWdx(3), deltasig(3);

	// svp=start velocity point, evp=end velocity point
	for (i=gv.sph_svp; i<=gv.sph_evp; i++)
	{
		if (!par[i].ToCompute() || par[i].lennardJones) continue;
		MParticle &pari = par[i];
		pari.a.clear();

		for (j=0; j<pari.nnbr; j++)
		{
			id = pari.m_nbrlist(j);
			MCommonParticle &parj = (id>-1) ? *(MCommonParticle *)&par[id] :
											  *(MCommonParticle *)&gpar[-id];

			// do not consider Lennard-Jones neighbours when computing acceleration
			if (parj.lennardJones) continue;
			Volj = parj.mass/parj.rho;
			havg = AVG(pari.h, parj.h);
			//MKernel *pkernel = GenerateKernel(havg);
			//pkernel->GradW(dWdx, dWdr, pari.x, parj.x);
			//delete pkernel;
			MKernelBSpline::GradW(dWdx, dWdr, pari.x, parj.x, gv.sph_ndim, havg);

			deltasig.clear();
			for (n=0; n<gv.sph_ndim; n++) for (m=0; m<gv.sph_ndim; m++)
				deltasig(n) += ( (parj.sigma(m,n)-parj.q(m,n)) / SQR(parj.rho) +
		                       (pari.sigma(m,n)-pari.q(m,n)) / SQR(pari.rho) ) * dWdx(m);

			pari.a += parj.mass * deltasig;
		}
	}

	return OK;
}

int MSolution::MorrisAcceleration()
{
	MSimulationData &sd = *m_simdata;
	MGlobalVars &gv = sd.m_globvars;
	MParticleData &par = sd.par;
	MMaterialData &mat = sd.mat;
	MGhostParticleData &gpar = sd.gpar;
	int j, id_j, r;
	real Volj, havg, dWdr, rab2, product, viscosity_i, viscosity_j;
	RealVector dWdx(3);

	for (int i=1; i<=gv.sph_np; i++)
	{

		par[i].a.clear();
		if (!par[i].ToCompute() || par[i].lennardJones) continue;
		//if (InRigidBody(par[i])) continue;
		MParticle &pari = par[i];
		viscosity_i = mat[pari.mat].av_l;

		for (j=0; j<pari.nnbr; j++)
		{
			id_j = pari.m_nbrlist(j);
			MCommonParticle &parj = (id_j>-1) ? *(MCommonParticle *)&par[id_j] :
											    *(MCommonParticle *)&gpar[-id_j];

			// do not consider Lennard-Jones neighbours when computing acceleration
			if (parj.lennardJones) continue;
			viscosity_j = mat[parj.mat].av_l;
			Volj = parj.mass/parj.rho;
			havg = AVG(pari.h, parj.h);
			MKernelBSpline::GradW(dWdx, dWdr, pari.x, parj.x, gv.sph_ndim, havg);

			for (r=0, rab2=0.0, product=0.0; r<gv.sph_ndim; r++) {
		        rab2 += SQR(pari.x(r) - parj.x(r));
		        product += (pari.x(r) - parj.x(r)) * dWdx(r);
			}

			// First term expresses local pressure change, second expresses viscosity forces
			pari.a += -parj.mass * ( parj.p/SQR(parj.rho) + pari.p/SQR(pari.rho) ) * dWdx;
			if (InRigidBody(parj) && InRigidBody(pari)) continue;
			pari.a +=  parj.mass * (viscosity_i*pari.rho + viscosity_j*parj.rho) * product /
						   (pari.rho*parj.rho) / (rab2+0.01*SQR(havg)) * (pari.v-parj.v);
		} // for j
		// Additional surface tension member, coefficient of material applied to particle i
		if (sd.m_optvars.sph_surface_tension)
			pari.a += -mat[pari.mat].surfaceTensionCoeff/pari.rho * pari.normalizedNormalDiv*pari.normal;
	} // for i

	return OK;
}

int MSolution::UpdateVelocity()
{
	MSimulationData &sd = *m_simdata;
	MGlobalVars &gv = sd.m_globvars;
	MOptionVars &ov = sd.m_optvars;
	MParticleData &par = sd.par;
	MGhostParticleData &gpar = sd.gpar;
	int j, id_j;
	real dtn, epsilon, W, Volj;

	// Calculate delta t at time n
	dtn = AVG(gv.sph_dt, gv.sph_dtold);

	// Update velocities
	for (int i=1; i<=gv.sph_np; i++)
	{
		if (!par[i].ToCompute() || par[i].lennardJones || InRigidBody(par[i])) continue;

		par[i].v += dtn * par[i].a;
		par[i].CalculateVabs(gv.sph_ndim);
		//////////////////// Impose cylinder velocity bound conditions/////////////////
		//if (ov.bConstrainVelocityCylinder) ConstrainVelocityCylinder(par[i]);
		//////////////////// Impose cylinder velocity bound conditions/////////////////

		// Temporary for Couette flow - Decuzzi example
		//if (par[i].x(1)>6.01e-06)	par[i].v = 6.13e-03, 0, 0;
		//if (par[i].x(0)>6.01e-06)	par[i].v = 0, 6.13e-03, 0;

		// Temporary for Morris example, prescribed velocities at the entrance, March 2010.

		//if ( par[i].x(1)>1.e-03 && par[i].x(1)<1.5e-03 ) {
		//	real R=0.5e-03, v0=1.25e-05;
		//	real vy = (1 - par[i].x(0)*par[i].x(0) / (R*R))*v0;
		//	par[i].v.clear();
		//	par[i].v(1)=vy;
		//}

		// Temporary for shear cavity example, 
		//prescribed velocities at the top, May 2012.

		//if ( par[i].x(2)>=0.96e-03 && par[i].x(2)<1.e-03) {
		//	
                  //      par[i].v.clear();
                    //    par[i].v(0) = 1.e-04;
                      //  par[i].v(1) = 1.e-04;

		//}
        /*     if ( par[i].x(1)>=1.e-03) {
			
                        par[i].v.clear();
                        par[i].v(0) = 1.e-03;
                     //   par[i].v(1) = 1.e-04;

		}
		*/


	}

	// Rigid body velocity update
	if (ov.sph_rigid_body) RigidBodyVelocity();

	// XSPH option. Correct velocity. Note that currently v and x are held at different times
	if (ov.sph_veloc_opt != 1) return OK;

	epsilon = 0.05;

	for (int i=1; i<=gv.sph_np; i++)
	{
		if (!par[i].ToCompute() || par[i].lennardJones || InRigidBody(par[i]) ) continue;
		MParticle &pari = par[i];
		pari.smooth_v.clear();

		for (j=0; j<pari.nnbr; j++)
		{
			id_j = pari.m_nbrlist(j);
			MCommonParticle &parj = (id_j>-1) ? *(MCommonParticle *)&par[id_j] :
											  	*(MCommonParticle *)&gpar[-id_j];

			MKernel *kernel = GenerateKernel(parj.h);
			W = kernel->Kernel(pari.x, parj.x);
			delete kernel;

			// do not consider Lennard-Jones and rigid body particles in XSPH
			if (parj.lennardJones || InRigidBody(parj)) continue;
			Volj = parj.mass / AVG(parj.rho,pari.rho);
			pari.smooth_v += Volj * (parj.v-pari.v) * W;
		}

		// Finally, update the velocity
		pari.v += epsilon*pari.smooth_v;
		// Impose displacement boundary condition
		ConstrainQuantity(pari.v, pari.dispbc);
	}

	return OK;
}

int MSolution::BoundCond()
{
	MSimulationData &sd = *m_simdata;
	MGlobalVars &gv = sd.m_globvars;
	MBaseAccelVars &bav = sd.m_baseaccelvars;
	MParticleData &par = sd.par;

	// Apply base accelerations
	if (bav.sph_baseaccel)
	{
		for (int i=1; i<=gv.sph_np; i++)
		{
			if (!par[i].ToCompute() || par[i].lennardJones) continue;
			// Privremeno dodato da se kuglica ne bi ubrzavala (ili da se samo ona ubrzava)
			if (!InRigidBody(par[i])) par[i].a += bav.sph_base_a;
		}
	}

	// Apply nodal displacement boundary conditions
	for (int i=1; i<=gv.sph_np; i++)
	{
		if (!par[i].ToCompute() || par[i].lennardJones) continue;
		ConstrainQuantity(par[i].a, par[i].dispbc);
	}

	return OK;
}

// Impose displacement boundary condition on quantity (v or a)
int MSolution::ConstrainQuantity(RealVector &quantity, int code)
{
 	switch (code)
 	{
		case 1:
			// constrained x displacement
			quantity(0) = 0.0;
			break;
		case 2:
			// constrained y displacement
			quantity(1) = 0.0;
			break;
		case 3:
			// constrained z displacement
			quantity(2) = 0.0;
			break;
		case 4:
			// constrained x and y displacement
			quantity(0) = 0.0;
			quantity(1) = 0.0;
			break;
		case 5:
			// constrained y and z displacement
			quantity(1) = 0.0;
			quantity(2) = 0.0;
			break;
		case 6:
			// constrained z and x displacement
			quantity(2) = 0.0;
			quantity(0) = 0.0;
			break;
		case 7:
			// constrained x, y and z displacement
			quantity.clear();
			break;
 	}

 	return OK;
}

int MSolution::MoveParticles()
{
	MSimulationData &sd = *m_simdata;
	MGlobalVars &gv = sd.m_globvars;
	MParticleData &par = sd.par;
	MSymPerVars &spv = sd.m_sympervars;
	real overhead;
	int i;

	for (i=0; i<3; i++) {
		gv.sph_coord_minmax(0,i) =  1.0e20;
		gv.sph_coord_minmax(1,i) = -1.0e20;
	}

	// Update positions and store max and min values for linked list
	for (i=1; i<=gv.sph_np; i++)
		if (par[i].ToCompute() && !par[i].lennardJones)
			par[i].x += par[i].v * gv.sph_dt;

	// Treat periodic boundary and create new boundary box
	for (i=1; i<=gv.sph_np; i++)
	{
		for (int m=0; m<gv.sph_ndim; m++)
		{
			// Check periodic boundary conditions
			if ( spv.sph_boundary_code(1,m)==2 && par[i].x(m)>=spv.sph_boundary_x(1,m) ) {
				overhead = par[i].x(m)-spv.sph_boundary_x(1,m);
				par[i].x(m) = spv.sph_boundary_x(0,m) + overhead;
			}
			if ( spv.sph_boundary_code(0,m)==2 && par[i].x(m)<=spv.sph_boundary_x(0,m) ) {
				overhead = spv.sph_boundary_x(0,m) - par[i].x(m);
				par[i].x(m) = spv.sph_boundary_x(0,m) - overhead;
			}
			// Create new boundbox
			gv.sph_coord_minmax(0,m) = MIN(gv.sph_coord_minmax(0,m), par[i].x(m));
			gv.sph_coord_minmax(1,m) = MAX(gv.sph_coord_minmax(1,m), par[i].x(m));
		}
	}

	return OK;
}

int MSolution::GhostSetup()
{
	MSimulationData &sd = *m_simdata;
	MNeighbourVars &nv = sd.m_neighbourvars;
	MSymPerVars &spv = sd.m_sympervars;
	IntVector up_limit(3), down_limit(3);
	RealMatrix &minadd=spv.sph_mincoord_add, &maxadd=spv.sph_maxcoord_add;
	RealVector xmin(3),xmax(3),ymin(3),ymax(3),zmin(3),zmax(3);

	for (int i=0; i<3; i++) {
		xmin(i) = minadd(i,0);
		xmax(i) = maxadd(i,0);
		ymin(i) = minadd(i,1);
		ymax(i) = maxadd(i,1);
		zmin(i) = minadd(i,2);
		zmax(i) = maxadd(i,2);
	}

	sd.gpar.DeleteContents();

	// Xmin
	if (spv.sph_boundary_code(0,0)==2)
	{
		down_limit(0)=1; down_limit(1)=0; down_limit(2)=0;
		up_limit(0)=1; up_limit(1)=nv.sph_gridlim(1)-1; up_limit(2)=nv.sph_gridlim(2)-1;
		if ( GenerateGhostParticles(down_limit,up_limit,xmin) != OK ) return ERROR;
	}
	// Xmax
	if (spv.sph_boundary_code(0,0)==2)
	{
		down_limit(0)=nv.sph_gridlim(0)-2; down_limit(1)=0; down_limit(2)=0;
		up_limit(0)=nv.sph_gridlim(0)-2; up_limit(1)=nv.sph_gridlim(1)-1; up_limit(2)=nv.sph_gridlim(2)-1;
		if ( GenerateGhostParticles(down_limit,up_limit,xmax) != OK ) return ERROR;
	}
	// Ymin
	if (spv.sph_boundary_code(0,1)==2)
	{
		down_limit(0)=0; down_limit(1)=1; down_limit(2)=0;
		up_limit(0)=nv.sph_gridlim(0)-1; up_limit(1)=1; up_limit(2)=nv.sph_gridlim(2)-1;
		if ( GenerateGhostParticles(down_limit,up_limit,ymin) != OK ) return ERROR;
	}
	// Ymax
	if (spv.sph_boundary_code(1,1)==2)
	{
		down_limit(0)=0; down_limit(1)=nv.sph_gridlim(1)-2; down_limit(2)=0;
		up_limit(0)=nv.sph_gridlim(0)-1; up_limit(1)=nv.sph_gridlim(1)-2; up_limit(2)=nv.sph_gridlim(2)-1;
		if ( GenerateGhostParticles(down_limit,up_limit,ymax) != OK ) return ERROR;
	}
	// Zmin
	if (spv.sph_boundary_code(0,2)==2)
	{
		down_limit(0)=0; down_limit(1)=0; down_limit(2)=1;
		up_limit(0)=nv.sph_gridlim(0)-1; up_limit(1)=nv.sph_gridlim(1)-1; up_limit(2)=1;
		if ( GenerateGhostParticles(down_limit,up_limit,zmin) != OK ) return ERROR;
	}
	// Zmax
	if (spv.sph_boundary_code(1,2)==2)
	{
		down_limit(0)=0; down_limit(1)=0; down_limit(2)=nv.sph_gridlim(2)-2;
		up_limit(0)=nv.sph_gridlim(0)-1; up_limit(1)=nv.sph_gridlim(1)-1; up_limit(2)=nv.sph_gridlim(2)-2;
		if ( GenerateGhostParticles(down_limit,up_limit,zmax) != OK ) return ERROR;
	}

	return OK;
}

// Generate ghost particles from particles in cells sent in limit vars
int MSolution::GenerateGhostParticles(IntVector &down_limit, IntVector &up_limit, RealVector &addVec)
{
	MSimulationData &sd = *m_simdata;
	MNeighbourVars &nv = sd.m_neighbourvars;
	MParticleData &par = sd.par;
	MGhostParticleData &gpar = sd.gpar;
	int i, j, k, id, n, min_lim, max_lim;
	IntVector offset(3);
	offset.clear();

	// check if periodic boundary is on X (domain decomposition made on X)
	if (sd.m_sympervars.sph_boundary_code(0,0) != 2) {
		min_lim = MAX( 0, box_procmin(sd.my_id) );
		max_lim = MIN( nv.sph_gridlim(0)-1, box_procmax(sd.my_id)+1 );
	}
	else {
		min_lim = down_limit(0);
		max_lim = up_limit(0);
	}

	// Move cell to the other end of the model
	for (n=0; n<3; n++)
		if (addVec(n)!=0.0) offset(n) = SGN(addVec(n)) * (nv.sph_gridlim(n)-2);

	for (k=down_limit(2); k<=up_limit(2); k++)
	for (j=down_limit(1); j<=up_limit(1); j++)
	for (i=min_lim; i<=max_lim; i++)
	{
		id = nv.GetGridValue(i,j,k);
		if (id<0) continue; // if cell does not contain real particles, continue

		while (id>0)
		{
			gpar.CreateGhostFromReal(par[id], id, addVec,
				                     nv.GetGridValue( i+offset(0),j+offset(1),k+offset(2) ) );
			id = par[id].llpointer;
		}
	}

	return OK;
}

// Function to export quantity graph temporarily, hardwired, called by Solution()
int MSolution::PlotOutput(const std::string &filename)
{
	MSimulationData &sd = *m_simdata;
	MGlobalVars &gv = sd.m_globvars;
	MOutputVars &outv = sd.m_outvars;
	MParticleData &par = sd.par;
	int i, j, num_parts=101;
	real W, Volj, Volk, quantity, quantity1;
	RealVector xi(3);
	RealVector &start_point = outv.sph_plot_start_x;
	RealVector &end_point = outv.sph_plot_end_x;
	FILE *file;

	// Flow between cylinders example - start_point_x=0.5e-03_d, end_point_x=2.5e-03_d

	if ( (file = fopen(filename.c_str(), "wt")) == NULL )
    {
    	printf("ERROR: Cannot open %s for plot output.\n", filename.c_str() );
    	return ERROR;
    }

	// The velocity distribution over the line
	for (i=0; i<=num_parts; i++)
	{
	    xi = start_point + i * (end_point-start_point)/num_parts;

	    for (j=gv.sph_ssp, quantity=0.0; j<=gv.sph_evp; j++)
	    {
			MParticle &parj = par[j];

			MKernel *pkernel = GenerateKernel(parj.h);
			W = pkernel->Kernel(xi, parj.x);
			delete pkernel;

			Volj = parj.mass / parj.rho;
	        quantity += Volj * parj.rho * W;
	    }

	    fprintf(file, "%13.5e%13.5e\n", quantity, xi(0) );
	}

	// Contour plot for "flow past cylinders" example
	/*
	for (i=0; i<num_parts; i++)
		for (int j=0; j<num_parts; j++)
		{
		    RealVector offset(i*1.0e-01/num_parts, j*1.0e-01/num_parts);
			xi = offset;

		    quantity = 0.;
		    quantity1 = 0.;
			for (int k=gv.sph_ssp; k<=gv.sph_evp; k++)
		    {
				MParticle &park = par[k];

				MKernel *pkernel = GenerateKernel(park.h);
				W = pkernel->Kernel(xi, park.x);
				delete pkernel;

				Volk = park.mass / park.rho;
		        quantity += Volk * sqrt( SQR(park.v(0)) + SQR(park.v(1)) ) * W;

				pkernel = GenerateKernel(3*park.h);
				W = pkernel->Kernel(xi, park.x);
				delete pkernel;

		        quantity1 += Volk * park.p * W;
		    }
		    fprintf(file, "%13.5e%13.5e%13.5e%13.5e\n", xi(0), xi(1), quantity, quantity1);
		}

	// Putanja 1
	for (i=0; i<num_parts; i++)
	{
	    RealVector offset(i*1.0e-01/num_parts, 1.0e-01/2);
		xi = offset;

	    quantity = 0.;
		for (int k=gv.sph_ssp; k<=gv.sph_evp; k++)
	    {
			MParticle &park = par[k];

			MKernel *pkernel = GenerateKernel(park.h);
			W = pkernel->Kernel(xi, park.x);
			delete pkernel;

			Volk = park.mass / park.rho;
	        quantity += Volk * sqrt( SQR(park.v(0)) + SQR(park.v(1)) ) * W;
	    }
	    fprintf(file, "%13.5e%13.5e\n", xi(0), quantity);
	}

	// Putanja 2
	for (i=0; i<num_parts; i++)
	{
	    RealVector offset(i*1.0e-01/num_parts, 0.95e-01);
		xi = offset;

	    quantity = 0.;
		for (int k=gv.sph_ssp; k<=gv.sph_evp; k++)
	    {
			MParticle &park = par[k];

			MKernel *pkernel = GenerateKernel(park.h);
			W = pkernel->Kernel(xi, park.x);
			delete pkernel;

			Volk = park.mass / park.rho;
	        quantity += Volk * sqrt( SQR(park.v(0)) + SQR(park.v(1)) ) * W;
	    }
	    fprintf(file, "%13.5e%13.5e\n", xi(0), quantity);
	}
	*/

	PrintScreenLog("Plot written into file.\n");

	fclose(file);
	return OK;
}

// Hardwired for flow between cylinders, called by UpdateVelocity() for each particle
int MSolution::ConstrainVelocityCylinder(MParticle &p)
{
	real rho, phi, R1=0.5e-03, R2=2.5e-03, tolerance=1.e-08;
	real W1 = m_simdata->m_optvars.InnerW;
	real W2 = m_simdata->m_optvars.OuterW;

	rho = sqrt( SQR(p.x(0)) + SQR(p.x(1)) );
	phi = atan2( p.x(1), p.x(0) );
	real w_interpolated = W1+(W2-W1)/(R2-R1) * (rho-R1);

    if ( rho < R1+tolerance )
    {
    	p.v(0) =  w_interpolated * rho * sin(phi);
    	p.v(1) = -w_interpolated * rho * cos(phi);
    }
    else if ( rho > R2-tolerance )
    {
    	p.v(0) =  w_interpolated * rho * sin(phi);
    	p.v(1) = -w_interpolated * rho * cos(phi);
	}

	return OK;
}

int MSolution::InitializeMpiProcess()
{
	MSimulationData &sd = *m_simdata;
	MNeighbourVars &nv = sd.m_neighbourvars;
	int nmin, nleft, nboxes, i, first, last;

	box_procmin.resize(sd.numprocs);
	box_procmax.resize(sd.numprocs);

	nmin = nv.sph_gridlim(0) / sd.numprocs;
	nleft = nv.sph_gridlim(0) % sd.numprocs;
	for (i = 0, last=-1; i < sd.numprocs; i++)
	{
		nboxes = (i < nleft) ? nmin+1 : nmin;
		first = last+1;
		last = first+nboxes-1;
		box_procmin(i) = first;
		box_procmax(i) = last;
	}

	return OK;
}


int MSolution::ExchangeData(int phase)
{
	MSimulationData &sd = *m_simdata;
	MNeighbourVars &nv = sd.m_neighbourvars;
	MParticleData &par = sd.par;
	int id,l,n,count,r,boxcrd[3],send_line,my_id=sd.my_id;
	int proc_to_send, proc_to_receive;

	int left_proc = my_id!=0 ? my_id-1 : NOPROC;
	int right_proc = my_id!=sd.numprocs-1 ? my_id+1 : NOPROC;

	// In phase 0, the data is sent to the right process and received from the left process
	// In phase 1, the data is sent to the left process and received from the right process
	if (phase==0) {
		proc_to_send = right_proc;
		send_line = box_procmax(my_id);
		proc_to_receive = left_proc;
	} else {
		proc_to_send = left_proc;
		send_line = box_procmin(my_id);
		proc_to_receive = right_proc;
	}

	// If there is a process to be sent data to
	if (proc_to_send != NOPROC)
	{
		// Prepare data to be sent
		for (l=0,count=0; l<nv.sph_gridlim(1); l++)
		for (n=0; n<nv.sph_gridlim(2); n++)
		{
		if ( (id=nv.GetGridValue(send_line,l,n)) < 0 ) continue;
		while ( id != LINKEDLISTEND)
		{
			MParticle &part = par[id];
			mpiExchangeTable[count++] = id;
			mpiExchangeTable[count++] = send_line;
			mpiExchangeTable[count++] = l;
			mpiExchangeTable[count++] = n;
			PackMpiQuantities(true, count, part);
			id = part.llpointer;
		}
		}
		mpiExchangeTable[count++] = LINKEDLISTEND;
		TransferMpiBuffer(true, count, proc_to_send, my_id);
	}

	// If we expect data to be received from any process
	if (proc_to_receive != NOPROC)
	{
		TransferMpiBuffer(false, count, proc_to_receive, proc_to_receive);
		// Unpack received data
		count=0;
		while ( (id=(int)mpiExchangeTable[count++]) != LINKEDLISTEND )
		{
			MParticle &part = par[id];
			for (r=0; r<3; r++) boxcrd[r]=(int) mpiExchangeTable[count++];
			PackMpiQuantities(false, count, part);
			part.active = true;
			part.fromOtherProc = true;
			part.llpointer = nv.GetGridValue( boxcrd[0], boxcrd[1], boxcrd[2] );
			nv.SetGridValue( boxcrd[0], boxcrd[1], boxcrd[2], id );
		}
	}

	return OK;
}


int MSolution::TransferParticles(int phase)
{
	MSimulationData &sd = *m_simdata;
	MNeighbourVars &nv = sd.m_neighbourvars;
	MParticleData &par = sd.par;
	int id,l,n,count,r,boxcrd[3],send_line,my_id=sd.my_id;

	int proc_to_send, proc_to_receive;

	int left_proc = my_id!=0 ? my_id-1 : NOPROC;
	int right_proc = my_id!=sd.numprocs-1 ? my_id+1 : NOPROC;

	// In phase 0, the data is sent to the right process and received from the left process
	// In phase 1, the data is sent to the left process and received from the right process
	if (phase==0) {
		proc_to_send = right_proc;
		send_line = box_procmax(my_id)+1;
		proc_to_receive = left_proc;
	} else {
		proc_to_send = left_proc;
		send_line = box_procmin(my_id)-1;
		proc_to_receive = right_proc;
	}

	// If there is a process to be sent data to
	if (proc_to_send != NOPROC)
	{
		// Prepare data to be sent
		for (l=0,count=0; l<nv.sph_gridlim(1); l++)
		for (n=0; n<nv.sph_gridlim(2); n++)
		{
		if ( (id=nv.GetGridValue(send_line,l,n)) < 0 ) continue;
		while ( id != LINKEDLISTEND)
		{
			MParticle &part = par[id];
			mpiExchangeTable[count++] = id;
			mpiExchangeTable[count++] = send_line;
			mpiExchangeTable[count++] = l;
			mpiExchangeTable[count++] = n;
			PackMpiQuantities(true, count, part);
			// Particle is inactive in this process
			part.active = false;
			id = part.llpointer;
		}
		}
		mpiExchangeTable[count++] = LINKEDLISTEND;
		TransferMpiBuffer(true, count, proc_to_send, my_id);
	}

	// If we expect data to be received from any process
	if (proc_to_receive != NOPROC)
	{
		TransferMpiBuffer(false, count, proc_to_receive, proc_to_receive);
		// Unpack received data
		count=0;
		while ( (id=(int)mpiExchangeTable[count++]) != LINKEDLISTEND )
		{
			MParticle &part = par[id];
			for (r=0; r<3; r++) boxcrd[r]=(int) mpiExchangeTable[count++];
			PackMpiQuantities(false, count, part);
			part.llpointer = nv.GetGridValue( boxcrd[0], boxcrd[1], boxcrd[2] );
			nv.SetGridValue( boxcrd[0], boxcrd[1], boxcrd[2], id );
			part.active = true;
			part.fromOtherProc = false;
		}
	}

	// Cleanup cells from the line which was just sent to the other process
	if (proc_to_send!=NOPROC)
	{
		for (l=0; l<nv.sph_gridlim(1); l++)
		for (n=0; n<nv.sph_gridlim(2); n++)
			nv.SetGridValue(send_line, l, n, LINKEDLISTEND);
	}

	return OK;
}

/**
 * This function is used for state plot and to transfer rigid body particles.
 * If rigidBodyParticles==true, only rigid body particles are transfered, otherwise all of them.
 * Rigid body particle locations have to be known to all processes because, since
 * center of mass could not be computed properly.
 * ROOT_PROCESS then computes rigid body mass center and sends data to all others.
 */
int MSolution::TransferParticlesToRoot(bool rigidBodyParticles)
{
	MSimulationData &sd = *m_simdata;
	MParticleData &par = sd.par;
	MNeighbourVars &nv = sd.m_neighbourvars;
	int count, r, l , n, m, id, boxcrd[3], my_id=sd.my_id;

	// if process is ROOT_PROCESS, collect data, otherwise send it
	if (my_id == ROOT_PROCESS)
	{
		for (int proc=1; proc<sd.numprocs; proc++)
		{
			TransferMpiBuffer(false, count, proc, proc);

			// Clean cells
			for (m=box_procmin(proc); m<=box_procmax(proc); m++)
			for (l=0; l<nv.sph_gridlim(1); l++)
			for (n=0; n<nv.sph_gridlim(2); n++)
				nv.SetGridValue(m, l, n, LINKEDLISTEND);

			count=0;
			while ( (id=(int)mpiExchangeTable[count++]) != LINKEDLISTEND )
			{
				MParticle &part = par[id];
				for (r=0; r<3; r++) boxcrd[r]=(int) mpiExchangeTable[count++];
				PackMpiQuantities(false, count, part);
				part.llpointer = nv.GetGridValue( boxcrd[0], boxcrd[1], boxcrd[2] );
				nv.SetGridValue( boxcrd[0], boxcrd[1], boxcrd[2], id );
				part.active = false;
				part.fromOtherProc = true;
			}
		}
	}
	else // process is not ROOT_PROCESS
	{
		for (m=box_procmin(my_id),count=0; m<=box_procmax(my_id); m++)
		for (l=0; l<nv.sph_gridlim(1); l++)
		for (n=0; n<nv.sph_gridlim(2); n++)
		{
		if ( (id=nv.GetGridValue(m,l,n)) < 0 ) continue;
		while (id != LINKEDLISTEND)
		{
			MParticle &part = par[id];
			if (!part.active) { id = part.llpointer; continue; }
			if (rigidBodyParticles && !InRigidBody(part) ) { id = part.llpointer; continue; }
			mpiExchangeTable[count++] = id;
			mpiExchangeTable[count++] = m;
			mpiExchangeTable[count++] = l;
			mpiExchangeTable[count++] = n;
			PackMpiQuantities(true, count, part);
			id = part.llpointer;
		}
		}
		mpiExchangeTable[count++] = LINKEDLISTEND;
		TransferMpiBuffer(true, count, ROOT_PROCESS, my_id);
	}

	return OK;
}

// Method used for packing physical particle parameters into mpiExchangeTable and vice-versa
int MSolution::PackMpiQuantities(bool packing, int &count, MParticle &part)
{
	int r,s;

	if (packing) {
		for (r=0; r<3; r++) mpiExchangeTable[count++]=part.x(r);
		for (r=0; r<3; r++) mpiExchangeTable[count++]=part.v(r);
		for (r=0; r<3; r++) mpiExchangeTable[count++]=part.normalizedNormal(r);
		for (r=0; r<3; r++) for (s=0; s<3; s++) mpiExchangeTable[count++]=part.sigma(r,s);
		for (r=0; r<3; r++) for (s=0; s<3; s++) mpiExchangeTable[count++]=part.q(r,s);
		mpiExchangeTable[count++] = part.c;
		mpiExchangeTable[count++] = part.p;
		mpiExchangeTable[count++] = part.rho;
		mpiExchangeTable[count++] = part.tracerod;
		mpiExchangeTable[count++] = part.mindist;
		mpiExchangeTable[count++] = part.N;
		mpiExchangeTable[count++] = part.normalizedNormalDiv;
	}
	else {
		for (r=0; r<3; r++) part.x(r)=mpiExchangeTable[count++];
		for (r=0; r<3; r++) part.v(r)=mpiExchangeTable[count++];
		for (r=0; r<3; r++) part.normalizedNormal(r)=mpiExchangeTable[count++];
		for (r=0; r<3; r++) for (s=0; s<3; s++) part.sigma(r,s)=mpiExchangeTable[count++];
		for (r=0; r<3; r++) for (s=0; s<3; s++) part.q(r,s)=mpiExchangeTable[count++];
		part.c =    	 mpiExchangeTable[count++];
		part.p =    	 mpiExchangeTable[count++];
		part.rho =  	 mpiExchangeTable[count++];
		part.tracerod =  mpiExchangeTable[count++];
		part.mindist =   mpiExchangeTable[count++];
		part.N 	= 	 mpiExchangeTable[count++];
		part.normalizedNormalDiv = mpiExchangeTable[count++];
	}

	return OK;
}

int MSolution::HandleMpiMemory(bool reserving)
{
	int nparticles = m_simdata->m_globvars.sph_np;

	if (reserving)
		mpiExchangeTable = new real[ nparticles*50 ];
	else
		delete[] mpiExchangeTable;

	return OK;
}

// Transfers mpiExchangeBuffer contents from one process to another
// count aprameter is used only in case of sending
int MSolution::TransferMpiBuffer(bool sending, int &count, int proc, int tag)
{
	int my_id = m_simdata->my_id, transferredDoubles=0, packetSize, TRANSFER_END=1000000;
	const int MAX_PACKET_SIZE = 256;
	real dummy;
	MPI_Status status;

	if (sending) {
		while (transferredDoubles<count)
		{
			packetSize = transferredDoubles+MAX_PACKET_SIZE<count ?
					MAX_PACKET_SIZE : count-transferredDoubles;
			MPI_Send( &mpiExchangeTable[transferredDoubles], packetSize, MPI_DOUBLE, proc,
					my_id, MPI_COMM_WORLD );
			transferredDoubles+= packetSize;
		}
		MPI_Send( &dummy, 1, MPI_DOUBLE, proc, TRANSFER_END, MPI_COMM_WORLD );
	}
	else {
        do {
        	// Probe incomming message, check SIZE and TAG
    		MPI_Probe(proc, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            MPI_Get_count(&status, MPI_DOUBLE, &packetSize);
            // If it is the TRANDFER_END tag, dummy receive and stop, otherwise receive doubles
            if ( status.MPI_TAG != TRANSFER_END ) {
            	MPI_Recv( &mpiExchangeTable[transferredDoubles], packetSize, MPI_DOUBLE, proc, proc,
            			MPI_COMM_WORLD, &status );
            	transferredDoubles+= packetSize;
            }
            else
            	MPI_Recv( &dummy, 1, MPI_DOUBLE, proc, TRANSFER_END, MPI_COMM_WORLD, &status );
        } while(status.MPI_TAG != TRANSFER_END);
	}

	return OK;
}

// Routine -adds- acceleration caused by Lennard-Jones particles
int MSolution::LennardJonesAcceleration()
{
	MSimulationData &sd = *m_simdata;
	MGlobalVars &gv = sd.m_globvars;
	MParticleData &par = sd.par;
	MGhostParticleData &gpar = sd.gpar;
	int i, j, id;
	RealVector ort(3);
	real r, accelerationIntensity, D = gv.lennardJonesD, r0 = gv.lennardJonesR0;
    
	for (i=1; i<=gv.sph_np; i++)
	{
		if (!par[i].ToCompute() || par[i].lennardJones) continue;
		MParticle &pari = par[i];

		for (j=0; j<pari.nnbr; j++)
		{
			id = pari.m_nbrlist(j);
			MCommonParticle &parj = (id>-1) ? *(MCommonParticle *)&par[id] :
											  *(MCommonParticle *)&gpar[-id];

			// Continue only for neighbouring particles which are of Lennard-Jones type
			if (!parj.lennardJones) continue;

			// Acceleration computation
			r = pari.Distance(parj, gv.sph_ndim);
			ort = (pari.x - parj.x) * (1/r);
			// purely repulsive force
			accelerationIntensity = (r<r0) ? D/r/pari.mass * ( pow(r0/r,12) - pow(r0/r,6) ) : 0;
			pari.a += accelerationIntensity * ort;
		}
	}

	return OK;
}

// According to Morris 2000. used to calculate surface tension
int MSolution::CalculateSurfaceNormals()
{
	MSimulationData &sd = *m_simdata;
	MGlobalVars &gv = sd.m_globvars;
	MParticleData &par = sd.par;
	MGhostParticleData &gpar = sd.gpar;
	int i, j, k, id_j;
	real Volj, havg, dWdr, intensity, product;
	RealVector dWdx, substractNormals;

	// Obtain normals & normalized normals in this loop
	for (i=1; i<=gv.sph_np; i++)
	{
		if (!par[i].ToCompute()) continue;
		MParticle &pari = par[i];
		pari.normal = 0;

		for (j=0; j<pari.nnbr; j++)
		{
			id_j = pari.m_nbrlist(j);
			MCommonParticle &parj = (id_j>-1) ? *(MCommonParticle *)&par[id_j] :
											    *(MCommonParticle *)&gpar[-id_j];
			Volj = parj.mass/parj.rho;
			havg = AVG(pari.h, parj.h);
			MKernelBSpline::GradW(dWdx, dWdr, pari.x, parj.x, gv.sph_ndim, havg);

			pari.normal += Volj*(parj.color-pari.color)*dWdx;
		} // for j

		// Calculate normalized normals w/o disturbance, only considered > epsilon=0.01h
		for (k=0, intensity=0.0; k<gv.sph_ndim; ++k)
			intensity+= SQR(pari.normal(k));
		intensity = sqrt(intensity);
		pari.N = intensity > 1.e-2*pari.h ? 1 : 0;
		if (pari.N==1)
			pari.normalizedNormal = (1./intensity)*pari.normal;
		else
			pari.normalizedNormal = 0;
	} // for i

	// Following loop calculates normalized normal divergence needed to obtain s. tension force
	for (i=1; i<=gv.sph_np; i++)
	{
		if (!par[i].ToCompute()) continue;
		MParticle &pari = par[i];
		pari.normalizedNormalDiv = 0;

		for (j=0; j<pari.nnbr; j++)
		{
			id_j = pari.m_nbrlist(j);
			MCommonParticle &parj = (id_j>-1) ? *(MCommonParticle *)&par[id_j] :
											    *(MCommonParticle *)&gpar[-id_j];
			Volj = parj.mass/parj.rho;
			havg = AVG(pari.h, parj.h);
			MKernelBSpline::GradW(dWdx, dWdr, pari.x, parj.x, gv.sph_ndim, havg);

			substractNormals = parj.normalizedNormal-pari.normalizedNormal;
			for (k=0, product=0; k<gv.sph_ndim; ++k) product+=substractNormals(k)*dWdx(k);
			pari.normalizedNormalDiv += MIN(pari.N,parj.N) * Volj * product;
		} // for j

	} // for i

	return OK;
}

// Computes rigid body linear and angle acceleration
int MSolution::RigidBodyAcceleration()
{
	MSimulationData &sd = *m_simdata;
	MGlobalVars &gv = sd.m_globvars;
	MParticleData &par = sd.par;
	MGhostParticleData &gpar = sd.gpar;
	MRigidBody &rigidBody = sd.m_rigidbodydata[1];
	RealVector &A = rigidBody.A, &ALPHA = rigidBody.ALPHA, totalA, totalALPHA;
	int k, a, z, id;
	RealVector rka, fk, Rlocal;
	real distance, y, x, B;

	// First compute center of mass rigidBody.R
	TransferParticlesToRoot(true);
	rigidBody.CalculateCenterOfMass(par);
	MPI_Bcast(&rigidBody.R[0], 3, MPI_DOUBLE, ROOT_PROCESS, MPI_COMM_WORLD);

	// 'k' is a rigid body particle, 'a' is a fluid particle
	for (k=1, A=0, ALPHA=0; k<=gv.sph_np; k++)
	{
		if (!par[k].ToCompute() || !InRigidBody(par[k]) ) continue;
		MParticle &par_k = par[k];
		/*
		 * This part computes Monaghan B force exerted on rigid body, not used now
		 * should be used together with RigidBodyForce() which calculates oposit
		 * force exerted on fluid
		 *
		for (a=0, fk=0; a<(int)par_k.m_nbrlist.size(); a++)
		{
			id = par_k.m_nbrlist[a];
			MCommonParticle &par_a = (id>-1) ? *(MCommonParticle *)&par[id]  :
											   *(MCommonParticle *)&gpar[-id];

			// Continue only for fluid neighbours
			if (InRigidBody(par_a)) continue;
			// Acceleration computation
			distance = par_k.Distance(par_a, gv.sph_ndim);
			rka = par_a.x - par_k.x;
			for (z=0,y=0; z<gv.sph_ndim; z++) y+=par_a.normalizedNormal(z)*rka(z);
			y = fabs(y);
			x = sqrt( SQR(distance) - SQR(y) );
			B = rigidBody.B(x,y,par_a.c,par_a.h);
			fk += par_a.mass/(par_k.mass+par_a.mass) * B * par_a.normalizedNormal;
		}*/

		fk = 0;
		Rlocal = par_k.x - rigidBody.R;
		// add acceleration and angle acceleration between k and a
		A +=  par_k.mass * ( fk + par_k.a);
		fk = par_k.a;
		ALPHA(0) += par_k.mass*(Rlocal(1)*fk(2) - Rlocal(2)*fk(1));
		ALPHA(1) += par_k.mass*(Rlocal(2)*fk(0) - Rlocal(0)*fk(2));
		ALPHA(2) += par_k.mass*(Rlocal(0)*fk(1) - Rlocal(1)*fk(0));
	}

	A *= 1./rigidBody.totalMass;
	ALPHA *= 1./rigidBody.inertiaMoment;

	MPI_Reduce(&A[0], &totalA[0], 3, MPI_DOUBLE, MPI_SUM, ROOT_PROCESS, MPI_COMM_WORLD);
	MPI_Bcast(&totalA[0], 3, MPI_DOUBLE, ROOT_PROCESS, MPI_COMM_WORLD);
	MPI_Reduce(&ALPHA[0], &totalALPHA[0], 3, MPI_DOUBLE, MPI_SUM, ROOT_PROCESS, MPI_COMM_WORLD);
	MPI_Bcast(&totalALPHA[0], 3, MPI_DOUBLE, ROOT_PROCESS, MPI_COMM_WORLD);

	A = totalA;
	ALPHA = totalALPHA;
	if (sd.m_baseaccelvars.sph_baseaccel) A += sd.m_baseaccelvars.sph_base_a;

	return OK;
}

// Method to calculate velocity of each particle belonging to rigid body
int MSolution::RigidBodyVelocity()
{
	MSimulationData &sd = *m_simdata;
	MGlobalVars &gv = sd.m_globvars;
	MParticleData &par = sd.par;
	MRigidBody &rigidBody = sd.m_rigidbodydata[1];
	RealVector &V=rigidBody.V, &OMEGA=rigidBody.OMEGA;
	RealVector Rlocal;

	// Calculate delta t at time n
	real dtn = AVG(gv.sph_dt, gv.sph_dtold);
	V += rigidBody.A*dtn;
	OMEGA += rigidBody.ALPHA*dtn;

	for (int k=1; k<=gv.sph_np; k++)
	{
		if (!par[k].ToCompute() || !InRigidBody(par[k])) continue;
		Rlocal = par[k].x - rigidBody.R;
		par[k].v(0) = V(0) + OMEGA(1)*Rlocal(2) - OMEGA(2)*Rlocal(1);
		par[k].v(1) = V(1) + OMEGA(2)*Rlocal(0) - OMEGA(0)*Rlocal(2);
		par[k].v(2) = V(2) + OMEGA(0)*Rlocal(1) - OMEGA(1)*Rlocal(0);
	}

	return OK;
}

/*
 *  Computes fluid forces originated from the rigid body, used only with Monghan B force approach
 */
int MSolution::RigidBodyForce()
{
	MSimulationData &sd = *m_simdata;
	MGlobalVars &gv = sd.m_globvars;
	MParticleData &par = sd.par;
	MGhostParticleData &gpar = sd.gpar;
	int a, k, z, id;
	RealVector rak;
	real distance, y, x, B;

	for (a=1; a<=gv.sph_np; a++)
	{
		if (!par[a].ToCompute() || par[a].lennardJones || InRigidBody(par[a])) continue;
		MParticle &par_a = par[a];

		for (k=0; k<par_a.nnbr; k++)
		{
			id = par_a.m_nbrlist(k);
			MCommonParticle &par_k = (id>-1) ? *(MCommonParticle *)&par[id]  :
											   *(MCommonParticle *)&gpar[-id];

			// Continue only for neighbouring particles which are of Lennard-Jones type
			if (!InRigidBody(par_k)) continue;
			// Acceleration computation
			distance = par_a.Distance(par_k, gv.sph_ndim);
			rak = par_a.x - par_k.x;
			for (z=0,y=0.0; z<gv.sph_ndim; z++) y+=par_k.normalizedNormal(z)*rak(z);
			y = fabs(y);
			x = sqrt( SQR(distance) - SQR(y) );
			B = sd.m_rigidbodydata[1].B(x,y,par_a.c,par_a.h);
			// add acceleration between i and j
			par_a.a += -par_k.mass/(par_a.mass+par_k.mass) * B * par_k.normalizedNormal;
		}
	}

	return OK;
}

/**
 * Method for exporting values for dropping cylinder, hardwired, called by Solution()
 * Call this function only after TransferParticlesToRoot(false)
 */
int MSolution::PlotDroppingCylinder(const std::string &filename)
{
	MSimulationData &sd = *m_simdata;
	MGlobalVars &gv = sd.m_globvars;
	MParticleData &par = sd.par;
	int j, count;
	RealVector x;
	FILE* file;
	//MRigidBody &rigidBody = sd.m_rigidbodydata[1];

	// Open file for appending
	//if ( (file = fopen(filename.c_str(), "at+")) == NULL )
	if ( (file = fopen(filename.c_str(), "wt")) == NULL )
    {
    	printf("ERROR: Cannot open %s for plot output.\n", filename.c_str() );
    	return ERROR;
    }

	//for (j=1, x=0, count=0; j<gv.sph_np; j++)
	//	if (par[j].mat==2) { x+=par[j].x; count++; }


    //fprintf(file, "%13.5e%13.5e%13.5e%13.5e\n", gv.sph_ptime, rigidBody.R(0), rigidBody.R(1),
    //		rigidBody.OMEGA(2));

    for (j=1; j<gv.sph_np; j++)
    	fprintf(file, "%13.5e%13.5e\n", par[j].x(0), par[j].rho);

	fclose(file);
	return OK;
}

int MSolution::UpdateH()
{
	MSimulationData &sd = *m_simdata;
	MGlobalVars &gv = sd.m_globvars;
	MOptionVars &ov = sd.m_optvars;
	MNeighbourVars &nv = sd.m_neighbourvars;
	MParticleData &par = sd.par;
	real power, hmax;

	switch (ov.sph_h_opt) {
	case 1:
		nv.sph_hmax = 0.0;
		for (int i=1; i<=par.GetSize(); i++)
			if (par[i].ToCompute() && par[i].p<par[i].pcut)
			{
				par[i].hold = par[i].h;
				par[i].h = par[i].hold * (1 + par[i].tracerod*gv.sph_dt/gv.sph_ndim);
				nv.sph_hmax = MAX(nv.sph_hmax, par[i].h);
			}
		MPI_Allreduce(&nv.sph_hmax, &hmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
		nv.sph_hmax = hmax;
		break;
	case 2:
		nv.sph_hmax = 0.0;
		power = 1.0/gv.sph_ndim;
		for (int i=1; i<=par.GetSize(); i++)
			if (par[i].ToCompute() && par[i].p<par[i].pcut)
			{
				par[i].h = par[i].h0*(pow(par[i].rho0/par[i].rho, power));
				nv.sph_hmax = MAX(nv.sph_hmax, par[i].h);
			}
		MPI_Allreduce(&nv.sph_hmax, &hmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
		nv.sph_hmax = hmax;
		break;
	}

	return OK;
}
