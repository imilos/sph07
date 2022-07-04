#include "MSimulationData.h"
#include <math.h>
#include <stdio.h>

using namespace std;


// Constructor
MSimulationData::MSimulationData()
{
	pi = 4*atan(1.0);
}

// Destructor
MSimulationData::~MSimulationData()
{
}

MSimulationData* MSimulationData::_instance=0;
MSimulationData *MSimulationData::Instance()
{
	if (_instance == 0)
		_instance = new MSimulationData;
	return _instance;
}

void MSimulationData::DestroyInstance()
{
	delete _instance;
	_instance=0;
}

// Main startup routine
int MSimulationData::Startup(int argc, char **argv, bool &newproblem)
{
	string filename;
	MFileHandlingVars &fv = m_filevars;
	char tmp[300];

	if (argc!=1 && argc!=2 && argc!=4)
	{
		printf("\tError in subroutine startup\n");
		printf("\tWrong number of command line arguments.\n");
		return ERROR;
	}

	//  if no command line arguments are given then ask for them
	if (argc==1)
	{
		// new problem
		newproblem = true;
		printf(" Input file: ");
		scanf("%s", tmp);
		fv.sph_filein = tmp;
		printf(" Output file: ");
		scanf("%s", tmp);
		fv.sph_fileout = tmp;
		if (fv.sph_fileout.size()==0) fv.sph_fileout=fv.sph_filein;
	}
 	else
 	{
 		newproblem = true;
 		ParseCommandLine(argc, argv);
 	}

	// for a new problem
	if (newproblem)
	{
		if (fv.OpenInputFile()!=0) return ERROR;
		if (fv.OpenLogFile()!=0) return ERROR;
 		// write program and version information to log file
		if (my_id==ROOT_PROCESS) PrintVersion();
	}

	return OK;
}

/************************************************************************
		Purpose: Stop program cleanly
		Called by: Lots
		Date: 27-06-2007
		Errors: None
		Notes: exitcode 0 = normal termination
             exitcode 1 = error termination
*************************************************************************/
void MSimulationData::Shutdown(int exitcode)
{
	time_t now;
	float cputime;
	char realtime[30];
	MFileHandlingVars &fv=m_filevars;
	MGlobalVars &gv = m_globvars;

	if (my_id!=ROOT_PROCESS) return;

	//  write date and time of end of run to log file
	if (exitcode == OK)
	{
		printf("\tN o r m a l   T e r m i n a t i o n\n");
 		if (fv.f_logfile)
 			fprintf(fv.f_logfile, "\tN o r m a l   T e r m i n a t i o n\n");
	}
	else
	{
		printf("\tE R R O R   T E R M I N A T I O N\n");
 		if (fv.f_logfile)
 			fprintf(fv.f_logfile, "\tE R R O R   T E R M I N A T I O N\n");
	}

	now = time((time_t *)NULL);
	sprintf(realtime, "%s", ctime(&now));
	cputime = now - gv.sph_starttime;

	printf("At timestep: %6d      Problem time: %10.4f\n", gv.sph_timestep, gv.sph_ptime);
	printf("Total CPU time: %12.5e\n", cputime);

	if (fv.f_logfile)
	{
		fprintf(fv.f_logfile, "Analysis stopped: %s\n", realtime);
		fprintf(fv.f_logfile, "At timestep: %6d      Problem time: %10.4f\n", gv.sph_timestep, gv.sph_ptime);
		fprintf(fv.f_logfile, "Total CPU time: %12.5e\n", cputime);
	}
}

/************************************************************************
    Purpose: Get filenames from command line
    Called by: Startup()
    Date: 29-06-2007
    Errors:
    Notes:
************************************************************************/
int MSimulationData::ParseCommandLine(int argc, char **argv)
{
	m_filevars.sph_filein = argv[1];
	// Only temporary, velocity constraints for cylinders example
	if ( argc==4)
	{
		m_optvars.bConstrainVelocityCylinder = true;
		sscanf(argv[2], "%lf", &m_optvars.InnerW);
		sscanf(argv[3], "%lf", &m_optvars.OuterW);
	}

	return OK;
}

/************************************************************************
    Purpose: Prints program info
    Called by: Startup()
    Date: 29-06-2007
    Errors:
    Notes:
************************************************************************/

void MSimulationData::PrintVersion()
{
	MFileHandlingVars &fv = m_filevars;
	MGlobalVars &gv = m_globvars;

	// Console output of version information
	printf("%s\n%s\n%s\n%s\n%s\n",
		   "\tSPH07",
           "\tA 1D/2D/3D Smooth Particle Hydrodynamics Code, based on MCM, Cranfield University, UK",
           "\tUniversity of Kragujevac, Serbia",
           "\tVersion: 0.1",
		   "\tCode Date: 30-06-2007");

	// Write version and general information into log file
	fprintf(fv.f_logfile, "%s\n%s\n%s\n%s\n%s\n",
		   "\tSPH07",
           "\tA 1D/2D/3D Smooth Particle Hydrodynamics Code, based on MCM, Cranfield University, UK",
           "\tUniversity of Kragujevac, Serbia",
           "\tVersion: 0.1",
		   "\tCode Date: 30-06-2007");

	fprintf(fv.f_logfile, "\n\tInput filename: %16s.mcm\n", fv.sph_filein.c_str());
	fprintf(fv.f_logfile, "\tLog filename:   %16s.log\n", fv.sph_filein.c_str());
	fprintf(fv.f_logfile, "\tAnalysis started: %s\n", ctime(&gv.sph_starttime));
}

void MSimulationData::Fatal(const std::string &message)
{
	printf("\t%s\n", message.c_str());
	fprintf(m_filevars.f_logfile, "\t%s\n", message.c_str());
}

int MSimulationData::GetInput()
{
	string line, filename;
	char tmp[300];
	int version_number;
	MFileHandlingVars &fv = m_filevars;

 	// read and write  title problem title & input version
	PassInputComments(line);
	sscanf(line.c_str(), "%78c%2d", tmp, &version_number);
	tmp[78] = '\0'; fv.sph_title = tmp;

  	if ( ControlInput(version_number)!= OK) return ERROR;
  	//sph_initialise_global();
	if ( MaterialInput()!= OK ) return ERROR;
    if ( NodeInput()!= OK ) return ERROR;
    // Time history not implemented yet ////////////////////////
    //if (sph_thnode>0) sph_nodehist();
    // Transducers not implemented yet /////////////////////////
    //if(mcm_num_transducer>0) call mcm_transducers_input();
    if ( VelocityInput()!= OK ) return ERROR;
 	if (m_optvars.sph_boundary)
 		if ( BoundaryInput()!= OK ) return ERROR;
	// Contact not implemented yet ////////////////////////////
	//if (sph_contacttype>0) sph_contact_input();
	// Base accelerations
	if (m_baseaccelvars.sph_baseaccel) BaseAccelerationInput();
	if (PlotInput()!=OK) return ERROR;
	if (m_optvars.sph_rigid_body)
		if (RigidBodyInput()!=OK) return ERROR;
	if (LennardJonesInput()!=OK) return ERROR;

	return OK;
}

int MSimulationData::PassInputComments(string &ret)
{
	char *tmp, line[300];

	do
	{
		tmp = fgets(line, 300, m_filevars.f_inputfile );
	}
	while ( (line[0]=='*' || line[0]=='$') && tmp!=NULL);

	ret = line;

	if (tmp)
	{
		return OK;
	}
	else
	{
		printf("\tError reading input in routine sph_pass_comments\n");
		return ERROR;
	}
}

int MSimulationData::ControlInput(int version_number)
{
	int keep_files, mirr_opt, surface_tension_opt, rigid_body_opt;
	MGlobalVars &gv = m_globvars;
	MOptionVars &ov = m_optvars;
	MFileHandlingVars &fv = m_filevars;
	MOutputVars &outv = m_outvars;
	MSymPerVars &spv = m_sympervars;
	MBaseAccelVars &bav = m_baseaccelvars;
	std::string line;

	// default do not overwrite output files
	keep_files = 0;
	mirr_opt = 0;
	surface_tension_opt = 0;
	rigid_body_opt = 0;
	//
	// Read control cards
	//
	// Control card 1: problem definition
	PassInputComments(line);
	sscanf(line.c_str(), "%5d%5d%10d%10d%10d%10d%10d", &gv.sph_axopt,&gv.sph_disctype,
		&gv.sph_nummat,&gv.sph_np,&gv.sph_max_np, &gv.sph_nstressp,&gv.sph_nvelocp);

	// Control card 2: time control
	PassInputComments(line);
	sscanf(line.c_str(), "%10lf%10lf%10lf%10lf", &gv.sph_endtime,&gv.sph_tssf,&gv.sph_itss,&ov.sph_drelax_scale);

	// Control card 3: output file control
	PassInputComments(line);
	sscanf(line.c_str(), "%10lf%5d%10lf%10d%5d%10d%5d%5d", &outv.sph_stpltime,&outv.sph_state_opt,
		&outv.sph_thpltime,&outv.sph_thnode,&outv.sph_num_transducer,&gv.sph_status_interval, &gv.sph_restart_interval, &gv.sph_run_restart);

    // Control card 4: input and initialization options
	PassInputComments(line);
	sscanf(line.c_str(), "%5d%5d%5d%5d%5d", &keep_files,&ov.sph_init_v_opt,&ov.sph_massopt,&ov.sph_init_h_opt,&ov.sph_init_rhoe);
    //
    // Control card 5: analysis options
	PassInputComments(line);
	sscanf(line.c_str(), "%5d%5d%5d%5d%5d%5d%5d%5d%5d%5d", &ov.sph_contacttype, &ov.sph_tcrit_opt, &ov.sph_veloc_opt,
			&mirr_opt, &ov.sph_nlcur, &ov.sph_nthpx, &ov.sph_nthpy, &ov.sph_nthpz, &surface_tension_opt, &rigid_body_opt);
	//
	// Control card 6: Interpolation options
	PassInputComments(line);
	sscanf(line.c_str(), "%5d%5d%5d", &ov.sph_h_opt,&ov.sph_krtype,&gv.sph_maxnbr);
	//
	// Control card 7: Empty at this time, blank line
	PassInputComments(line);
	//sscanf(line.c_str(), "", );

	//
	//=====================================================================================
	// Initialize values if required
	//
	if (surface_tension_opt!=0) ov.sph_surface_tension = true;
	if (rigid_body_opt!=0) ov.sph_rigid_body = true;
	if (gv.sph_tssf==0.0) gv.sph_tssf=0.8;
	if (ov.sph_krtype==0) ov.sph_krtype=1;
	if (gv.sph_maxnbr==0) gv.sph_maxnbr=40;
	if (gv.sph_status_interval<1) gv.sph_status_interval = 100;
	if (mirr_opt==1) ov.sph_boundary = true;
	spv.sph_g_maxnbr = gv.sph_maxnbr/2;

	// If maximum number of particles is smaller than np
	if (gv.sph_max_np<gv.sph_np) gv.sph_max_np = gv.sph_np;
	// Set number of dimensions
	gv.sph_ndim=gv.sph_axopt;
	// Dynamic realxation
	if (ov.sph_drelax_scale>0.1 && ov.sph_drelax_scale<1.0) ov.sph_drelax = true;
	// Set up loop variables for discretisation options
	// Collocated SPH is the only option so far
	if (gv.sph_disctype==0 || gv.sph_disctype==1 || gv.sph_disctype==2)
	{
	 	gv.sph_ssp=1;
	  	gv.sph_svp=1;
	  	gv.sph_esp=gv.sph_np;
	  	gv.sph_evp=gv.sph_np;
	}
	// Check for base accelerations
	if (ov.sph_nthpx==1 || ov.sph_nthpy==1 || ov.sph_nthpz==1) bav.sph_baseaccel=true;

	//=====================================================================================
	// Write control data into log file
	//
	// Problem title and version number
	fprintf(fv.f_logfile,"\n%s\n%78s%2d\n%s\n",
	"**********************************************************************",
	fv.sph_title.c_str(), version_number,
	"**********************************************************************");

	fprintf(fv.f_logfile,"\n\n%s\n\n", "\tCONTROL DATA");

	// Control card 1
	fprintf(fv.f_logfile,"\n%s\n%s\n%s\n",
	"**********************************************************************",
	"*               Control card 1                                       *",
	"**********************************************************************");

	fprintf(fv.f_logfile, "\n%s%7d\n%s%7d\n%s%11d\n%s%11d\n%s%11d\n%s%7d\n",
	  "axis option.................................    ",gv.sph_axopt,
      "discretisation type.........................    ",gv.sph_disctype,
      "number of materials.........................",gv.sph_nummat,
      "number of nodes.............................",gv.sph_np,
	  "maximum number of nodes.....................",gv.sph_max_np,
      "number of dimensions........................    ",gv.sph_ndim);

	// Control card 2
	fprintf(fv.f_logfile,"\n%s\n%s\n%s\n",
	"**********************************************************************",
	"*               Control card 2                                       *",
	"**********************************************************************");

	fprintf(fv.f_logfile, "\n%s%10.2e\n%s%10.2e\n%s%10.2e\n%s\n",
      "termination time............................ ",gv.sph_endtime,
      "time step scale factor...................... ",gv.sph_tssf,
      "initial time step size...................... ",gv.sph_itss,
      "\t==0.0,  program picks initial step size        ");

    if (ov.sph_drelax)
    	fprintf(fv.f_logfile, "%s%10.2e\n",
    		"Dynamic relaxation active with scale factor  ", ov.sph_drelax_scale);
   	else
    	fprintf(fv.f_logfile, "%s\n", "Dynamic relaxation not active.");

	// Control card 3
	fprintf(fv.f_logfile,"\n%s\n%s\n%s\n",
	"**********************************************************************",
	"*               Control card 3                                       *",
	"**********************************************************************");

	fprintf(fv.f_logfile, "\n%s%10.2e\n%s%10d\n%s%10.2e\n%s%10d\n%s%5d\n%s%10d\n%s%10d\n%s%10d\n",
      "time interval between state plots........... ",outv.sph_stpltime,
	  "output file fomat........................... ",outv.sph_state_opt,
      "time interval between time history plots.... ",outv.sph_thpltime,
      "number of time history nodes................ ",outv.sph_thnode,
	  "Number of pressure transducers..............      ",outv.sph_num_transducer,
	  "Number of steps between status reports...... ",gv.sph_status_interval,
	  "Number of steps between restart files....... ",gv.sph_restart_interval,
	  "Number of steps between running restart files",gv.sph_run_restart);

	// Control card 4
	fprintf(fv.f_logfile,"\n%s\n%s\n%s\n",
	"**********************************************************************",
	"*               Control card 4                                       *",
	"**********************************************************************");

	fprintf(fv.f_logfile,"\n%s%10d\n%s\n%s%10d\n%s%10d\n%s%10d\n%s%10d\n",
      "Preserve output files flag.................. ",keep_files,
	  "\t==1, output files will be overwritten      ",
	  "Prescribed initial velocity flag............ ",ov.sph_init_v_opt,
	  "Mass initialisation option.................. ",ov.sph_massopt,
	  "Smoothing length initialisation option...... ",ov.sph_init_h_opt,
	  "Particle density and energy initialisation   ",ov.sph_init_rhoe );

	// Control card 5
	fprintf(fv.f_logfile,"\n%s\n%s\n%s\n",
	"**********************************************************************",
	"*               Control card 5                                       *",
	"**********************************************************************");

	fprintf(fv.f_logfile,"\n%s%7d\n%s%7d\n%s%7d\n%s%7d\n%s%7d\n%s%7d\n%s\n%s\n%s%7d\n%s\n%s\n%s%7d\n%s\n%s\n",
      "contact type................................    ",ov.sph_contacttype,
	  "critical timestep calculation option........    ",ov.sph_tcrit_opt,
	  "Velocity smoothing option...................    ",ov.sph_veloc_opt,
      "Symmetry planes flag........................    ",mirr_opt,
	  "Number of load curves.......................    ",ov.sph_nlcur,
      "X-dir base acceleration.....................    ",ov.sph_nthpx,
      "\t==0,  no                                       ",
      "\t==1,  yes                                      ",
      "Y-dir base acceleration.....................    ",ov.sph_nthpy,
      "\t==0,  no                                       ",
      "\t==1,  yes                                      ",
      "Z-dir base acceleration.....................    ",ov.sph_nthpz,
      "\t==0,  no                                       ",
      "\t==1,  yes                                      ");

	// Control card 6
	fprintf(fv.f_logfile,"\n%s\n%s\n%s\n",
	"**********************************************************************",
	"*               Control card 6                                       *",
	"**********************************************************************");

	fprintf(fv.f_logfile,"\n%s%7d\n%s%11d\n%s%7d\n",
      "Variable smoothing length option............    ",ov.sph_h_opt,
      "kernel type.................................",ov.sph_krtype,
	  "Maximum number of neighbours................    ",gv.sph_maxnbr);

	// Control card 7
	fprintf(fv.f_logfile,"\n%s\n%s\n%s\n",
	"**********************************************************************",
	"*               Control card 7                                       *",
	"**********************************************************************");
	fprintf(fv.f_logfile,"\nCard is not used in this version.\n");

	return OK;
}

int MSimulationData::MaterialInput()
{
	MGlobalVars &gv= m_globvars;
	int m, ihq, j;
	std::string line, msg;
	real qh;

	// Pass throough all materials
	for (m=1; m<=gv.sph_nummat; m++)
	{
		MMaterial tmp;

		PassInputComments(line);
		sscanf(line.c_str(), "%5d%5d%10lf%5d%5d%10lf%5d%10lf%10lf%10lf",&tmp.n,&tmp.model,&tmp.rho,
			&tmp.eos,&ihq,&qh,&tmp.visc_type,&tmp.av_q,&tmp.av_l, &tmp.surfaceTensionCoeff);

		// Read material options card
		PassInputComments(line);
		sscanf(line.c_str(), "%10lf%10lf%10lf%10lf",&tmp.mass,&tmp.h,&tmp.rho_min,&tmp.rho_max);

		// Set default quantities if necessary
		tmp.CheckDefaults();

		// if supplied material number is outside the allowed range
		if (tmp.n<0 || tmp.n>gv.sph_nummat)
		{
			Fatal("Error in input - material number is out of range.");
			return ERROR;
		}

		// read material title card
		PassInputComments(line);
		sscanf(line.c_str(), "%6s%6s%6s%6s%6s%6s%6s%6s%6s%6s%6s%6s",
			tmp.head[0][0],tmp.head[0][1],tmp.head[0][2],tmp.head[0][3],
			tmp.head[0][4],tmp.head[0][5],tmp.head[0][6],tmp.head[0][7],
			tmp.head[0][8],tmp.head[0][9],tmp.head[0][10],tmp.head[0][11]);

	 	// What to be read depends on material model
	 	switch (tmp.model)
	 	{
		 	case 1: case 3: case 9:
		   		for (j=0; j<6; j++)
		   		{
		   			PassInputComments(line);
		   			if (j==0)
		   				sscanf(line.c_str(),"%10lf%10lf%10lf%10lf%10lf%10lf%10lf%10lf",
						&tmp.strinput(0),&tmp.strinput(1),&tmp.strinput(2),&tmp.strinput(3),
						&tmp.strinput(4),&tmp.strinput(5),&tmp.strinput(6),&tmp.strinput(7));
		   		}
		   		break;

		   	case 10:
		   		for (j=0; j<6; j++)
		   		{
		   			PassInputComments(line);
		   			sscanf(line.c_str(), "%10lf%10lf%10lf%10lf%10lf%10lf%10lf%10lf",
						&tmp.strinput(0),&tmp.strinput(1),&tmp.strinput(2),&tmp.strinput(3),
						&tmp.strinput(4),&tmp.strinput(5),&tmp.strinput(6),&tmp.strinput(7));
		   		}
		   		break;
	 	}

	 	// Read equation of state title
	 	if (tmp.eos != 0)
	 	{
			// Check EOS support
			if ( tmp.CheckMaterialModel(msg) != OK )
			{
				Fatal(msg);
				return ERROR;
			}

	 		PassInputComments(line);
	 		sscanf(line.c_str(), "%6s%6s%6s%6s%6s%6s%6s%6s%6s%6s%6s%6s",
	 			tmp.head[1][0],tmp.head[1][1],tmp.head[1][2],tmp.head[1][3],
	 			tmp.head[1][4],tmp.head[1][5],tmp.head[1][6],tmp.head[1][7],
	 			tmp.head[1][8],tmp.head[1][9],tmp.head[1][10],tmp.head[1][11]);
	 	}

		//  read in equation of state constants,
		switch (tmp.eos)
		{
			case 28:
				PassInputComments(line);
			    sscanf(line.c_str(), "%10lf%10lf", &tmp.eosinput(0), &tmp.eosinput(1));
			    break;
			case 29:
				PassInputComments(line);
			    sscanf(line.c_str(), "%10lf", &tmp.eosinput(0));
			    break;
			case 13:
				PassInputComments(line);
			    sscanf(line.c_str(), "%10lf%10lf%10lf%10lf", &tmp.eosinput(0),&tmp.eosinput(1),
			    		&tmp.eosinput(2),&tmp.eosinput(3));
			    break;
		}

		// write material data into log file
		PrintMaterialProperties(tmp);
		// Add read material into material data
		mat.AddMaterial(tmp);
	} // material loop

	return OK;
}

// Print out material properties into log file
int MSimulationData::PrintMaterialProperties(MMaterial &tmp)
{
	int j;
	MFileHandlingVars &fv = m_filevars;
	MOptionVars &ov = m_optvars;
	//
	// If first material then print header into log
	// default viscosity values hardwired
	//
	if (tmp.n==1)
	{
		fprintf(fv.f_logfile, "\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n",
				 " m a t e r i a l   d e f i n i t i o n ",
                 " material models                       ",
                 "     ==1   isotropic elastic          ",
                 "     ==3   kinematic/isotropic elastic-plastic ",
                 "     ==9   fluid material              ",
                 "     ==10  isotropic elastic-plastic hydrodynamic",
                 "     ==15  johnson/cook elastic-plastic ");

		fprintf(fv.f_logfile, "\n%s\n%s\n%s\n%s\n%s\n%s\n",
                 " equation-of-state models             ",
                 "     ==4   gruneisen                 ",
                 "     ==13  perfect gas               ",
				 "     ==28  monaghan  (quasi-incompr.)",
				 "     ==29  morris    (quasi-incompr.)",
				 "     ==41  mie-gruneisen             ");

		fprintf(fv.f_logfile, "\n%s\n%s\n",
                 " bulk viscosity models                ",
                 "     ==1   standard                  ");

		fprintf(fv.f_logfile, "\n%s\n%s%5d\n%s%12.4e\n%s%12.4e\n",
                 " default viscosities                  ",
                 "     bulk viscosity type...................",1,
                 "     linear bulk viscosity coefficient.....",0.06,
                 "     quadratic bulk viscosity coefficient..",1.5);
	}

	// Head output
	for (j=0; j<12; j++)
		fprintf(fv.f_logfile, "%6s", tmp.head[0][j]);
	fprintf(fv.f_logfile, "\n");

	// General material data
	fprintf(fv.f_logfile, "\n%s%5d\n%s%5d\n%s%5d\n%s%5d\n%s%12.4e\n%s%12.4e\n%s%12.4e\n%s%12.4e\n",
		" material constants set number ........ ",tmp.n,
        "\t\tmaterial model ............. ",tmp.model,
        "\t\tequation-of-state model .... ",tmp.eos,
        "\t\tbulk viscosity model ....... ",tmp.visc_type,
        "\tden .............................. =", tmp.rho,
        "\tquadratic bulk viscosity ......... =", tmp.av_q,
        "\tlinear bulk viscosity ............ =", tmp.av_l,
        "\tsurface tension coefficient ...... =", tmp.surfaceTensionCoeff);

	//
	// write initial total material mass and h if given
	//
	if (ov.sph_massopt == 0)
		fprintf(fv.f_logfile, "\n%s%12.4e\n",
			"\tTotal material mass............... =", tmp.mass);

	if (ov.sph_init_h_opt == 0)
		fprintf(fv.f_logfile, "\n%s%12.4e\n",
			"\tInitial h for material............ =", tmp.h);

	fprintf(fv.f_logfile, "\n%s%12.4e\n%s%12.4e\n",
		"\tMin density limit factor for mat.. =", tmp.rho_min,
        "\tMax density limit factor for mat.. =", tmp.rho_max);

    //
    // Print material properties and constants
    //
    switch (tmp.model)
    {
    	// model == 9, FLUID
    	case 9:
			fprintf(fv.f_logfile, "\n%s%12.4e\n%s%12.4e\n",
        		"\tpressure cutoff .................. =", tmp.strinput(0),
        		"\tviscosity coefficient ............ =", tmp.strinput(1));
		//
		// other material models - to be implemented
		//
    }

	// Equation of state header
	if (tmp.eos > 0)
		for (j=0; j<12; j++) fprintf(fv.f_logfile, "%s", tmp.head[1][j]);
	fprintf(fv.f_logfile, "\n");

	// Equation of state description
	switch (tmp.eos)
	{
		case 28:
			fprintf(fv.f_logfile, "\n%s%12.4e\n%s%12.4e\n",
	        	"\tB ................................ =", tmp.eosinput(0),
	        	"\tGamma ............................ =", tmp.eosinput(1));
			break;
		case 29:
			fprintf(fv.f_logfile, "\n%s%12.4e\n",
	        	"\tC ................................ =", tmp.eosinput(0));
			break;
	}

	return OK;
}

int MSimulationData::NodeInput()
{
	MGlobalVars &gv = m_globvars;
	MFileHandlingVars &fv = m_filevars;
	MOptionVars &ov = m_optvars;
	int j;
	std::string line;
	char msg[300];
	real pardata[4];

	for (j=1; j<=gv.sph_np; j++)
	{
		MParticle tmp;

		PassInputComments(line);
		// Number of dimensions determines the format
		switch (gv.sph_ndim)
		{
			case 1:
				sscanf(line.c_str(), "%8d%5d%20lf%7d%7d", &tmp.n,&tmp.dispbc,&tmp.x(0),&tmp.mat,&tmp.color);
				break;
			case 2:
				sscanf(line.c_str(), "%8d%5d%20lf%20lf%7d%7d", &tmp.n,&tmp.dispbc,&tmp.x(0),&tmp.x(1),&tmp.mat,&tmp.color);
				break;
			case 3:
				sscanf(line.c_str(), "%8d%5d%20lf%20lf%20lf%7d%7d", &tmp.n,&tmp.dispbc,&tmp.x(0),&tmp.x(1),&tmp.x(2),&tmp.mat,&tmp.color);
				break;
		}

		// Check node ordering, zero based
		if (j != tmp.n)
		{
			sprintf(msg, "\tCard for node %8d is out of order. Should be %8d\n", tmp.n, j);
			Fatal(msg);
			return ERROR;
		}

		// Check if material is out of range
		if (tmp.mat==LENNARDJONES) {
			tmp.lennardJones = true;
			// Converts Lennard-Jones particle into first material, it must be there
			tmp.mat = 1;
		}
		else if (tmp.mat<1 || tmp.mat>gv.sph_nummat) {
				sprintf(msg,"Node %8d has incorrect material type.", tmp.n);
				Fatal(msg);
				return ERROR;
		}

		// Add particle to array
		par.AddParticle(tmp);
	} // node loop

	//
	// now write out nodal coordinates to log file
	//
	switch (gv.sph_disctype)
	{
		// Colocated
		case 0: case 1: case 2:
			fprintf(fv.f_logfile, "\nN O D E   I N P U T   C A R D S\n\n");
			// Number of dimensions
			switch (gv.sph_ndim)
			{
				case 1:
					 fprintf(fv.f_logfile, " Node id Node bc    x coordinate Node material\n");
					 for (j=1; j<=gv.sph_np; j++)
					 	fprintf(fv.f_logfile,"%8d%8d%16.8e%5d\n",j,par[j].dispbc,par[j].x(0),par[j].mat);
					 break;

				case 2:
					 fprintf(fv.f_logfile, " Node id Node bc    x coordinate    y coordinate Node material\n");
					 for (j=1; j<=gv.sph_np; j++)
					 	fprintf(fv.f_logfile,"%8d%8d%16.8e%16.8e%5d\n",j,par[j].dispbc,par[j].x(0),par[j].x(1),par[j].mat);
					 break;

				case 3:
					 fprintf(fv.f_logfile, " Node id Node bc    x coordinate    y coordinate    z coordinate Node material\n");
					 for (j=1; j<=gv.sph_np; j++)
					 	fprintf(fv.f_logfile,"%8d%8d%16.8e%16.8e%16.8e%5d\n",j,par[j].dispbc,par[j].x(0),par[j].x(1),par[j].x(2),par[j].mat);
					 break;
			}
	}
	//
	// If additional node cards are present
	//
	if (ov.sph_massopt==1 || ov.sph_init_h_opt==1 || ov.sph_init_rhoe==1)
	{
		for (j=1; j<=gv.sph_np; j++)
		{
			PassInputComments(line);
			sscanf(line.c_str(), "%8d%20lf%20lf%20lf%20lf", &par[j].n,&pardata[0],&pardata[1],&pardata[2],&pardata[3]);
			// Check node ordering, zero based
			if (j != par[j].n)
			{
				sprintf(msg,"\tAdditional card for node %8d is out of order. Should be %8d\n", par[j].n, j+1);
				Fatal(msg);
				return ERROR;
			}
			// Copy data into structure
			if (ov.sph_massopt==1) par[j].mass = pardata[0];
			if (ov.sph_init_h_opt==1) par[j].h = pardata[1];
			if (ov.sph_init_rhoe==1)
			{
				par[j].rho = pardata[2];
				par[j].e = pardata[3];
			}
		}
	} // if clause for additional node cards

	return OK;
}

int MSimulationData::VelocityInput()
{
	int i, j, nid;
	bool genflag;
	std::string line;
	char msg[300];
	MOptionVars &ov = m_optvars;
	MGlobalVars &gv = m_globvars;
	MFileHandlingVars &fv = m_filevars;
	// Start and end velocity points
	int svp = gv.sph_svp;
	int evp = gv.sph_evp;

	// If velocities are not defined in input file return
	if (ov.sph_init_v_opt!=1) return OK;

	PassInputComments(line);

	// Read initial velocities from the input file
	switch (gv.sph_ndim)
	{
		case 1:
			sscanf(line.c_str(), "%8d%10lf", &nid,&par[svp].v(0));
			break;
		case 2:
			sscanf(line.c_str(), "%8d%10lf%10lf", &nid,&par[svp].v(0),&par[svp].v(1));
			break;
		case 3:
			sscanf(line.c_str(), "%8d%10lf%10lf%10lf", &nid,&par[svp].v(0),&par[svp].v(1),&par[svp].v(2));
			break;
	}

	// error if first velocity is not for first velocity node
	if (nid != svp)
	{
		sprintf(msg,"\tError, first nodal velocity card is not for node %d\n", svp);
		Fatal(msg);
		return ERROR;
	}

	genflag = false;
	// One based index
	for (i=svp+1; i<=evp; i++)
	{
		// The particles are either copied or generated
		if (!genflag)
		{
			PassInputComments(line);

			// Read velocity line
			switch (gv.sph_ndim)
			{
				case 1: sscanf(line.c_str(),"%8d%10lf",&nid,&par[i].v(0)); break;
				case 2: sscanf(line.c_str(),"%8d%10lf%10lf",&nid,&par[i].v(0),&par[i].v(1)); break;
				case 3: sscanf(line.c_str(),"%8d%10lf%10lf%10lf",&nid,&par[i].v(0),&par[i].v(1),&par[i].v(2)); break;
			}
  			if (i < nid)
			{
				// Error if non existant node
				if (nid > gv.sph_np)
				{
					Fatal("\tParticle id in velocity input exceeds number of particles\n");
					return ERROR;
				}
				// Error if velocities are not constant in range
				if ( par[i].v(0)!=par[i-1].v(0) || par[i].v(1)!=par[i-1].v(1) ||
					par[i].v(2)!=par[i-1].v(2) )
				{
					Fatal("\tParticle velocity for auto generation not constant.\n");
					return ERROR;
				}
				genflag=true;
			}
			// Error: particle order incorrect
			else if (i > nid)
			{
				Fatal("\tVelocities must be specified in increasing order.\n");
				return ERROR;
			}
		}
		else // genflag = true
		{
			// if incrementing then assign correct velocity to the node
			for (j=0; j<gv.sph_ndim; j++) par[i].v(j)=par[i-1].v(j);
			if (i==nid) genflag=false;
		}
	} // velocity point loop

	// Writting velocities into log file
	fprintf(fv.f_logfile, "\nI N I T I A L   P A R T I C L E  V E L O C I T I E S\n\n");

	switch (gv.sph_ndim)
	{
		case 1:
			fprintf(fv.f_logfile," Node id      x velocity\n");
			for (i=svp; i<=evp; i++)
				fprintf(fv.f_logfile,"%8d%16.8e\n",i,par[i].v(0));
			break;

		case 2:
			fprintf(fv.f_logfile," Node id      x velocity      y velocity\n");
			for (i=svp; i<=evp; i++)
				fprintf(fv.f_logfile,"%8d%16.8e%16.8e\n",i,par[i].v(0),par[i].v(1));
			break;

		case 3:
			fprintf(fv.f_logfile," Node id      x velocity      y velocity      z velocity\n");
			for (i=svp; i<=evp; i++)
				fprintf(fv.f_logfile,"%8d%16.8e%16.8e%16.8e\n",i,par[i].v(0),par[i].v(1),par[i].v(2));
			break;
	}

	return OK;
}

int MSimulationData::BoundaryInput()
{
	std::string line;
	bool found;
	int i;
	MGlobalVars &gv = m_globvars;
	MSymPerVars &spv = m_sympervars;
	MOptionVars &ov = m_optvars;
	MFileHandlingVars &fv = m_filevars;

	// Read symmetry flags
	PassInputComments(line);
	switch (gv.sph_ndim)
	{
		case 1:
			sscanf(line.c_str(),"%10d%10d", &spv.sph_boundary_code(0,0),&spv.sph_boundary_code(1,0));
			break;
		case 2:
			sscanf(line.c_str(),"%10d%10d%10d%10d", &spv.sph_boundary_code(0,0),&spv.sph_boundary_code(1,0),&spv.sph_boundary_code(0,1),&spv.sph_boundary_code(1,1));
			break;
		case 3:
			sscanf(line.c_str(),"%10d%10d%10d%10d%10d%10d",
				&spv.sph_boundary_code(0,0),&spv.sph_boundary_code(1,0),&spv.sph_boundary_code(0,1),&spv.sph_boundary_code(1,1),&spv.sph_boundary_code(0,2),&spv.sph_boundary_code(1,2));
			break;
	}

	// Read symmetry planes (lines)
	PassInputComments(line);
	switch (gv.sph_ndim)
	{
		case 1:
			sscanf(line.c_str(),"%10lf%10lf", &spv.sph_boundary_x(0,0),&spv.sph_boundary_x(1,0));
			break;
		case 2:
			sscanf(line.c_str(),"%10lf%10lf%10lf%10lf", &spv.sph_boundary_x(0,0),&spv.sph_boundary_x(1,0),&spv.sph_boundary_x(0,1),&spv.sph_boundary_x(1,1));
			break;
		case 3:
			sscanf(line.c_str(),"%10lf%10lf%10lf%10lf%10lf%10lf",
				&spv.sph_boundary_x(0,0),&spv.sph_boundary_x(1,0),&spv.sph_boundary_x(0,1),&spv.sph_boundary_x(1,1),&spv.sph_boundary_x(0,2),&spv.sph_boundary_x(1,2));
			break;
	}

	// Error handling

	// Check for boundary condition type
	// Type 1: symmetry
	for (i=0,found=false; i<3; i++)
		if (spv.sph_boundary_code(0,i) || spv.sph_boundary_code(1,i)) found=true;

	// Type 2: periodic
	for (i=0; i<3; i++)
		if (spv.sph_boundary_code(0,i)==2)
		{
			found = true;
			if (spv.sph_boundary_code(1,i)!=2)
			{
				Fatal("Error in definition of periodic boundary conditions.");
				return ERROR;
			}
		}

	for (i=0; i<3; i++)
		if (spv.sph_boundary_code(1,i)==2)
		{
			found = true;
			if (spv.sph_boundary_code(0,i)!=2)
			{
				Fatal("Error in definition of periodic boundary conditions.");
				return ERROR;
			}
		}

	// If all codes are zero, assure that spv.sph_boundary == 0
	if (!found) ov.sph_boundary = false;

	// Write input to log file
	fprintf(fv.f_logfile, "\nB O U N D A R Y   P L A N E S\n\n");

	// Symmetry planes
	if (spv.sph_boundary_code(0,0)==1)
		fprintf(fv.f_logfile, " Xmin symmetry plane active at x coordinate:%11.3e\n", spv.sph_boundary_x(0,0));
	if (spv.sph_boundary_code(1,0)==1)
		fprintf(fv.f_logfile, " Xmax symmetry plane active at x coordinate:%11.3e\n", spv.sph_boundary_x(1,0));
	if (spv.sph_boundary_code(0,1)==1)
		fprintf(fv.f_logfile, " Ymin symmetry plane active at y coordinate:%11.3e\n", spv.sph_boundary_x(0,1));
	if (spv.sph_boundary_code(1,1)==1)
		fprintf(fv.f_logfile, " Ymax symmetry plane active at y coordinate:%11.3e\n", spv.sph_boundary_x(1,1));
	if (spv.sph_boundary_code(0,2)==1)
		fprintf(fv.f_logfile, " Zmin symmetry plane active at z coordinate:%11.3e\n", spv.sph_boundary_x(0,2));
	if (spv.sph_boundary_code(1,2)==1)
		fprintf(fv.f_logfile, " Zmax symmetry plane active at z coordinate:%11.3e\n", spv.sph_boundary_x(1,2));

	// Periodic planes
	if (spv.sph_boundary_code(0,0)==2)
		fprintf(fv.f_logfile, " Xmin periodic plane active at x coordinate:%11.3e\n", spv.sph_boundary_x(0,0));
	if (spv.sph_boundary_code(1,0)==2)
		fprintf(fv.f_logfile, " Xmax periodic plane active at x coordinate:%11.3e\n", spv.sph_boundary_x(1,0));
	if (spv.sph_boundary_code(0,1)==2)
		fprintf(fv.f_logfile, " Ymin periodic plane active at y coordinate:%11.3e\n", spv.sph_boundary_x(0,1));
	if (spv.sph_boundary_code(1,1)==2)
		fprintf(fv.f_logfile, " Ymax periodic plane active at y coordinate:%11.3e\n", spv.sph_boundary_x(1,1));
	if (spv.sph_boundary_code(0,2)==2)
		fprintf(fv.f_logfile, " Zmin periodic plane active at z coordinate:%11.3e\n", spv.sph_boundary_x(0,2));
	if (spv.sph_boundary_code(1,2)==2)
		fprintf(fv.f_logfile, " Zmax periodic plane active at z coordinate:%11.3e\n", spv.sph_boundary_x(1,2));

	// 1-both free, 2-min fixed, 3-max fixed, 4-both fixed
	for (i=0; i<3; i++)
	{
		if (spv.sph_boundary_code(0,i)>0) spv.sph_boundary_type(i) += 1;
		if (spv.sph_boundary_code(1,i)>0) spv.sph_boundary_type(i) += 2;
	}

	return OK;
}

int MSimulationData::BaseAccelerationInput()
{
	std::string line;
	MOptionVars &ov = m_optvars;
	MBaseAccelVars &bav = m_baseaccelvars;

	if (bav.sph_baseaccel)
	{
		PassInputComments(line);
		sscanf(line.c_str(), "%20lf%20lf%20lf",
			&bav.sph_base_a(0), &bav.sph_base_a(1), &bav.sph_base_a(2) );
	}

	if (ov.sph_nthpx != 1) bav.sph_base_a(0) = 0.0;
	if (ov.sph_nthpy != 1) bav.sph_base_a(1) = 0.0;
	if (ov.sph_nthpz != 1) bav.sph_base_a(2) = 0.0;

	return OK;
}

int MSimulationData::PlotInput()
{
	std::string line;
	MOutputVars &outv = m_outvars;

	PassInputComments(line);
	sscanf(line.c_str(), "%10lf", &outv.sph_plottime);

	if (outv.sph_plottime != 0.0)
	{
		PassInputComments(line);
		sscanf(line.c_str(), "%20lf%20lf%20lf",
			&outv.sph_plot_start_x(0), &outv.sph_plot_start_x(1), &outv.sph_plot_start_x(2) );
		PassInputComments(line);
		sscanf(line.c_str(), "%20lf%20lf%20lf",
			&outv.sph_plot_end_x(0), &outv.sph_plot_end_x(1), &outv.sph_plot_end_x(2) );
	}

	return OK;
}

int MSimulationData::LennardJonesInput()
{
	std::string line;
	MGlobalVars &gv = m_globvars;
	MOptionVars &ov = m_optvars;

	// if any particle is of lennard-jones type?
	for (int i=1; i<gv.sph_np; ++i)
		if (par[i].lennardJones) ov.sph_lennardjones=true;

	if (ov.sph_lennardjones) {
		PassInputComments(line);
		sscanf(line.c_str(), "%10lf%10lf", &gv.lennardJonesD, &gv.lennardJonesR0);
	}

	return OK;
}

int MSimulationData::RigidBodyInput()
{
	std::string line;
	MRigidBodyData &rbd = m_rigidbodydata;

	PassInputComments(line);
	sscanf(line.c_str(), "%5d%10lf%10lf%10lf%10lf", &rbd.ID, &rbd.totalMass, &rbd.inertiaMoment,
			&rbd.particleMass, &rbd.deltap);

	return OK;
}
