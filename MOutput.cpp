#include "MOutput.h"
#include <vector>
#include <stdio.h>

MOutput::MOutput(MSimulationData *simdata)
{
	m_simdata = simdata;
	#ifdef MICROSOFT_PATCH //added for correcting Windows output
	unsigned int old_exponent_format;
	old_exponent_format = _set_output_format(_TWO_DIGIT_EXPONENT);
    #endif
}

MOutput::~MOutput()
{
}

MOutput* MOutput::_instance=0;
MOutput *MOutput::Instance(MSimulationData *simdata)
{
	if (_instance == 0)
		_instance = new MOutput(simdata);
	return _instance;
}

void MOutput::DestroyInstance()
{
	delete _instance;
	_instance=0;
}


void MOutput::PrintScreenLog(const std::string &msg)
{
	printf("%s\n", msg.c_str());
	fprintf(m_simdata->m_filevars.f_logfile, "%s\n", msg.c_str());
}

int MOutput::InitStateOutput()
{
	MSimulationData &sd = *m_simdata;
	MGlobalVars &gv = sd.m_globvars;
	MOutputVars &outv = sd.m_outvars;

	switch (outv.sph_state_opt)
	{
		// Ensight 5.0 - not implemented yet
		case 1:
			PrintScreenLog("ERROR: Ensight 5.0 format not supported yet.");
			return ERROR;
			break;
		// Ensight CASE format
		case 2:
			// Calc time for second state plot
			outv.sph_nextsttime = gv.sph_ptime + outv.sph_stpltime;
			if (WriteCaseGeo()!=OK) return ERROR;
			break;
		// LS-DYNA format - not implemented yet
		case 3:
			PrintScreenLog("ERROR: LS-DYNA format not supported yet.");
			return ERROR;
			break;
		// MCMGUI format - not implemented yet
		default:
			PrintScreenLog("ERROR: MCMGUI format not supported yet.");
			return ERROR;
			break;
	}

	return OK;
}

int MOutput::StateOutput()
{
	MSimulationData &sd = *m_simdata;
	MGlobalVars &gv = sd.m_globvars;
	MOutputVars &outv = sd.m_outvars;

	switch (outv.sph_state_opt)
	{
		// Ensight CASE
		case 2:
			switch (gv.sph_ndim)
			{
				case 1:
#ifdef PARAVIEW_PATCH
					WriteCaseGeo();
#endif
					WriteCase();
					WriteCase1D();
					IncreaseStateCounter();
					break;
				case 2:
// Due to the known Paraview bug, the geometry has to be written in each time step
#ifdef PARAVIEW_PATCH
					WriteCaseGeo();
#endif
					WriteCase();
					WriteCase2D();
					IncreaseStateCounter();
					break;
				case 3:
#ifdef PARAVIEW_PATCH
					WriteCaseGeo();
#endif
					WriteCase();
					WriteCase3D();
					IncreaseStateCounter();
					break;
			} break;
	}

	return OK;
}

int MOutput::WriteCase()
{
	MSimulationData &sd = *m_simdata;
	MFileHandlingVars &fv = sd.m_filevars;
	MGlobalVars &gv = sd.m_globvars;
	MOutputVars &outv = sd.m_outvars;
	const std::string &base = fv.sph_filein;

	// status="0" means to open new file
	if (fv.OpenOutputFile(".case", &fv.f_case, 0)!=OK) return ERROR;

	fprintf(fv.f_case, "FORMAT\n");
	fprintf(fv.f_case, "type:   ensight\n\n");
	fprintf(fv.f_case, "GEOMETRY\n");
	fprintf(fv.f_case, "model:                                    %26s\n\n", (base+".geo").c_str());
	fprintf(fv.f_case, "VARIABLES\n");
	fprintf(fv.f_case, "scalar per node:         1 1 mass         %26s\n", (base+"_mas.ens").c_str());
	fprintf(fv.f_case, "scalar per node:         1 1 boundary     %26s\n", (base+"_bnd.ens").c_str());
	fprintf(fv.f_case, "scalar per node:         1 1 nnbr         %26s\n", (base+"_nbr.ens").c_str());
	fprintf(fv.f_case, "vector per node:         1 1 displacement %26s\n", (base+"_dis.ens").c_str());
	fprintf(fv.f_case, "vector per node:         1 1 velocity     %26s\n", (base+"_vel.ens").c_str());
	fprintf(fv.f_case, "vector per node:         1 1 acceleration %26s\n", (base+"_acc.ens").c_str());
	fprintf(fv.f_case, "scalar per node:         1 1 efps         %26s\n", (base+"_eps.ens").c_str());
	fprintf(fv.f_case, "scalar per node:         1 1 density      %26s\n", (base+"_rho.ens").c_str());
	fprintf(fv.f_case, "scalar per node:         1 1 energy       %26s\n", (base+"_egy.ens").c_str());
	fprintf(fv.f_case, "scalar per node:         1 1 soundspeed   %26s\n", (base+"_snd.ens").c_str());
	fprintf(fv.f_case, "scalar per node:         1 1 temperature  %26s\n", (base+"_tmp.ens").c_str());
	fprintf(fv.f_case, "tensor symm per node:    1 1 stress       %26s\n\n", (base+"_sig.ens").c_str());

	fprintf(fv.f_case, "TIME\n");
	fprintf(fv.f_case, "time set:                1\n");
	fprintf(fv.f_case, "number of steps:     %5d\n", outv.sph_istate+1);
	fprintf(fv.f_case, "filename start number:   1\n");
	fprintf(fv.f_case, "filename increment:      1\n");

	if (outv.sph_istate==0)
		fprintf(fv.f_case, "time values:  %12.5e\n", gv.sph_ptime);
	else
	{
		fprintf(fv.f_case, "time values:  %12.5e\n", 0.0);
		for (int i=1; i<(int)outv.sph_timenum.size(); i++)
			fprintf(fv.f_case, "              %12.5e\n", outv.sph_timenum[i]);
		fprintf(fv.f_case, "              %12.5e\n\n", gv.sph_ptime);
	}

	fprintf(fv.f_case, "FILE\n");
	fprintf(fv.f_case, "file set:                1\n");
	fprintf(fv.f_case, "filename index:          1\n");
	fprintf(fv.f_case, "number of steps:     %5d\n", outv.sph_istate+1);

	return OK;
}

int MOutput::WriteCase1D()
{
	MSimulationData &sd = *m_simdata;
	MFileHandlingVars &fv = sd.m_filevars;
	MGlobalVars &gv = sd.m_globvars;
	MOutputVars &outv = sd.m_outvars;
	MParticleData &par = sd.par;
	int i, j, count;

	// Open time case file; outv.sph_istate describes state, new or append
	if (fv.OpenOutputFile("_time.txt", &fv.f_time, outv.sph_istate)!=OK) return ERROR;
	fprintf(fv.f_time, "%13.6e\n", gv.sph_ptime);
	fv.CloseLastOpenedFile();

	// Check state - weather to open new file or to append
	// true - new, false - append
	fv.OpenCaseFiles(outv.sph_istate);

	fprintf(fv.f_sig, "BEGIN TIME STEP\n");
	fprintf(fv.f_dis, "BEGIN TIME STEP\n");
	fprintf(fv.f_vel, "BEGIN TIME STEP\n");
	fprintf(fv.f_acc, "BEGIN TIME STEP\n");
	fprintf(fv.f_eps, "BEGIN TIME STEP\n");
	fprintf(fv.f_rho, "BEGIN TIME STEP\n");
	fprintf(fv.f_egy, "BEGIN TIME STEP\n");
	fprintf(fv.f_snd, "BEGIN TIME STEP\n");
	fprintf(fv.f_mas, "BEGIN TIME STEP\n");
	fprintf(fv.f_bnd, "BEGIN TIME STEP\n");
	fprintf(fv.f_nbr, "BEGIN TIME STEP\n");
	fprintf(fv.f_tmp, "BEGIN TIME STEP\n");

	fprintf(fv.f_sig, "stress\n");
	fprintf(fv.f_dis, "displacement\n");
	fprintf(fv.f_vel, "velocity\n");
	fprintf(fv.f_acc, "acceleration\n");
	fprintf(fv.f_eps, "efps\n");
	fprintf(fv.f_rho, "density\n");
	fprintf(fv.f_egy, "energy\n");
	fprintf(fv.f_snd, "speed of sound\n");
	fprintf(fv.f_mas, "mass\n");
	fprintf(fv.f_bnd, "boundary particles\n");
	fprintf(fv.f_nbr, "number of neighbours\n");
	fprintf(fv.f_tmp, "temperaure\n");

	//int counter=1;
	for (i=1; i<=gv.sph_np; ++i)
	{
		RealVector dis(3);
		dis = par[i].x - par[i].xzero;

		//if (par[i].dispbc != 0) continue;

		fprintf(fv.f_dis, "%12.3e%12.3e%12.3e", dis(0), 0.0, 0.0);
		fprintf(fv.f_vel, "%12.3e%12.3e%12.3e", par[i].v(0), 0.0, 0.0);
		fprintf(fv.f_acc, "%12.3e%12.3e%12.3e", par[i].a(0), 0.0, 0.0);
		fprintf(fv.f_mas, "%12.3e", par[i].mass);
		fprintf(fv.f_bnd, "%12.3e", (real)par[i].boundary);
		fprintf(fv.f_nbr, "%12.3e", (real)par[i].nnbr);

		if (i%2==0)
		{
			fprintf(fv.f_dis, "\n");
			fprintf(fv.f_vel, "\n");
			fprintf(fv.f_acc, "\n");
		}
		if (i%6==0)
		{
			fprintf(fv.f_mas, "\n");
			fprintf(fv.f_bnd, "\n");
			fprintf(fv.f_nbr, "\n");
		}
		//counter++;
	}

	//if ((counter-1)%2 != 0)
	if (gv.sph_np%2 != 0)
	{
		fprintf(fv.f_dis, "\n");
		fprintf(fv.f_vel, "\n");
		fprintf(fv.f_acc, "\n");
	}
	//if ((counter-1)%6 != 0)
	if (gv.sph_np%6 != 0)
	{
		fprintf(fv.f_mas, "\n");
		fprintf(fv.f_bnd, "\n");
		fprintf(fv.f_nbr, "\n");
	}

	// ----------------------------------------------------
	// write the element variables in each file
	// ----------------------------------------------------
	for (j=1; j<=gv.sph_nummat; j++)
	{
		fprintf(fv.f_rho, "part %8d\n", j);
		fprintf(fv.f_sig, "part %8d\n", j);
		fprintf(fv.f_eps, "part %8d\n", j);
		fprintf(fv.f_egy, "part %8d\n", j);
		fprintf(fv.f_snd, "part %8d\n", j);
		fprintf(fv.f_tmp, "part %8d\n", j);

		fprintf(fv.f_rho, "point\n");
		fprintf(fv.f_sig, "point\n");
		fprintf(fv.f_eps, "point\n");
		fprintf(fv.f_egy, "point\n");
		fprintf(fv.f_snd, "point\n");
		fprintf(fv.f_tmp, "point\n");

		for (i=1,count=0; i<=gv.sph_np; i++)
		{
			if (par[i].mat==j)
			{
				count++;
				fprintf(fv.f_sig, "%12.4e%12.4e%12.4e%12.4e%12.4e%12.4e\n",
						par[i].sigma(0,0), 0.0, 0.0, 0.0, 0.0, 0.0 );

				fprintf(fv.f_rho, "%12.4e", par[i].rho);
				fprintf(fv.f_eps, "%12.4e", par[i].efps);
				fprintf(fv.f_egy, "%12.4e", par[i].e / par[i].mass);
				fprintf(fv.f_snd, "%12.4e", par[i].c);
				fprintf(fv.f_tmp, "%12.4e", par[i].temper);

				if (count%6==0)
				{
					fprintf(fv.f_rho, "\n");
					fprintf(fv.f_eps, "\n");
					fprintf(fv.f_egy, "\n");
					fprintf(fv.f_snd, "\n");
					fprintf(fv.f_tmp, "\n");
				}
			}
		} // end particle loop

		if (count%6 != 0)
		{
			fprintf(fv.f_rho, "\n");
			fprintf(fv.f_eps, "\n");
			fprintf(fv.f_egy, "\n");
			fprintf(fv.f_snd, "\n");
			fprintf(fv.f_tmp, "\n");
		}
	} // end material loop

	fprintf(fv.f_sig, "END TIME STEP\n");
	fprintf(fv.f_dis, "END TIME STEP\n");
	fprintf(fv.f_vel, "END TIME STEP\n");
	fprintf(fv.f_acc, "END TIME STEP\n");
	fprintf(fv.f_eps, "END TIME STEP\n");
	fprintf(fv.f_rho, "END TIME STEP\n");
	fprintf(fv.f_egy, "END TIME STEP\n");
	fprintf(fv.f_snd, "END TIME STEP\n");
	fprintf(fv.f_mas, "END TIME STEP\n");
	fprintf(fv.f_bnd, "END TIME STEP\n");
	fprintf(fv.f_nbr, "END TIME STEP\n");
	fprintf(fv.f_tmp, "END TIME STEP\n");

	fv.CloseAllButInputAndLog();

	return OK;
}

int MOutput::WriteCase2D()
{
	MSimulationData &sd = *m_simdata;
	MFileHandlingVars &fv = sd.m_filevars;
	MGlobalVars &gv = sd.m_globvars;
	MOutputVars &outv = sd.m_outvars;
	MParticleData &par = sd.par;
	int i, j, count;

	// Open time case file; outv.sph_istate describes state, new or append
	if (fv.OpenOutputFile("_time.txt", &fv.f_time, outv.sph_istate)!=OK) return ERROR;
	fprintf(fv.f_time, "%13.6e\n", gv.sph_ptime);
	fv.CloseLastOpenedFile();

	// Check state - weather to open new file or to append
	// true - new, false - append
	fv.OpenCaseFiles(outv.sph_istate);

	fprintf(fv.f_sig, "BEGIN TIME STEP\n");
	fprintf(fv.f_dis, "BEGIN TIME STEP\n");
	fprintf(fv.f_vel, "BEGIN TIME STEP\n");
	fprintf(fv.f_acc, "BEGIN TIME STEP\n");
	fprintf(fv.f_eps, "BEGIN TIME STEP\n");
	fprintf(fv.f_rho, "BEGIN TIME STEP\n");
	fprintf(fv.f_egy, "BEGIN TIME STEP\n");
	fprintf(fv.f_snd, "BEGIN TIME STEP\n");
	fprintf(fv.f_mas, "BEGIN TIME STEP\n");
	fprintf(fv.f_bnd, "BEGIN TIME STEP\n");
	fprintf(fv.f_nbr, "BEGIN TIME STEP\n");
	fprintf(fv.f_tmp, "BEGIN TIME STEP\n");

	fprintf(fv.f_sig, "stress\n");
	fprintf(fv.f_dis, "displacement\n");
	fprintf(fv.f_vel, "velocity\n");
	fprintf(fv.f_acc, "acceleration\n");
	fprintf(fv.f_eps, "efps\n");
	fprintf(fv.f_rho, "density\n");
	fprintf(fv.f_egy, "energy\n");
	fprintf(fv.f_snd, "speed of sound\n");
	fprintf(fv.f_mas, "mass\n");
	fprintf(fv.f_bnd, "boundary particles\n");
	fprintf(fv.f_nbr, "number of neighbours\n");
	fprintf(fv.f_tmp, "temperaure\n");

	//int counter=1;
	for (i=1; i<=gv.sph_np; ++i)
	{
		RealVector dis(3);
		dis = par[i].x - par[i].xzero;

		//if (par[i].dispbc != 0) continue;

		fprintf(fv.f_dis, "%12.3e%12.3e%12.3e", dis(0), dis(1), 0.0);
		fprintf(fv.f_vel, "%12.3e%12.3e%12.3e", par[i].v(0), par[i].v(1), 0.0);
		fprintf(fv.f_acc, "%12.3e%12.3e%12.3e", par[i].a(0), par[i].a(1), 0.0);
		fprintf(fv.f_mas, "%12.3e", par[i].mass);
		fprintf(fv.f_bnd, "%12.3e", (real)par[i].boundary);
		fprintf(fv.f_nbr, "%12.3e", (real)par[i].nnbr);

		if (i%2==0)
		{
			fprintf(fv.f_dis, "\n");
			fprintf(fv.f_vel, "\n");
			fprintf(fv.f_acc, "\n");
		}
		if (i%6==0)
		{
			fprintf(fv.f_mas, "\n");
			fprintf(fv.f_bnd, "\n");
			fprintf(fv.f_nbr, "\n");
		}
		//counter++;
	}

	//if ((counter-1)%2 != 0)
	if (gv.sph_np%2 != 0)
	{
		fprintf(fv.f_dis, "\n");
		fprintf(fv.f_vel, "\n");
		fprintf(fv.f_acc, "\n");
	}
	//if ((counter-1)%6 != 0)
	if (gv.sph_np%6 != 0)
	{
		fprintf(fv.f_mas, "\n");
		fprintf(fv.f_bnd, "\n");
		fprintf(fv.f_nbr, "\n");
	}

	// ----------------------------------------------------
	// write the element variables in each file
	// ----------------------------------------------------
	for (j=1; j<=gv.sph_nummat; j++)
	{
		fprintf(fv.f_rho, "part %8d\n", j);
		fprintf(fv.f_sig, "part %8d\n", j);
		fprintf(fv.f_eps, "part %8d\n", j);
		fprintf(fv.f_egy, "part %8d\n", j);
		fprintf(fv.f_snd, "part %8d\n", j);
		fprintf(fv.f_tmp, "part %8d\n", j);

		fprintf(fv.f_rho, "point\n");
		fprintf(fv.f_sig, "point\n");
		fprintf(fv.f_eps, "point\n");
		fprintf(fv.f_egy, "point\n");
		fprintf(fv.f_snd, "point\n");
		fprintf(fv.f_tmp, "point\n");

		for (i=1,count=0; i<=gv.sph_np; i++)
		{
			if (par[i].mat==j)
			{
				count++;
				fprintf(fv.f_sig, "%12.4e%12.4e%12.4e%12.4e%12.4e%12.4e\n",
						par[i].sigma(0,0), par[i].sigma(1,1), 0.0,
						par[i].sigma(0,1), 0.0, 0.0 );

				fprintf(fv.f_rho, "%12.4e", par[i].rho);
				fprintf(fv.f_eps, "%12.4e", par[i].efps);
				fprintf(fv.f_egy, "%12.4e", par[i].e / par[i].mass);
				fprintf(fv.f_snd, "%12.4e", par[i].c);
				fprintf(fv.f_tmp, "%12.4e", par[i].temper);

				if (count%6==0)
				{
					fprintf(fv.f_rho, "\n");
					fprintf(fv.f_eps, "\n");
					fprintf(fv.f_egy, "\n");
					fprintf(fv.f_snd, "\n");
					fprintf(fv.f_tmp, "\n");
				}
			}
		} // end particle loop

		if (count%6 != 0)
		{
			fprintf(fv.f_rho, "\n");
			fprintf(fv.f_eps, "\n");
			fprintf(fv.f_egy, "\n");
			fprintf(fv.f_snd, "\n");
			fprintf(fv.f_tmp, "\n");
		}
	} // end material loop

	fprintf(fv.f_sig, "END TIME STEP\n");
	fprintf(fv.f_dis, "END TIME STEP\n");
	fprintf(fv.f_vel, "END TIME STEP\n");
	fprintf(fv.f_acc, "END TIME STEP\n");
	fprintf(fv.f_eps, "END TIME STEP\n");
	fprintf(fv.f_rho, "END TIME STEP\n");
	fprintf(fv.f_egy, "END TIME STEP\n");
	fprintf(fv.f_snd, "END TIME STEP\n");
	fprintf(fv.f_mas, "END TIME STEP\n");
	fprintf(fv.f_bnd, "END TIME STEP\n");
	fprintf(fv.f_nbr, "END TIME STEP\n");
	fprintf(fv.f_tmp, "END TIME STEP\n");

	fv.CloseAllButInputAndLog();

	return OK;
}

int MOutput::WriteCase3D()
{
	MSimulationData &sd = *m_simdata;
	MFileHandlingVars &fv = sd.m_filevars;
	MGlobalVars &gv = sd.m_globvars;
	MOutputVars &outv = sd.m_outvars;
	MParticleData &par = sd.par;
	int i, j;

	// Open time case file; outv.sph_istate describes state, new or append
	if (fv.OpenOutputFile("_time.txt", &fv.f_time, outv.sph_istate)!=OK) return ERROR;
	fprintf(fv.f_time, "%13.6e\n", gv.sph_ptime);
	fv.CloseLastOpenedFile();

	// Check state - weather to open new file or to append
	// true - new, false - append
	fv.OpenCaseFiles(outv.sph_istate);

	fprintf(fv.f_sig, "BEGIN TIME STEP\n");
	fprintf(fv.f_dis, "BEGIN TIME STEP\n");
	fprintf(fv.f_vel, "BEGIN TIME STEP\n");
	fprintf(fv.f_acc, "BEGIN TIME STEP\n");
	fprintf(fv.f_eps, "BEGIN TIME STEP\n");
	fprintf(fv.f_rho, "BEGIN TIME STEP\n");
	fprintf(fv.f_egy, "BEGIN TIME STEP\n");
	fprintf(fv.f_snd, "BEGIN TIME STEP\n");
	fprintf(fv.f_mas, "BEGIN TIME STEP\n");
	fprintf(fv.f_bnd, "BEGIN TIME STEP\n");
	fprintf(fv.f_nbr, "BEGIN TIME STEP\n");
	fprintf(fv.f_tmp, "BEGIN TIME STEP\n");

	fprintf(fv.f_sig, "stress\n");
	fprintf(fv.f_dis, "displacement\n");
	fprintf(fv.f_vel, "velocity\n");
	fprintf(fv.f_acc, "acceleration\n");
	fprintf(fv.f_eps, "efps\n");
	fprintf(fv.f_rho, "density\n");
	fprintf(fv.f_egy, "energy\n");
	fprintf(fv.f_snd, "speed of sound\n");
	fprintf(fv.f_mas, "mass\n");
	fprintf(fv.f_bnd, "boundary particles\n");
	fprintf(fv.f_nbr, "number of neighbours\n");
	fprintf(fv.f_tmp, "temperaure\n");

	//int counter=1;
	for (i=1; i<=gv.sph_np; ++i)
	{
		RealVector dis(3);
		dis = par[i].x - par[i].xzero;

		//if (par[i].dispbc != 0) continue;

		fprintf(fv.f_dis, "%12.3e%12.3e%12.3e", dis(0), dis(1), dis(2));
		fprintf(fv.f_vel, "%12.3e%12.3e%12.3e", par[i].v(0), par[i].v(1), par[i].v(2));
		fprintf(fv.f_acc, "%12.3e%12.3e%12.3e", par[i].a(0), par[i].a(1), par[i].a(2));
		fprintf(fv.f_mas, "%12.3e", par[i].mass);
		fprintf(fv.f_bnd, "%12.3e", (real)par[i].boundary);
		fprintf(fv.f_nbr, "%12.3e", (real)par[i].nnbr);

		if (i%2==0)
		{
			fprintf(fv.f_dis, "\n");
			fprintf(fv.f_vel, "\n");
			fprintf(fv.f_acc, "\n");
		}
		if (i%6==0)
		{
			fprintf(fv.f_mas, "\n");
			fprintf(fv.f_bnd, "\n");
			fprintf(fv.f_nbr, "\n");
		}
		//counter++;
	}

	//if ((counter-1)%2 != 0)
	if (gv.sph_np%2 != 0)
	{
		fprintf(fv.f_dis, "\n");
		fprintf(fv.f_vel, "\n");
		fprintf(fv.f_acc, "\n");
	}
	//if ((counter-1)%6 != 0)
	if (gv.sph_np%6 != 0)
	{
		fprintf(fv.f_mas, "\n");
		fprintf(fv.f_bnd, "\n");
		fprintf(fv.f_nbr, "\n");
	}

	// ----------------------------------------------------
	// write the element variables in each file
	// ----------------------------------------------------
	for (j=1; j<=gv.sph_nummat; j++)
	{
		fprintf(fv.f_rho, "part %8d\n", j);
		fprintf(fv.f_sig, "part %8d\n", j);
		fprintf(fv.f_eps, "part %8d\n", j);
		fprintf(fv.f_egy, "part %8d\n", j);
		fprintf(fv.f_snd, "part %8d\n", j);
		fprintf(fv.f_tmp, "part %8d\n", j);

		fprintf(fv.f_rho, "point\n");
		fprintf(fv.f_sig, "point\n");
		fprintf(fv.f_eps, "point\n");
		fprintf(fv.f_egy, "point\n");
		fprintf(fv.f_snd, "point\n");
		fprintf(fv.f_tmp, "point\n");

		for (i=1; i<=gv.sph_np; i++)
		{
			if (par[i].mat==j)
			{
				fprintf(fv.f_sig, "%12.4e%12.4e%12.4e%12.4e%12.4e%12.4e\n",
						par[i].sigma(0,0), par[i].sigma(1,1), 0.0,
						par[i].sigma(0,1), 0.0, 0.0 );

				fprintf(fv.f_rho, "%12.4e", par[i].rho);
				fprintf(fv.f_eps, "%12.4e", par[i].efps);
				fprintf(fv.f_egy, "%12.4e", par[i].e / par[i].mass);
				fprintf(fv.f_snd, "%12.4e", par[i].c);
				fprintf(fv.f_tmp, "%12.4e", par[i].temper);

				if (i%6==0)
				{
					fprintf(fv.f_rho, "\n");
					fprintf(fv.f_eps, "\n");
					fprintf(fv.f_egy, "\n");
					fprintf(fv.f_snd, "\n");
					fprintf(fv.f_tmp, "\n");
				}
			}
		} // end particle loop

		if (i%6 != 0)
		{
			fprintf(fv.f_rho, "\n");
			fprintf(fv.f_eps, "\n");
			fprintf(fv.f_egy, "\n");
			fprintf(fv.f_snd, "\n");
			fprintf(fv.f_tmp, "\n");
		}
	} // end material loop

	fprintf(fv.f_sig, "END TIME STEP\n");
	fprintf(fv.f_dis, "END TIME STEP\n");
	fprintf(fv.f_vel, "END TIME STEP\n");
	fprintf(fv.f_acc, "END TIME STEP\n");
	fprintf(fv.f_eps, "END TIME STEP\n");
	fprintf(fv.f_rho, "END TIME STEP\n");
	fprintf(fv.f_egy, "END TIME STEP\n");
	fprintf(fv.f_snd, "END TIME STEP\n");
	fprintf(fv.f_mas, "END TIME STEP\n");
	fprintf(fv.f_bnd, "END TIME STEP\n");
	fprintf(fv.f_nbr, "END TIME STEP\n");
	fprintf(fv.f_tmp, "END TIME STEP\n");

	fv.CloseAllButInputAndLog();

	return OK;
}

int MOutput::WriteCaseGeo()
{
	MSimulationData &sd = *m_simdata;
	MFileHandlingVars &fv = sd.m_filevars;
	MGlobalVars &gv = sd.m_globvars;
	MParticleData &par = sd.par;
	int i, j;

	std::vector<int> par_per_mat;

	// Opens geometry file either for in 'new' or 'append' mode
	fv.OpenOutputFile(".geo", &fv.f_geo, sd.m_outvars.sph_istate);

	// Init particles per material array, one based
	par_per_mat.push_back(0);
	for (j=1; j<=gv.sph_nummat; j++) par_per_mat.push_back(0);

	fprintf(fv.f_geo, "BEGIN TIME STEP\n");
	fprintf(fv.f_geo, "Title1\n");
	fprintf(fv.f_geo, "Title2\n");
	fprintf(fv.f_geo, "node id given\n");
	fprintf(fv.f_geo, "element id given\n");
	fprintf(fv.f_geo, "coordinates\n");
	fprintf(fv.f_geo, "%d\n", gv.sph_np);  //write number of nodes

	// Next loop writes the coordinates of all nodes and counts particles in materials
	//int counter=0;
	for (j=1; j<=gv.sph_np; j++)
	{
		int imaterial = par[j].mat;
		//if (par[j].dispbc==0)
		//{
			//fprintf(fv.f_geo, "%8d%12.3e%12.3e%12.3e\n", j, par[j].x(0)*1e3, par[j].x(1)*1e3, par[j].x(2)*1e3);
			//fprintf(fv.f_geo, "%8d%12.5e%12.5e%12.5e\n", j, par[j].x(0), par[j].x(1), par[j].x(2));
            fprintf(fv.f_geo, "%8d %12.5e %12.5e %12.5e\n", j, par[j].x(0), par[j].x(1), par[j].x(2));
			par_per_mat[imaterial]++;
		//}
	}

	for (j=1; j<=gv.sph_nummat; j++)
	{
		fprintf(fv.f_geo, "part %8d\n", j);
		fprintf(fv.f_geo, "SPH particles in part %8d\n", j);
		fprintf(fv.f_geo, "point\n");
		fprintf(fv.f_geo, "%8d\n", par_per_mat[j]);

		for (i=1; i<=gv.sph_np; i++)
			if (par[i].mat==j) fprintf(fv.f_geo, "%8d%8d\n", i, i);
	}

	fprintf(fv.f_geo, "END TIME STEP\n");

	fv.CloseLastOpenedFile();



	return OK;
}

void MOutput::IncreaseStateCounter()
{
	MSimulationData &sd = *m_simdata;
	MGlobalVars &gv = sd.m_globvars;
	MOutputVars &outv = sd.m_outvars;

	outv.sph_istate++;
	outv.sph_timenum.push_back(gv.sph_ptime);
}

