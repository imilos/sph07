#include <stdio.h>
#include "global.h"
#include "MSimulationData.h"
#include "MSimulationInit.h"
#include "MSolution.h"
#include "MOutput.h"
#include <mpi.h>

int main(int argc, char **argv)
{
	MSimulationData &simdata = *MSimulationData::Instance();
	MOutput &output = *MOutput::Instance(&simdata);
	MSolution &solution = *MSolution::Instance(&simdata, &output);
	MSimulationInit &init = *MSimulationInit::Instance(&simdata, &solution, &output);

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &simdata.numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &simdata.my_id);

	bool newproblem;

	// Startup procedure
	if ( simdata.Startup(argc, argv, newproblem) != OK )
	{
		simdata.Shutdown(ERROR);
		return ERROR;
	}

	// Get input from file
	if ( simdata.GetInput() != OK )
	{
		simdata.Shutdown(ERROR);
		return ERROR;
	}

	// Initialize everything in order to start simulation
	if ( init.Initialize() != OK )
	{
		simdata.Shutdown(ERROR);
		return ERROR;
	}

	// Simulation itself
	if ( solution.Solution()!=OK )
	{
		simdata.Shutdown(ERROR);
		return ERROR;
	}

	simdata.Shutdown(OK);
	init.DestroyInstance();
	solution.DestroyInstance();
	output.DestroyInstance();
	simdata.DestroyInstance();

	MPI_Finalize();

	return OK;
}
