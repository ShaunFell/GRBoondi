/*
Main method to run simulation
*/
#include "ProcaFieldLevel.hpp" // new Level class
#include "SetupFunctions.hpp"  //For setting up MPI processes
#include "runGRBoondi.hpp"     //main run function

int main(int argc, char *argv[])
{
    mainSetup(argc, argv); // setup MPI processes

    int status = runGRBoondi<ProcaFieldLevel>(
        argc, argv); // run simulation with modified level class

    if (status == 0)
    {
        pout() << "GRBoondi finished." << std::endl;
    }
    else
    {
        pout() << "GRBoondi failed with return code " << status << std::endl;
    }

    mainFinalize(); // cleanup MPI processes
    return 0.;
}