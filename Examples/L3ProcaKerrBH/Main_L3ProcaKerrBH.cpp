/*
Main method to run simulation
*/
#include "runGRChombo.hpp" //main run function
#include "L3ProcaFieldLevel.hpp" // new Level class
#include "SetupFunctions.hpp" //For setting up MPI processes



int main(int argc, char *argv[])
{
    mainSetup(argc,argv); //setup MPI processes

    int status = runGRChombo<ProcaFieldLevel>(argc, argv); //run simulation with modified level class

    if (status==0)
    {
        pout() << "GRChombo finished." << std::endl;
    } else {
        pout() << "GRChombo failed with return code " << status << std::endl;
    }

    mainFinalize(); //cleanup MPI processes
    return 0.;
}