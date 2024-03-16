/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

//main run function
#include "runGRChombo.hpp"

//level function
#include "BaseProcaFieldLevel.hpp"



/**
 * Main function for running GRChombo simulation.
 *
 * @param argc number of command line arguments
 * @param argv array of command line argument strings
 *
 * @return status code indicating success or failure of the simulation
 *
 * @throws None
 */
int main(int argc, char *argv[])
{
    //setup MPI processes
    mainSetup(argc, argv);

    //run simulation
    int status = runGRChombo(argc, argv);

    if (status == 0)
        pout() << "GRChombo finished." << std::endl;
    else
        pout() << "GRChombo failed with return code " << status << std::endl;

    //finalize MPI processes
    mainFinalize();

    //done
    return status;
}
