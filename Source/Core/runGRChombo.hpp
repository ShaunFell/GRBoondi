// Chombo includes
#include "parstream.H" //Gives us pout()

// System includes
#include <iostream>

// Our general includes
#include "ProcaLevelFactory.hpp"
#include "AMRInterpolator.hpp"
#include "GRAMR.hpp"
#include "GRParmParse.hpp"
#include "MultiLevelTask.hpp"
#include "SetupFunctions.hpp"

//defaults
#include "DefaultBackground.hpp"
#include "BaseProcaField.hpp"

// Problem specific includes. These must be defined in the problems root directory. E.g. EMKerrBH defines these
#include "BaseProcaFieldLevel.hpp" 
#include "SimulationParameters.hpp"

// Main run function templated over level class

template <class level_t>
int runGRChombo(int argc, char *argv[])
{
    // Load the parameter file and construct the SimulationParameter class
    // To add more parameters edit the SimulationParameters file.
    char *in_file = argv[1];
    GRParmParse pp(argc - 2, argv + 2, NULL, in_file);
    SimulationParameters sim_params(pp);
    pout() << "Symmetry factor of computational grid: " << sim_params.SymmetryFactor << endl;

    
    
    if (sim_params.just_check_params)
        return 0;

    // DefaultLevelFactor is templated over level_t, itself a template parameter
    // Setup the AMR object and initialize the grid
    GRAMR gr_amr;
    ProcaLevelFactory<level_t> problem_level_factory(gr_amr, sim_params);
    setupAMRObject(gr_amr, problem_level_factory);

    //setup interpolating object
    AMRInterpolator<Lagrange<4>> interpolator(
        gr_amr, sim_params.origin, sim_params.dx, sim_params.boundary_params,
        sim_params.verbosity);
    gr_amr.set_interpolator(&interpolator); // also sets puncture_tracker interpolator

    //setup timing routine
    using Clock = std::chrono::steady_clock;
    using Minutes = std::chrono::duration<double, std::ratio<60, 1>>;
    std::chrono::time_point<Clock> start_time = Clock::now();

    // We want to calculate the charges and fluxes at t = 0 
    //call the PostTimeStep right now!!!!
    pout() << "Running initial PostTimeStep" << endl;
    auto task = [](GRAMRLevel *level)
    {
        if (level->time() == 0.)
        {
            level->specificPostTimeStep();
        }
    };

    // call 'now' really now
    MultiLevelTaskPtr<> call_task(task);
    call_task.execute(gr_amr);  


    //go go go !!!!!! Run simulation
    gr_amr.run(sim_params.stop_time, sim_params.max_steps);

    auto now = Clock::now();
    auto duration = std::chrono::duration_cast<Minutes>(now - start_time);
    pout() << "Total simulation time (mins): " << duration.count() << ".\n";

    gr_amr.conclude();

    CH_TIMER_REPORT(); // Report results when running with Chombo timers.

    return 0;
}