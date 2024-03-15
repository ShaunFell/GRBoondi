/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

// Chombo includes
#include "parstream.H" //Gives us pout()

// System includes
#include <iostream>

// Our general includes
#include "DefaultLevelFactory.hpp"
#include "AMRInterpolator.hpp"
#include "GRAMR.hpp"
#include "GRParmParse.hpp"
#include "MultiLevelTask.hpp"
#include "SetupFunctions.hpp"
#include "SimulationParameters.hpp"

// Problem specific includes:
#include "BaseProcaFieldLevel.hpp"


template <class level_t>
int runGRChombo(int argc, char *argv[])
{
    // Load the parameter file and construct the SimulationParameter class
    // To add more parameters edit the SimulationParameters file.
    char *in_file = argv[1];
    GRParmParse pp(argc - 2, argv + 2, NULL, in_file);
    SimulationParameters sim_params(pp);

    if (sim_params.just_check_params)
        return 0;

    // The line below selects the problem that is simulated
    // (To simulate a different problem, define a new child of AMRLevel
    // and an associated LevelFactory)
    GRAMR bh_amr;
    DefaultLevelFactory<level_t> proca_field_level_fact(bh_amr, sim_params);
    setupAMRObject(bh_amr, proca_field_level_fact);

    //setup interpolating object
    AMRInterpolator<Lagrange<4>> interpolator(
        bh_amr, sim_params.origin, sim_params.dx, sim_params.boundary_params,
        sim_params.verbosity);
    bh_amr.set_interpolator(&interpolator); // also sets puncture_tracker interpolator

    using Clock = std::chrono::steady_clock;
    using Minutes = std::chrono::duration<double, std::ratio<60, 1>>;

    std::chrono::time_point<Clock> start_time = Clock::now();

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
    call_task.execute(bh_amr); 


    //go go go !!!!!! Run simulation
    bh_amr.run(sim_params.stop_time, sim_params.max_steps);

    auto now = Clock::now();
    auto duration = std::chrono::duration_cast<Minutes>(now - start_time);
    pout() << "Total simulation time (mins): " << duration.count() << ".\n";

    bh_amr.conclude();

    CH_TIMER_REPORT(); // Report results when running with Chombo timers.

    return 0;
}