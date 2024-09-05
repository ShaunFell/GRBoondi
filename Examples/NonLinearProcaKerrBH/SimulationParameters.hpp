#ifndef SIMULATIONPARAMETERS_HPP_
#define SIMULATIONPARAMETERS_HPP_

#include "InitialConditions.hpp"
#include "KerrSchildNew.hpp"
#include "L2_simp.hpp"
#include "ProcaField.hpp"
#include "ProcaSimulationParameters.hpp"

class SimulationParameters : public ProcaSimulationParameters
{

  public:
    SimulationParameters(GRParmParse &pp) : ProcaSimulationParameters(pp)
    {
        read_params(pp);
        check_params();
    }

    void read_params(GRParmParse &pp)
    {
        // Initial Kerr data
        pp.load("kerr_mass", background_params.mass);
        pp.load("kerr_spin", background_params.spin);
        pp.load("kerr_center", background_params.center, center);

        // Initial EM field params
        pp.load("initial_amplitude", initial_conditions_params.init_amplitude);

        // L2 lagrangian params
        pp.load("L2_alpha", matter_params.alpha2);

        // G2 function params
        pp.load("proca_mass", matter_params.mass);

        // constraint violation damping
        pp.load("z_damping", matter_params.vector_damping);

        // self interaction coupling
        pp.load("self_interaction", matter_params.self_interaction);

        // tagging criteria
        pp.load("activate_gnn_tagging", activate_gnn_tagging, false);

        // relaxation time of matter field into superradiant growth
        pp.load("relaxation_time", relaxation_time, 0.0);

        // turn on black hole evolution
        pp.load("bh_evolution", evolve_bh, false);
    }

    // parameters of kerr bh
    KerrSchildNew::params_t background_params;

    // initial conditions parameters
    Initial_Proca_Conditions::params_t initial_conditions_params;

    // Proca parameters
    ProcaField::params_t matter_params;

    // tagging criteria
    bool activate_gnn_tagging;

    double relaxation_time;
    bool evolve_bh;
};

#endif /* SIMULATIONPARAMETERS_HPP_ */