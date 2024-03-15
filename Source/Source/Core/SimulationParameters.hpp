/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SIMULATIONPARAMETERS_HPP_
#define SIMULATIONPARAMETERS_HPP_

// General includes
#include "ChomboParameters.hpp"
#include "GRParmParse.hpp"
/* #include "SimulationParametersBase.hpp" */

// Problem specific includes:
#include "KerrSchild.hpp"
#include "Potential.hpp"
#include "GeneralizedProcaField.hpp"
#include "InitialProcaData.hpp"
#include "SphericalExtraction.hpp"

class SimulationParameters : public ChomboParameters
{
  public:
    SimulationParameters(GRParmParse &pp) : ChomboParameters(pp)
    {
        read_params(pp);
        check_params();
    }

    /// Read parameters from the parameter file
    void read_params(GRParmParse &pp)
    {
        //filenames
        pp.load("integrals_filename", integrals_filename);

        // Initial Kerr data
        pp.load("kerr_mass", kerrSchild_params.mass);
        pp.load("kerr_spin", kerrSchild_params.spin);
        pp.load("kerr_center", kerrSchild_params.center, center);
        //pp.load("kerr_spindir", kerr_params.spin_direction);

        //proca data
        pp.load("proca_mass", potential_params.mass);
        pp.load("proca_self_interaction", potential_params.self_interaction);

        pp.load("proca_damping", proca_params.vector_damping);

        pp.load("initial_proca_amplitude",initialdata_params.amplitude);

        //constants
        pp.load("G_Newton", G_Newton);

        pp.load("activate_integration", activate_integration, false);

        //extraction params
        pp.load("inner_r", inner_r, 0.0);
        pp.load("outer_r", outer_r, 200.0);

        pp.load("activate_extraction", activate_extraction, false);

        //grid parameters
        pp.load("grid_scaling", grid_scaling, 1.);
        pp.load("nan_check", nan_check, 1);
        pp.load("sigma", sigma, 0.1);

        //excision 
        pp.load("excision_width", excision_width, 1.);

    }

    void check_params()
    {
        warn_parameter("kerr_mass", kerrSchild_params.mass, kerrSchild_params.mass >= 0.0,
                       "should be >= 0.0");
 
    
        warn_parameter("inner_r", inner_r, inner_r != 0.0, "set to default parameters (0.0)");
        warn_parameter("outer_r", outer_r, outer_r != 200.0, "set to default parameter (200.0)");

        check_parameter("kerr_spin", kerrSchild_params.spin,
                        std::abs(kerrSchild_params.spin) <= kerrSchild_params.mass,
                        "must satisfy |a| <= M = " +
                            std::to_string(kerrSchild_params.mass));


        check_parameter("excision_with_AH", excise_with_AH, !(excise_with_AH && !AH_activate), "AH Finder must be turned on to use dynamical excision");

        check_parameter("grid_scaling", grid_scaling, grid_scaling>0, "Grid scaling parameter must be greater than zero");

        FOR(idir)
        {
            std::string name = "kerr_center[" + std::to_string(idir) + "]";
            warn_parameter(
                name, kerrSchild_params.center[idir],
                (kerrSchild_params.center[idir] >= 0) &&
                    (kerrSchild_params.center[idir] <= (ivN[idir] + 1) * coarsest_dx),
                "should be within the computational domain");
            
        }
    }

    double G_Newton, outer_r, inner_r;
    double excision_width;
    double sigma;
    int nan_check;
    std::string integrals_filename;


    KerrSchild::params_t kerrSchild_params;
    ProcaPotential::params_t potential_params;
    ProcaField<ProcaPotential>::params_t proca_params;
    spherical_extraction_params_t extraction_params;

    init_params_t  initialdata_params;

    bool activate_integration;
    bool activate_ham_tagging;
    bool activate_gauss_tagging;
    bool activate_extraction;
    double grid_scaling;

    double AH_initial_guess;
    bool excise_with_AH;
    bool AH_activate;


};

#endif /* SIMULATIONPARAMETERS_HPP_ */
