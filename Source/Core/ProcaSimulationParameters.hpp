/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef PROCASIMULATIONPARAMETERS_HPP_
#define PROCASIMULATIONPARAMETERS_HPP_

// General includes
#include "ChomboParameters.hpp"
#include "GRParmParse.hpp"
/* #include "SimulationParametersBase.hpp" */

// Proca specific includes:
#include "BaseProcaField.hpp"
#include "SphericalExtraction.hpp"

class ProcaSimulationParameters : public ChomboParameters
{
  public:
    ProcaSimulationParameters(GRParmParse &pp) : ChomboParameters(pp)
    {
        read_params(pp);
        check_params();
    }

    /// Read parameters from the parameter file
    void read_params(GRParmParse &pp)
    {
        //filenames
        pp.load("integrals_filename", integrals_filename);


        //extraction params
        pp.load("activate_integration", activate_integration, false);
        pp.load("inner_r", inner_r, 0.0);
        pp.load("outer_r", outer_r, 200.0);

        pp.load("activate_extraction", activate_extraction, false);

        //grid parameters
        pp.load("nan_check", nan_check, 1);
        pp.load("sigma", sigma, 0.1);

        //excision 
        pp.load("excision_width", excision_width, 1.);

        //tagging
        pp.load("initial_ratio", initial_ratio, 0.25);
        pp.load("activate_ham_tagging", activate_ham_tagging, false);
        pp.load("grid_scaling", grid_scaling, 1.);


    }

    void check_params()
    {
        warn_parameter("inner_r", inner_r, inner_r != 0.0, "set to default parameters (0.0)");
        warn_parameter("outer_r", outer_r, outer_r != 200.0, "set to default parameter (200.0)");

        check_parameter("grid_scaling", grid_scaling, grid_scaling>0, "Grid scaling parameter must be greater than zero");
            
    }

    double outer_r, inner_r;
    double excision_width;
    double sigma;
    int nan_check;
    std::string integrals_filename;

    spherical_extraction_params_t extraction_params;

    bool activate_integration;
    bool activate_ham_tagging;
    bool activate_extraction;
    double grid_scaling;

    double initial_ratio; 



};

#endif /* PROCASIMULATIONPARAMETERS_HPP_ */
