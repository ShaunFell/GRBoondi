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
        pp.load("num_extraction_radii", extraction_params.num_extraction_radii, 0);
        pp.load("extraction_radii", extraction_params.extraction_radii, extraction_params.num_extraction_radii);
        pp.load("extraction_level", extraction_params.extraction_levels, 1, 0);
        pp.load("num_points_phi", extraction_params.num_points_phi, 25);
        pp.load("num_points_theta", extraction_params.num_points_theta, 37);
        pp.load("extraction_center", extraction_params.center, center);
        pp.load("write_extraction", extraction_params.write_extraction, false);
        pp.load("fluxintegrals_filename", extraction_params.integral_file_prefix, std::string("flux_integrals"));
        std::string extraction_path;
        if (pp.contains("extraction_subpath"))
        {
            pp.load("extraction_subpath", extraction_path);
            if (!extraction_path.empty() && extraction_path.back() != '/')
                extraction_path += "/";
            if (output_path != "./" && !output_path.empty())
                extraction_path = output_path + extraction_path;
        }
        else
            extraction_path = data_path;
        
        extraction_params.data_path = data_path;
        extraction_params.extraction_path = extraction_path;
        pp.load("extraction_file_prefix",
                    extraction_params.extraction_file_prefix,
                    std::string("Weyl4_extraction_"));
        

        pp.load("activate_extraction", activate_extraction, false);
        pp.load("activate_integration", activate_integration, false);

        //grid parameters
        pp.load("nan_check", nan_check, 1);
        pp.load("sigma", sigma, 0.1);

        //excision 
        pp.load("excision_width", excision_width, 1.);

        //tagging
        pp.load("initial_ratio", initial_ratio, 0.25);
        pp.load("activate_ham_tagging", activate_ham_tagging, false);
        pp.load("activate_extraction_tagging", activate_extraction_tagging, false);
        pp.load("grid_scaling", grid_scaling, 1.);
        pp.load("activate_tagging_diagnostic", tagging_diagnostic, false);




    }

    void check_params()
    {
        warn_parameter("inner_r", inner_r, inner_r != 0.0, "set to default parameters (0.0)");
        warn_parameter("outer_r", outer_r, outer_r != 200.0, "set to default parameter (200.0)");
        warn_parameter("activate_tagging_diagnostic", !activate_tagging_diagnostic, "Diagnostic tagging turned on. Tagging criteria will be written to c_Tagging_Diagnostic");

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
    bool activate_extraction_tagging;
    bool activate_extraction;
    double grid_scaling;

    double initial_ratio; 

    bool tagging_diagnostic;



};

#endif /* PROCASIMULATIONPARAMETERS_HPP_ */
