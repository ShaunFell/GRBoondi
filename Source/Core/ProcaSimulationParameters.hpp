/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef PROCASIMULATIONPARAMETERS_HPP_
#define PROCASIMULATIONPARAMETERS_HPP_

// General includes
#include "ChomboParameters.hpp"
#include "GRParmParse.hpp"
#include "UserVariables.hpp"
#include "VariableType.hpp"
#include "BoundaryConditions.hpp"

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
        pp.load("diagnostic_excision_width", diagnostic_excision_width, 0.97);
        pp.load("evolution_excision_width", evolution_excision_width, 0.97);

        //tagging
        pp.load("initial_ratio", initial_ratio, 0.25);
        pp.load("activate_extraction_tagging", activate_extraction_tagging, false);
        pp.load("grid_scaling", grid_scaling, 1.);
        pp.load("regrid_threshold", regrid_threshold, 0.5);

        //linear momentum direction
        pp.load("linear_momentum_dir", linear_momentum_dir, 0);

        //load extraction variables
        UserVariables::load_vars_to_vector(pp, "extraction_vars", "num_extraction_vars", extraction_vars,num_extraction_vars);

        //load integration variables
        UserVariables::load_vars_to_vector(pp, "integration_vars", "num_integration_vars", integration_vars,num_integration_vars);

        //boundary conditions
        boundary_params.read_params(pp);
        //symmetry factor
        FOR(dir)
        {
            if (boundary_params.lo_boundary[dir] == BoundaryConditions::REFLECTIVE_BC || 
            boundary_params.hi_boundary[dir] == BoundaryConditions::REFLECTIVE_BC)
            {
                SymmetryFactor *= 2;
            }
        }
    };

    void check_params()
    {

        check_parameter("grid_scaling", grid_scaling, grid_scaling>0, "Grid scaling parameter must be greater than zero");
        check_parameter("evolution_excision_width", evolution_excision_width, evolution_excision_width <= diagnostic_excision_width, "Evolution excision width must be less than or equal to diagnostic excision width");
    }

    double diagnostic_excision_width;
    double evolution_excision_width;
    double sigma;
    int nan_check;
    std::string integrals_filename;

    spherical_extraction_params_t extraction_params;
    int num_extraction_vars;
    std::vector<std::pair<int, VariableType>> extraction_vars;

    int num_integration_vars;
    std::vector<std::pair<int, VariableType>> integration_vars;

    bool activate_integration;
    bool activate_ham_tagging;
    bool activate_extraction_tagging;
    bool activate_extraction;
    double grid_scaling;

    double initial_ratio; 

    double regrid_threshold;

    //linear momentum direction
    int linear_momentum_dir;


    // Boundary conditions
    BoundaryConditions::params_t boundary_params; // set boundaries in each dir
    double SymmetryFactor {1.0}; // 1.0 for non-reflective boundary conditions, 2.0 for 1 symmetry, 4.0 for 2 symmetries, and 8.0 for 3



};

#endif /* PROCASIMULATIONPARAMETERS_HPP_ */
