#ifndef SIMULATIONPARAMETERS_HPP_
#define SIMULATIONPARAMETERS_HPP_

#include "InitialConditions.hpp"
#include "KerrdeSitter.hpp"
#include "L2_simp.hpp"
#include "ProcaField.hpp"
#include "ProcaSimulationParameters.hpp"

//For interpolating the BH horizon
#include "SimpleDataReader.hpp"
#include "DataManipulation.hpp"

class SimulationParameters : public ProcaSimulationParameters
{

  public:
    SimulationParameters(GRParmParse &pp) : ProcaSimulationParameters(pp)
    {
        pout() << "Reading simulation parameters" << std::endl;
        read_params(pp);
        check_params();

        pout() << "Computing BH horizon" << std::endl;
        compute_bh_horizon();
    }

    void read_params(GRParmParse &pp)
    {
        // Initial Kerr data
        pp.load("kds_mass", background_params.mass);
        pp.load("kds_spin", background_params.spin);
        pp.load("kds_center", background_params.center, center);
        pp.load("kds_CC", background_params.cosmo_constant);

        // Initial EM field params
        pp.load("initial_amplitude", initial_conditions_params.init_amplitude);

        // L2 lagrangian params
        pp.load("L2_alpha", matter_params.alpha2);

        // G2 function params
        pp.load("proca_mass", matter_params.mass);

        // constraint violation damping
        pp.load("z_damping", matter_params.vector_damping);

        // relaxation time of matter field into superradiant growth
        pp.load("relaxation_time", relaxation_time, 0.0);

        // turn on black hole evolution
        pp.load("bh_evolution", evolve_bh, false);
    }

    void check_params()
    {
        const std::vector<int> vars_to_extract =
            DiagnosticVariables::convert_pairs_to_enum(extraction_vars);
        check_parameter("bh_evolution", evolve_bh,
                        DiagnosticVariables::is_variable_to_extract(
                            c_Edot, vars_to_extract) &&
                            DiagnosticVariables::is_variable_to_extract(
                                c_Jdot, vars_to_extract),
                        "If you want to evolve the black hole, you must turn "
                        "on extraction of Edot and Jdot");
    }

    void compute_bh_horizon()
    {
        SimpleDataReader<double> reader {"KerrdeSitter_rPlus.dat"};
        DataContainer<double> data = reader.get_data();

        pout() << "Interpolating BH horizon" << std::endl;

        std::vector<std::vector<double>> coords {data.get_coords()};
        std::vector<double> horizon_values {data.get_data()};

        std::vector<double> query_point { background_params.cosmo_constant, background_params.spin };

        std::pair< std::vector<double>, std::vector<std::vector<double>> > nearest_3_neighbors { DataManipulation::find_nearest_neighbors(coords, query_point, 3) };

        std::vector<double> point1 { nearest_3_neighbors.second[0][0], nearest_3_neighbors.second[0][1] , horizon_values[nearest_3_neighbors.first[0]] };
        std::vector<double> point2 { nearest_3_neighbors.second[1][0], nearest_3_neighbors.second[1][1] , horizon_values[nearest_3_neighbors.first[1]] };
        std::vector<double> point3 { nearest_3_neighbors.second[2][0], nearest_3_neighbors.second[2][1] , horizon_values[nearest_3_neighbors.first[2]] };

        double interpolated_horizon { DataManipulation::lin_interp_2d(point1, point2, point3, query_point) };

        pout() << "Interpolated horizon value: " << interpolated_horizon << std::endl;
        background_params.r_plus = interpolated_horizon;
    }

    // parameters of kerr bh
    KerrdeSitter::params_t background_params;

    // initial conditions parameters
    Initial_Proca_Conditions::params_t initial_conditions_params;

    // Proca parameters
    ProcaField::params_t matter_params;

    double relaxation_time;
    bool evolve_bh;
};

#endif /* SIMULATIONPARAMETERS_HPP_ */