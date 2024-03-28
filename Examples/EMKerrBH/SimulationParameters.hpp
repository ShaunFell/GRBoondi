#ifndef SIMULATIONPARAMETERS_HPP_
#define SIMULATIONPARAMETERS_HPP_

#include "ProcaSimulationParameters.hpp"
#include "KerrSchild.hpp"
#include "InitialConditions.hpp"
#include "ProcaField.hpp"


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

            //Initial EM field params
            pp.load("initial_em_field", initial_conditions_params.init_amplitude);
        }

        //parameters of kerr bh
        KerrSchild::params_t background_params;
        
        //initial conditions parameters
        Initial_EM_Conditions::params_t initial_conditions_params;

        //Proca field params. In this class, an empty struct
        ProcaField::params_t matter_params;


};





#endif /* SIMULATIONPARAMETERS_HPP_ */