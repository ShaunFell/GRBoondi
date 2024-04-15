#ifndef SIMULATIONPARAMETERS_HPP_
#define SIMULATIONPARAMETERS_HPP_

#include "ProcaSimulationParameters.hpp"
#include "KerrSchild.hpp"
#include "L2_simp.hpp"
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
        }

        void check_params()
        {
        }

};





#endif /* SIMULATIONPARAMETERS_HPP_ */