#ifndef SIMULATIONPARAMETERS_HPP_
#define SIMULATIONPARAMETERS_HPP_

#include "ProcaSimulationParameters.hpp"

class SimulationParameters : public ProcaSimulationParameters
{

  public:
    SimulationParameters(GRParmParse &pp) : ProcaSimulationParameters(pp)
    {
        read_params(pp);
        check_params();
    }

    void read_params(GRParmParse &pp) {}
};

#endif /* SIMULATIONPARAMETERS_HPP_ */