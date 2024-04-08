/* GRBoondi 2024
 */

#ifndef PROCALEVELFACTORY_H_INCLUDED
#define PROCALEVELFACTORY_H_INCLUDED

// Chombo includes:
#include "AMRLevelFactory.H"

// Our includes
#include "GRAMR.hpp"
#include "SimulationParameters.hpp"

// Chombo namespace
#include "UsingNamespace.H"

template <class level_t> class ProcaLevelFactory : public AMRLevelFactory
{
  public:
    ProcaLevelFactory(GRAMR &gr_amr, SimulationParameters &a_sim_params)
        : m_gr_amr(gr_amr), m_p(a_sim_params)
    {
    }

    virtual AMRLevel *new_amrlevel() const
    {
        level_t *level_ptr = new level_t(m_gr_amr, m_p, m_p.verbosity);
        level_ptr->initialDtMultiplier(m_p.dt_multiplier);
        return (static_cast<AMRLevel *>(level_ptr));
    }

    virtual ~ProcaLevelFactory() {}

  protected:
    GRAMR &m_gr_amr;
    SimulationParameters &m_p;
};
#endif /* PROCALEVELFACTORY_H_INCLUDED */