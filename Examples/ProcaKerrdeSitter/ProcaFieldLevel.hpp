#ifndef PROCAFIELDLEVEL_H_INCLUDED
#define PROCAFIELDLEVEL_H_INCLUDED

/*
Example of defining a level class that sets initial conditions
*/

#include "BaseProcaFieldLevel.hpp"
#include "BoxLoops.hpp"
#include "DumbExcisionEvolution.hpp"
#include "InitialConditions.hpp"
#include "KerrdeSitter.hpp"
#include "ProcaField.hpp"
#include "UserVariables.hpp"

// Inherits from BaseProcaFieldLevel with background = KerrdeSitter and matter =
// BaseProcaField
class ProcaFieldLevel : public BaseProcaFieldLevel<KerrdeSitter, ProcaField>
{
  public:
    // inherit constructor from base class
    using BaseProcaFieldLevel::BaseProcaFieldLevel;

    // override method to calculate initial data
    virtual void initialData() override
    {

        KerrdeSitter kerr_schild(m_p.background_params, m_dx);

        // Initialize the initial conditions class
        Initial_Proca_Conditions initial_conditions(
            m_dx, m_p.initial_conditions_params, m_p.matter_params,
            m_p.background_params, kerr_schild);

        // Loop over box cells and assign initial EM field
        BoxLoops::loop(initial_conditions, m_state_new, m_state_new,
                       INCLUDE_GHOST_CELLS);

        // Excise within horizon
        DumbExcisionEvolution excisor(
            m_dx, m_p.center, m_p.evolution_excision_inner,
            m_p.evolution_excision_outer, m_p.background_params.r_plus);

        // Loop over box cells and excise cells within horizon
        BoxLoops::loop(excisor, m_state_new, m_state_new, SKIP_GHOST_CELLS,
                       disable_simd()); // disable SIMD for this loop since
                                        // excision doesnt use SIMD

        // Done!
    }

    // override default method for new excision routine
    virtual void specificEvalRHS(GRLevelData &a_soln, GRLevelData &a_rhs,
                                 const double a_time) override
    {
        // Calculate right hand side
        KerrdeSitter background_init(m_p.background_params, m_dx);
        ProcaField proca_field(background_init, m_p.matter_params);

        MatterEvolution<ProcaField, KerrdeSitter> matter_class(
            proca_field, background_init, m_p.sigma, m_dx, m_p.center);
        DumbExcisionEvolution excisor(
            m_dx, m_p.center, m_p.evolution_excision_inner,
            m_p.evolution_excision_outer, m_p.background_params.r_plus);

        BoxLoops::loop(matter_class, a_soln, a_rhs, SKIP_GHOST_CELLS);
        BoxLoops::loop(excisor, a_soln, a_rhs, SKIP_GHOST_CELLS,
                       disable_simd());
    }
};

#endif // PROCAFIELDLEVEL_H_INCLUDED