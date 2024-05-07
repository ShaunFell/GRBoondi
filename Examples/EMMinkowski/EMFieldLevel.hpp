#ifndef EMFIELDLEVEL_HPP_
#define EMFIELDLEVEL_HPP_

/*
Example of defining a level class that sets simple initial conditions
*/

#include "BaseProcaFieldLevel.hpp"
#include "BoxLoops.hpp"
#include "DefaultBackground.hpp"
#include "ExcisionEvolution.hpp"
#include "InitialConditions.hpp"
#include "ProcaField.hpp"
#include "UserVariables.hpp"

// Inherits from BaseProcaFieldLevel with background = DefaultBackground and
// matter = Default
class EMFieldLevel : public BaseProcaFieldLevel<DefaultBackground, ProcaField>
{
  public:
    // inherit constructor from base class
    using BaseProcaFieldLevel::BaseProcaFieldLevel;

    // override method to calculate initial data
    virtual void initialData() override
    {

        // Initialize the initial conditions class
        Initial_EM_Conditions initial_conditions(
            m_dx, m_p.initial_conditions_params, m_p.background_params);

        // Loop over box cells and assign initial EM field
        BoxLoops::loop(initial_conditions, m_state_new, m_state_new,
                       INCLUDE_GHOST_CELLS);

        // Excise within horizon
        DefaultBackground kerr_schild(m_p.background_params, m_dx);
        ExcisionEvolution<BaseProcaField<DefaultBackground, ProcaField>,
                          DefaultBackground>
            excisor(m_dx, m_p.center, m_p.evolution_excision_width,
                    kerr_schild);

        // Loop over box cells and excise cells within horizon
        BoxLoops::loop(excisor, m_state_new, m_state_new, SKIP_GHOST_CELLS,
                       disable_simd()); // disable SIMD for this loop since
                                        // excision doesnt use SIMD

        // Done!
    }
};

#endif // EMFIELDLEVEL_HPP_