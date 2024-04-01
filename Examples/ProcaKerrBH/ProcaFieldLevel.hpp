#ifndef PROCAFIELDLEVEL_H_INCLUDED
#define PROCAFIELDLEVEL_H_INCLUDED

/*
Example of defining a level class that sets initial conditions
*/

#include "BaseProcaFieldLevel.hpp"
#include "ProcaField.hpp"
#include "UserVariables.hpp"
#include "InitialConditions.hpp"
#include "ExcisionEvolution.hpp"
#include "KerrSchild.hpp"
#include "BoxLoops.hpp"



//Inherits from BaseProcaFieldLevel with background = KerrSchild and matter = BaseProcaField
class ProcaFieldLevel : public BaseProcaFieldLevel<KerrSchild, ProcaField>
{
public:
    //inherit constructor from base class (note: BaseProcaFieldLevel itself inherits from GRAMRLevel)
    using BaseProcaFieldLevel::BaseProcaFieldLevel;


    //override method to calculate initial data
    virtual void initialData() override {
        
        KerrSchild kerr_schild(m_p.background_params, m_dx);

        //Initialize the initial conditions class
        Initial_Proca_Conditions initial_conditions(m_dx, m_p.initial_conditions_params, m_p.matter_params, m_p.background_params, kerr_schild);

        //Loop over box cells and assign initial EM field
        BoxLoops::loop(initial_conditions, m_state_new, m_state_new, INCLUDE_GHOST_CELLS);

        //Excise within horizon
        ExcisionEvolution<ProcaField, KerrSchild> excisor(m_dx, m_p.center, m_p.evolution_excision_width, kerr_schild);

        //Loop over box cells and excise cells within horizon
        BoxLoops::loop(excisor, m_state_new, m_state_new, SKIP_GHOST_CELLS, disable_simd()); //disable SIMD for this loop since excision doesnt use SIMD
        
        //Done!
    }

    //override method to add black hole evolution
    virtual void additionalPostTimeStep() override
    {
        //if the matter field has finished its transient initial relaxation, then we start the black hole evolution
        if (m_time > m_p.relaxation_time && m_p.bh_evolution)
        {

        }

    };
};


#endif //PROCAFIELDLEVEL_H_INCLUDED