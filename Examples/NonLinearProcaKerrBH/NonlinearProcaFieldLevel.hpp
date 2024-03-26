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

//constraint tagging
#include "SecondClassConstraint.hpp"
#include "CustomTaggingCriterion.hpp"



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



    //Since we use a custom tagging criteria that uses the diagnostic data, 
    //  we need to override the preTagCells method to populate the diagnostic 
    //  data on the grid before tagging
    virtual void preTagCells() override 
    {

        if (m_p.activate_gnn_tagging)
        {
            //Initialize the background class
            KerrSchild kerr_schild(m_p.background_params, m_dx);

            //Initialize the Proca field class
            ProcaField proca_field(kerr_schild, m_p.matter_params);

            //Initialize the constraint class, to be computed on the grid
            SecondClassConstraint<KerrSchild> constraint(m_dx, m_p.center, kerr_schild, proca_field);

            //Loop over the box cells and compute the constraint on the grid
            BoxLoops::loop(constraint, m_state_new, m_state_diagnostics, SKIP_GHOST_CELLS);

            //Excise within horizon
            ExcisionDiagnostics<ProcaField, KerrSchild> excisor (kerr_schild, m_dx, m_p.center, m_p.evolution_excision_width);

            //Loop over the box cells and excise within the horizon
            BoxLoops::loop(excisor, m_state_new, m_state_diagnostics, SKIP_GHOST_CELLS, disable_simd()); //disable SIMD for this loop since excision doesnt use SIMD
        };
    }; //end of preTagCells method



    //compute the tagging criteria on the grid, overriding the Base function
    virtual void computeTaggingCriterion(FArrayBox &tagging_criterion, const FArrayBox &current_state, const FArrayBox &current_state_diagnostics) override
    {
        //Initialize the tagging class
        CustomTaggingCriterion tagger(m_dx, m_level, m_p.grid_scaling * m_p.L, m_p.initial_ratio, m_p.center, m_p.activate_gnn_tagging);

        //Loop over the box cells and compute the tagging criterion
        BoxLoops::loop(tagger, current_state_diagnostics, tagging_criterion);
    }
};


#endif //PROCAFIELDLEVEL_H_INCLUDED