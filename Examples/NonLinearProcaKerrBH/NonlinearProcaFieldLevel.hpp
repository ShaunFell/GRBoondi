#ifndef NONLINEARPROCAFIELDLEVEL_H_INCLUDED
#define NONLINEARPROCAFIELDLEVEL_H_INCLUDED


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



//Inherits from BaseProcaFieldLevel with background = KerrSchild and matter = ProcaField
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

        //Loop over box cells and assign initial Proca field
        BoxLoops::loop(initial_conditions, m_state_new, m_state_new, INCLUDE_GHOST_CELLS);

        //Excise within horizon
        ExcisionEvolution<ProcaField, KerrSchild> excisor(m_dx, m_p.center, m_p.evolution_excision_width, kerr_schild);

        //Loop over box cells and excise cells within horizon
        BoxLoops::loop(excisor, m_state_new, m_state_new, SKIP_GHOST_CELLS, disable_simd()); //disable SIMD for this loop since excision doesnt use SIMD
        
        //Done!
    }

    // we add the computation of the gnn diagnostic variable here
    virtual void additionalPrePlotLevel() override 
    {
        //Initialize the background class
        KerrSchild kerr_schild(m_p.background_params, m_dx);

        //Initialize the Proca field class
        ProcaField proca_field(kerr_schild, m_p.matter_params);

        //Initialize the constraint class, to be computed on the grid
        SecondClassConstraint<KerrSchild> constraint(m_dx, m_p.center, kerr_schild, proca_field);

        //Loop over the box cells and compute the constraint on the grid
        BoxLoops::loop(constraint, m_state_diagnostics, m_state_diagnostics, SKIP_GHOST_CELLS);
    }



    //Since we use a custom tagging criteria that uses the diagnostic data, 
    //  we need to override the preTagCells method to populate the diagnostic 
    //  data on the grid before tagging
    virtual void preTagCells() override 
    {

        //We only need to 
        if (m_p.activate_gnn_tagging && m_time > m_dt)
        {
            //We dont need any derivatives, so no need to populate ghost cells

            //Initialize the background class
            KerrSchild kerr_schild(m_p.background_params, m_dx);

            //Initialize the Proca field class
            ProcaField proca_field(kerr_schild, m_p.matter_params);

            //Initialize the constraint class, to be computed on the grid
            SecondClassConstraint<KerrSchild> constraint(m_dx, m_p.center, kerr_schild, proca_field);

            //Loop over the box cells and compute the constraint on the grid
            BoxLoops::loop(constraint, m_state_new, m_state_diagnostics, SKIP_GHOST_CELLS);

            const std::vector<int> vars_to_excise= DiagnosticVariables::convert_pairs_to_enum(m_p.diagnostic_excision_vars);
            ExcisionDiagnostics<proca_t,background_t> excisor(background_init, m_dx, m_p.center, m_p.diagnostic_inner_boundary, m_p.diagnostic_outer_boundary,vars_to_excise);

            //excise within the excision zone
            BoxLoops::loop(excisor,m_state_diagnostics, m_state_diagnostics, SKIP_GHOST_CELLS, disable_simd());
        };
    }; //end of preTagCells method



    //compute the tagging criteria on the grid, overriding the Base function
    virtual void computeTaggingCriterion(FArrayBox &tagging_criterion, const FArrayBox &current_state, const FArrayBox &current_state_diagnostics) override
    {
        //We dont need to activate gnn tagging at t=0
        bool activate_gnn_tagging = (m_time > 0.0 && m_p.activate_gnn_tagging) ? true : false;

        //Initialize the tagging class
        CustomTaggingCriterion tagger(m_dx, m_level, m_p.grid_scaling * m_p.L, m_p.initial_ratio, m_p.regrid_threshold, m_p.center, activate_gnn_tagging);

        //Loop over the box cells and compute the tagging criterion
        BoxLoops::loop(tagger, current_state_diagnostics, tagging_criterion);
    }
};


#endif //NONLINEARPROCAFIELDLEVEL_H_INCLUDED