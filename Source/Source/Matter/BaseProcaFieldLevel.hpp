#ifndef PROCAFIELDLEVEL_H_INCLUDED
#define PROCAFIELDLEVEL_H_INCLUDED

#include "DefaultLevelFactory.hpp"
#include "GRAMRLevel.hpp"
#include "Potential.hpp"
#include "GeneralizedProcaField.hpp"
#include "DefaultLevelAdditions.hpp"

/* #include "GRAMR.hpp"
 */

template <class background_t, class proca_t>
class ProcaFieldLevel : public GRAMRLevel
{

    friend class DefaultLevelFactory<ProcaFieldLevel<background_t, proca_t>>;
    //inherit constructors from GRAMRLevel;
    using GRAMRLevel::GRAMRLevel;
    
    virtual void specificAdvance(); //do things at end of advance step, after RK4 calculation

    //pure virtual function. Must be defined by derived class
    virtual void initialData() = 0; //initialize data for the field

#ifdef CH_USE_HDF5
    virtual void prePlotLevel();
#endif

    virtual void specificEvalRHS(GRLevelData &a_soln, GRLevelData &a_rhs,
                                    const double a_time); //RHS routines used at each RK4 step
    

    virtual void preTagCells() override; //things to do before tagging cells (e.g. filling ghosts)

    virtual void specificPostTimeStep() override;

    virtual void computeTaggingCriterion(
        FArrayBox &tagging_criterion, const FArrayBox &current_state,
        const FArrayBox &current_state_diagnostics) override;

    
};

#endif //PROCAFIELDLEVEL_H_INCLUDED