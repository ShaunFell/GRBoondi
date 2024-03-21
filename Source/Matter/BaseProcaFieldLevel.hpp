#ifndef BASEPROCAFIELDLEVEL_H_INCLUDED
#define BASEPROCAFIELDLEVEL_H_INCLUDED


//general includes common to most GR problems
#include "DefaultLevelFactory.hpp"
#include "GRAMRLevel.hpp"
#include "GRAMR.hpp"
#include "BaseProcaFieldLevel.hpp"
#include "AMRReductions.hpp"
#include "ComputePack.hpp"
#include "BoxLoops.hpp"
#include "NanCheck.hpp"
#include "SetValue.hpp"
#include "SmallDataIO.hpp"

//Flux Extraction
#include "FluxExtraction.hpp"
#include "ChargesFluxes.hpp"

//Diagnostics
#include "ProcaSquared.hpp"
#include "DampingField.hpp"

//RHS Update
#include "MatterEvolution.hpp"

//cell tagging
#include "TaggingCriterion.hpp"

//problem specific includes
#include "DiagnosticVariables.hpp"
#include "ExcisionEvolution.hpp"
#include "BaseProcaField.hpp"

//Excision
#include "ExcisionDiagnostics.hpp"


 

template <class background_t, class proca_t>
class BaseProcaFieldLevel : public GRAMRLevel
{

    friend class DefaultLevelFactory<BaseProcaFieldLevel<background_t, proca_t>>;
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


#include "BaseProcaFieldLevel.impl.hpp"

#endif //BASEPROCAFIELDLEVEL_H_INCLUDED