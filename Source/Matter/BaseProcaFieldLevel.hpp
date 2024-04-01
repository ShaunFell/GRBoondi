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
#include "LinearMomConservation.hpp"
#include "DataContainer.hpp"

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

    // additional pre-plot routines from the user, e.g. writing custom diagnostic data
    virtual void additionalPrePlotLevel(){};
#endif

    virtual void specificEvalRHS(GRLevelData &a_soln, GRLevelData &a_rhs,
                                    const double a_time); //RHS routines used at each RK4 step
    

    virtual void preTagCells() override {}; //things to do before tagging cells (e.g. filling ghosts)

    virtual void specificPostTimeStep() override;

    // additional post-time-step routines from the user, e.g. computating additional fluxes
    virtual void additionalPostTimeStep() {};

    virtual void computeTaggingCriterion(
        FArrayBox &tagging_criterion, const FArrayBox &current_state,
        const FArrayBox &current_state_diagnostics) override;

    public:
     //container classes to store results of flux and density integrations
     //note: Since these are in the level class, each level instance will have its own copy of these
    //              but since the fluxes and integrals are calculated on only a single level (usually m_level=0), 
    //              only that level will contain the actual results
    DataContainer<std::vector<double>> m_flux_container{}; //default initialize
    DataContainer<double> m_integral_container{};

    
};


#include "BaseProcaFieldLevel.impl.hpp"

#endif //BASEPROCAFIELDLEVEL_H_INCLUDED