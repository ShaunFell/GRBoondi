/* GRBoondi 2024
 * Please refer to LICENSE in GRBoondi's root directory.
 */

/*
   Base level class for the Proca field.
*/

#ifndef BASEPROCAFIELDLEVEL_H_INCLUDED
#define BASEPROCAFIELDLEVEL_H_INCLUDED

// general includes common to most GR problems
#include "AMRReductions.hpp"
#include "BaseProcaFieldLevel.hpp"
#include "BoxLoops.hpp"
#include "ComputePack.hpp"
#include "DefaultLevelFactory.hpp"
#include "GRAMR.hpp"
#include "GRAMRLevel.hpp"
#include "NanCheck.hpp"
#include "SetValue.hpp"
#include "SimulationParameters.hpp"
#include "SmallDataIO.hpp"

// Flux Extraction
#include "ChargesFluxes.hpp"
#include "DataContainer.hpp"
#include "FluxExtraction.hpp"
#include "LinearMomConservation.hpp"

// Diagnostics
#include "DampingField.hpp"
#include "ProcaSquared.hpp"

// RHS Update
#include "MatterEvolution.hpp"

// cell tagging
#include "TaggingCriterion.hpp"

// problem specific includes
#include "BaseProcaField.hpp"
#include "DiagnosticVariables.hpp"
#include "ExcisionEvolution.hpp"

// Excision
#include "ExcisionDiagnostics.hpp"

template <class background_t, class proca_t>
class BaseProcaFieldLevel : public GRAMRLevel
{

    friend class DefaultLevelFactory<
        BaseProcaFieldLevel<background_t, proca_t>>;

    // Define own constructor with non-const SimulationParameters
    // Allows for dynamical modifications of parameters during execution
  public:
    BaseProcaFieldLevel(GRAMR &gr_amr, SimulationParameters &a_p,
                        int a_verbosity)
        : GRAMRLevel(gr_amr, a_p, a_verbosity)
    {
        if (m_verbosity)
            pout() << "BaseProcaField constructor" << endl;
    }

    BaseProcaFieldLevel(){};

    virtual void specificAdvance(); // do things at end of advance step, after
                                    // RK4 calculation

    // pure virtual function. Must be defined by derived class
    virtual void initialData() = 0; // initialize data for the field

#ifdef CH_USE_HDF5
    virtual void prePlotLevel();

    // additional pre-plot routines from the user, e.g. writing custom
    // diagnostic data
    virtual void additionalPrePlotLevel(){};
#endif

    virtual void
    specificEvalRHS(GRLevelData &a_soln, GRLevelData &a_rhs,
                    const double a_time); // RHS routines used at each RK4 step.
                                          // e.g. calculate the time derivatives

    virtual void preTagCells()
        override{}; // things to do before tagging cells (e.g. filling ghosts)

    virtual void specificPostTimeStep() override;

    // additional post-time-step routines from the user, e.g. computating
    // additional fluxes
    virtual void additionalPostTimeStep(){};

    virtual void computeTaggingCriterion(
        FArrayBox &tagging_criterion, const FArrayBox &current_state,
        const FArrayBox &current_state_diagnostics) override;

  public:
    // container classes to store results of flux and density integrations
    // note: Since these are in the level class, each level instance will have
    // its own copy of these
    //              but since the fluxes and integrals are calculated on only a
    //              single level (usually m_level=0), only that level will
    //              contain the actual results
    TimeDataContainer<std::vector<double>>
        m_flux_container{}; // default initialize
    TimeDataContainer<double> m_integral_container{};
};

#include "BaseProcaFieldLevel.impl.hpp"

#endif // BASEPROCAFIELDLEVEL_H_INCLUDED