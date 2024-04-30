#ifndef PROCAFIELDLEVEL_H_INCLUDED
#define PROCAFIELDLEVEL_H_INCLUDED

/*
Example of defining a level class that sets initial conditions
*/

#include "BaseProcaFieldLevel.hpp"
#include "BoxLoops.hpp"
#include "EnhancedExcisionEvolution.hpp"
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
        KerrdeSitter background_init{m_p.background_params, m_dx};
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

    /*     //override method to add black hole evolution
        virtual void additionalPostTimeStep() override
        {
            int min_level = m_p.extraction_params.min_extraction_level(); //get
       the minimum level for extraction, as specified in parameter file bool
       at_course_timestep_on_any_level = at_level_timestep_multiple(min_level);

            //if the matter field has finished its transient initial relaxation,
       then we start the black hole evolution
            //Note: The conditions for this code block must match exactly the
       conditions for the flux calculation code block in
       BaseProcaFieldLevel.impl.hpp if (m_time >= m_p.relaxation_time &&
       m_p.evolve_bh &&m_p.activate_extraction && m_level == min_level &&
       at_course_timestep_on_any_level)
            {
                pout() << "Performing BH evolution" << endl;

                //get enums of variables we want to extract
                const std::vector<int> vars_to_extract =
       DiagnosticVariables::convert_pairs_to_enum(m_p.extraction_vars);

                // If we're in this code block, We know c_Edot and c_Jdot exist
       in the Diagnostic variables because of the check done in
       ProcaSimulationParameters.hpp int Edot_inx =
       std::find(vars_to_extract.begin(), vars_to_extract.end(), c_Edot) -
       vars_to_extract.begin(); int Jdot_inx =
       std::find(vars_to_extract.begin(), vars_to_extract.end(), c_Jdot) -
       vars_to_extract.begin();

                //Extract the flux values at the horizon
                double Edot_horizon = m_flux_container[m_time][0][Edot_inx];
                double Jdot_horizon = m_flux_container[m_time][0][Jdot_inx];

                //Account for the finite resolution
                double deltaT = m_dt; //from GRAMRLevel
                Edot_horizon *= deltaT;
                Jdot_horizon *= deltaT;

                //Update the black hole mass
                //Note: These flux values are positive if they are outward
       flowing, and negative if they flow into the horizon
                //              Hence, we subtract the value
                double M = m_p.background_params.mass;
                double a = m_p.background_params.spin;
                double J = a * M;
                m_p.background_params.mass -= Edot_horizon;
                m_p.background_params.spin -= a * ( Jdot_horizon / J -
       Edot_horizon / M );

                if (m_verbosity)
                {
                pout() << "Mass update: " << M << " -> " <<
       m_p.background_params.mass << endl; pout() << "Spin update: " << a << "
       -> " << m_p.background_params.spin << endl;
                };

                //Check for naked singularities!
                //Note: This error check is already present in the constructor
       for the KerrdeSitter class.
                //          We check it here and use a different error message
       related to the black hole evolution code above if
       ((m_p.background_params.spin > m_p.background_params.mass) ||
       (m_p.background_params.spin < -m_p.background_params.mass))
                {
                    MayDay::Error(
                            " The dimensionless spin has exceeded physical
       bounds. A naked singularity has formed!!"
                    );
                }

            }

        };*/
};

#endif // PROCAFIELDLEVEL_H_INCLUDED