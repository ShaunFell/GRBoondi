#ifndef PROCAFIELDLEVEL_H_INCLUDED
#define PROCAFIELDLEVEL_H_INCLUDED

/*
Example of defining a level class that sets initial conditions
*/

#include "BaseProcaFieldLevel.hpp"
#include "BoxLoops.hpp"
#include "ExcisionEvolution.hpp"
#include "InitialConditions.hpp"
#include "KerrSchildNew.hpp"
#include "ProcaField.hpp"
#include "UserVariables.hpp"

// Inherits from BaseProcaFieldLevel with background = KerrSchildNew and matter
// = BaseProcaField
class ProcaFieldLevel : public BaseProcaFieldLevel<KerrSchildNew, ProcaField>
{
  public:
    // inherit constructor from base class
    using BaseProcaFieldLevel::BaseProcaFieldLevel;

    // override method to calculate initial data
    virtual void initialData() override
    {

        KerrSchildNew kerr_schild(m_p.background_params, m_dx);

        // Initialize the initial conditions class
        Initial_Proca_Conditions initial_conditions(
            m_dx, m_p.initial_conditions_params, m_p.matter_params,
            m_p.background_params, kerr_schild);

        // Loop over box cells and assign initial EM field
        BoxLoops::loop(initial_conditions, m_state_new, m_state_new,
                       INCLUDE_GHOST_CELLS);

        // Excise within horizon
        ExcisionEvolution<ProcaField, KerrSchildNew> excisor(m_dx, m_p.center,
                                                             0.95, kerr_schild);

        // Loop over box cells and excise cells within horizon
        BoxLoops::loop(excisor, m_state_new, m_state_new, SKIP_GHOST_CELLS,
                       disable_simd()); // disable SIMD for this loop since
                                        // excision doesnt use SIMD

        // Done!
    }

    // override method to add black hole evolution
    virtual void additionalPostTimeStep() override
    {
        int min_level =
            m_p.extraction_params
                .min_extraction_level(); // get the minimum level for
                                         // extraction, as specified in
                                         // parameter file
        bool at_course_timestep_on_any_level =
            at_level_timestep_multiple(min_level);

        // if the matter field has finished its transient initial relaxation,
        // then we start the black hole evolution Note: The conditions for this
        // code block must match exactly the conditions for the flux calculation
        // code block in BaseProcaFieldLevel.impl.hpp
        if (m_time >= m_p.relaxation_time && m_p.evolve_bh &&
            m_p.activate_extraction && m_level == min_level &&
            at_course_timestep_on_any_level)
        {
            pout() << "Performing BH evolution" << endl;

            // get enums of variables we want to extract
            const std::vector<int> vars_to_extract =
                DiagnosticVariables::convert_pairs_to_enum(m_p.extraction_vars);

            // If we're in this code block, We know c_Edot and c_Jdot exist in
            // the Diagnostic variables because of the check done in
            // ProcaSimulationParameters.hpp
            int Edot_inx = std::find(vars_to_extract.begin(),
                                     vars_to_extract.end(), c_Edot) -
                           vars_to_extract.begin();
            int Jdot_inx = std::find(vars_to_extract.begin(),
                                     vars_to_extract.end(), c_Jdot) -
                           vars_to_extract.begin();

            // Extract the flux values at the horizon
            double Edot_horizon = m_flux_container[m_time][0][Edot_inx];
            double Jdot_horizon = m_flux_container[m_time][0][Jdot_inx];

            // Account for the finite resolution
            double deltaT = m_dt; // from GRAMRLevel
            Edot_horizon *= deltaT;
            Jdot_horizon *= deltaT;

            // Update the black hole mass
            // Note: These flux values are positive if they are outward flowing,
            // and negative if they flow into the horizon
            //               Hence, we subtract the value
            double M = m_p.background_params.mass;
            double a = m_p.background_params.spin;
            double J = a * M;
            m_p.background_params.mass -= Edot_horizon;
            m_p.background_params.spin -=
                a * (Jdot_horizon / J - Edot_horizon / M);

            if (m_verbosity)
            {
                pout() << "Mass update: " << M << " -> "
                       << m_p.background_params.mass << endl;
                pout() << "Spin update: " << a << " -> "
                       << m_p.background_params.spin << endl;
            };

            // Check for naked singularities!
            // Note: This error check is already present in the constructor for
            // the KerrSchildNew class.
            //           We check it here and use a different error message
            //           related to the black hole evolution code above
            if ((m_p.background_params.spin > m_p.background_params.mass) ||
                (m_p.background_params.spin < -m_p.background_params.mass))
            {
                MayDay::Error(" The dimensionless spin has exceeded physical "
                              "bounds. A naked singularity has formed!!");
            }

            bool first_step = (m_time == 0.0);

            // Save new value to file
            SmallDataIO BHParamFile(m_p.data_path + "BHEvolution", m_dt, m_time,
                                    m_restart_time, SmallDataIO::APPEND,
                                    first_step);

            // remove duplicates
            BHParamFile.remove_duplicate_time_data();

            // if its the first time we are saving data, create a header line
            if (first_step)
            {
                std::vector<std::string> header = {"time", "mass", "spin"};
                BHParamFile.write_header_line(header);
            }

            // now write the data
            BHParamFile.write_time_data_line({m_time,
                                              m_p.background_params.mass,
                                              m_p.background_params.spin});
        }
    };
};

#endif // PROCAFIELDLEVEL_H_INCLUDED