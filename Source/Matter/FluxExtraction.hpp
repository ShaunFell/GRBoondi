/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef FLUXEXTRACTION_HPP_
#define FLUXEXTRACTION_HPP_

#include "SphericalExtraction.hpp"
#include "DiagnosticVariables.hpp"
//!  The class allows extraction of the values of the flux components on
//!  spheroidal shells at specified radii, and integration over those shells
/*!
   The class allows the user to extract data from the grid for the flux
   components over spheroidal shells at specified radii. The values may then be
   written to an output file, or integrated across the surfaces.
*/
class FluxExtraction : public SphericalExtraction
{
    protected:
        std::vector<int> m_vars_to_extract;// vector storing the variables we wish to integrate over the spheres

     public:
        //! The constructor
        FluxExtraction(spherical_extraction_params_t &a_params, std::vector<int> a_vars_to_extract, double a_dt,
                    double a_time, bool a_first_step,
                    double a_restart_time = 0.0)
            : SphericalExtraction(a_params, a_dt, a_time, a_first_step,
                                a_restart_time), m_vars_to_extract(a_vars_to_extract)
        {
            //iterate over variables to extract and add them to the extractor
            for (auto var: m_vars_to_extract)
            {
                add_var(var, VariableType::diagnostic);
            }
        }

        //! The old constructor which assumes it is called in specificPostTimeStep
        //! so the first time step is when m_time == m_dt
        FluxExtraction(spherical_extraction_params_t a_params, std::vector<int> a_vars_to_extract, double a_dt,
                    double a_time, double a_restart_time = 0.0)
            : FluxExtraction(a_params, a_vars_to_extract,a_dt, a_time, (a_dt == a_time),
                            a_restart_time)
        {
        }


        //! Execute the query
        void execute_query(AMRInterpolator<Lagrange<4>> *a_interpolator)
        {
            // extract the values of the Flux scalars on the spheres
            extract(a_interpolator);

            // this would write out the values at every point on the sphere
            if (m_params.write_extraction)
            {
                write_extraction("FluxExtraction");
            }

            // Setup to integrate user specified variables
            std::vector<std::vector<double>> flux_integrals(m_vars_to_extract.size());
            for (int var {0}; var < m_vars_to_extract.size(); var++)
            {
                add_var_integrand(var, flux_integrals[var], IntegrationMethod::simpson);
            }

            // do the integration over the surface
            integrate();

            // Create the header for the integral file
            std::vector<std::string> labels(m_vars_to_extract.size());
            for (int var {0}; var < m_vars_to_extract.size(); var++)
            {
                labels[var] = DiagnosticVariables::variable_names[m_vars_to_extract[var]];
            }

            //write out to file
            write_integrals(m_params.integral_file_prefix, flux_integrals, labels);
        }
};

#endif /* FLUXEXTRACTION_HPP_ */
