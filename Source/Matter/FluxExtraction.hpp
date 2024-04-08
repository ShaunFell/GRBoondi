/* GRBoondi 2024
 */

#ifndef FLUXEXTRACTION_HPP_
#define FLUXEXTRACTION_HPP_

#include "SphericalExtraction.hpp"
#include "DiagnosticVariables.hpp"
#include "DataContainer.hpp"
//!  The class allows extraction of the values of the flux components on
//!  spheroidal shells at specified radii, and integration over those shells
/*!
   The class allows the user to extract data from the grid for the flux
   components over spheroidal shells at specified radii. The values may then be
   written to an output file, or integrated across the surfaces.
*/
class FluxExtraction : public SphericalExtraction
{
    using Container = DataContainer<std::vector<double>>;
    protected:
        std::vector<int> m_vars_to_extract;// vector storing the variables we wish to integrate over the spheres
        Container& m_flux_container; //container 
        double m_symmetry_mult; //symmetry multiplier for integration

     public:

        // Constructor that allows passing a container object by reference
        FluxExtraction(Container& a_flux_container, spherical_extraction_params_t &a_params, std::vector<int> a_vars_to_extract, double a_dt,
                    double a_time, bool a_first_step,
                    double a_restart_time = 0.0, double a_symmetry_mult = 1.0)
            : SphericalExtraction(a_params, a_dt, a_time, a_first_step,
                                a_restart_time), m_vars_to_extract(a_vars_to_extract), m_flux_container(a_flux_container), m_symmetry_mult(a_symmetry_mult)
        {
            //iterate over variables to extract and add them to the extractor
            for (auto var: m_vars_to_extract)
            {
                add_var(var, VariableType::diagnostic);
            }
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
                bool broadcast_result { true }; //broadcast the result to the other MPI processes
                add_var_integrand(var, flux_integrals[var], IntegrationMethod::simpson, IntegrationMethod::simpson, broadcast_result);
            }

            // do the integration over the surface
            integrate();

            // Create the header for the integral file
            std::vector<std::string> labels(m_vars_to_extract.size());
            for (int var {0}; var < m_vars_to_extract.size(); var++)
            {
                labels[var] = DiagnosticVariables::variable_names[m_vars_to_extract[var]];
            }

            //add symmetry multiplier
            for (std::vector<double>& vector_results: flux_integrals)
            {
                for (double& integral_result: vector_results)
                {
                    integral_result *= m_symmetry_mult;
                }
            }

            //write out to file
            write_integrals(m_params.integral_file_prefix, flux_integrals, labels);

            //write data to container object
            m_flux_container.update(m_time, flux_integrals);
        }
};

#endif /* FLUXEXTRACTION_HPP_ */
