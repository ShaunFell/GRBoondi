#ifndef INITIALPROCADATA_H_INCLUDED
#define INITIALPROCADATA_H_INCLUDED

#include "ADMFixedBGVars.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "GeneralizedProcaField.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp"
#include "VarsTools.hpp"
#include "simd.hpp"
#include "Potential.hpp"
#include "KerrSchild.hpp"


struct init_params_t
    {
        double amplitude;
    };

template <class background_t>
class InitialProcaData: public background_t
{
    public:

        //typename all variables, ADM variables + Matter
        template <class data_t> 
        using MetricVars = ADMFixedBGVars::Vars<data_t>;
        
        //typename only matter variables
        template <class data_t>
        using MatterVars = typename ProcaField<ProcaPotential>::template Vars<data_t>;

        using PotentialParams = ProcaPotential::params_t;
        using KerrParams = typename background_t::params_t;

    protected:
        double m_dx;
        const init_params_t m_params;
        const PotentialParams m_paramsPotential;
        const KerrParams m_paramsKerr;
        const background_t m_background;

    public:
        
        //constructor
        InitialProcaData(init_params_t a_params, PotentialParams b_params, KerrParams c_params, double a_dx, const background_t a_background): 
            background_t(c_params, a_dx), m_dx{a_dx}, m_params{a_params}, m_paramsPotential{b_params}, m_paramsKerr{c_params}, m_background{a_background}
        {
        };

        template <class data_t>
        void compute(Cell<data_t> current_cell) const
        {
            //based off the initial conditions used in    http://arxiv.org/abs/1705.01544

            //location of cell
            Coordinates<data_t> coords(current_cell, m_dx, m_paramsKerr.center);
            
            //compute background variables
            MetricVars<data_t> metric_vars;
            m_background.compute_metric_background(metric_vars, coords);
            
            //flush all variables on cell
            MatterVars<data_t> mattervars;
            VarsTools::assign(mattervars,0.);
            
            const data_t kerrMass = m_paramsKerr.mass;
            const data_t kerrSpin = m_paramsKerr.spin;
            const data_t kerrSpin2 = kerrSpin*kerrSpin;      

            const data_t rho = coords.get_radius(); //x^2 + y^2 + z^2
            const data_t rho2 = rho * rho;

            const data_t coords_z = coords.z;

            //Calculate conformal factor
            data_t chi = TensorAlgebra::compute_determinant_sym(metric_vars.gamma);
            chi = pow(chi, -1.0 / 3.0);


            // the Kerr Schild radius r
            const data_t r2 = 0.5 * (rho2 - kerrSpin2) +
                                sqrt(0.25 * (rho2 - kerrSpin2) * (rho2 - kerrSpin2) + kerrSpin2 * coords_z*coords_z);
            const data_t radius = sqrt(r2);

            //initial profile
            data_t alpha = kerrMass * m_paramsPotential.mass;
            data_t r0_BL { 1.0 / (m_paramsPotential.mass * alpha) };

            //set non-zero grid variables
            mattervars.Avec[0] = m_params.amplitude * pow(chi, 3.) * exp(-radius / r0_BL);

            //export to grid
            current_cell.store_vars(mattervars);
        };
};







#endif //INITIALPROCADATA_H_INCLUDED